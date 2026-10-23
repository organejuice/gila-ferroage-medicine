# Local, additive correction of SOP035 matching. Frozen source artifacts remain read-only.
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(Matrix);library(yaml);library(jsonlite);library(digest)})
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)>=2)
cfg<-read_yaml(args[1]);mode<-args[2];ncfg<-read_yaml(cfg$network_config);rcfg<-read_yaml(cfg$retention_config)
out<-cfg$output_dir;interim<-cfg$interim_dir;dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(interim,recursive=TRUE,showWarnings=FALSE)
source('src/130_virtual_ko_cached_backend.R');source('src/matching_correction_036_helpers.R')
tsv<-function(x,p)write.table(x,p,sep='\t',row.names=FALSE,quote=FALSE,na='')
atomic_rds<-function(x,p){saveRDS(x,paste0(p,'.partial'));stopifnot(file.rename(paste0(p,'.partial'),p))}
sha<-function(p){con<-file(p,'rb');on.exit(close(con));digest(readBin(con,'raw',n=file.info(p)$size),algo='sha256',serialize=FALSE)}
scope<-read.delim(cfg$scope_path,check.names=FALSE,encoding='UTF-8');names(scope)[1]<-'analysis'
excluded<-unique(c(ncfg$targets_microglia,ncfg$targets_astrocyte))
origin<-function(a)if(a=='primary')'results/2026-10-04-stratified-virtual-ko-033'else'results/2026-10-05-virtual-ko-target-retention-sensitivity-033'
baseline<-function(a,c)readRDS(file.path(cfg$baseline_dir,paste0(a,'__',c,'__baseline.rds')))
ko<-function(wt,target,seed,path,network_hash){
 sig<-digest(list(network_sha256=network_hash,target=target,seed=seed,config=sha(args[1]),backend=sha('src/130_virtual_ko_cached_backend.R'),adapter=sha('src/correct_virtual_ko_matching_036.R')),algo='sha256')
 marker<-paste0(path,'.json')
 if(file.exists(marker)){m<-fromJSON(marker);stopifnot(m$signature==sig);if(m$status=='COMPLETED')stopifnot(sha(path)==m$sha256);return(m$status)}
 ans<-delete_knk_gene(wt,target,ncfg,seed,FALSE)
 if(ans$status=='COMPLETED'){d<-ans$diffRegulation;stopifnot(!anyDuplicated(d$gene),all(is.finite(d$distance)),all(d$distance>=0),all(is.finite(d$p.adj)))
  tsv(d,paste0(path,'.partial'));stopifnot(file.rename(paste0(path,'.partial'),path))}
 write_json(list(status=ans$status,signature=sig,sha256=if(ans$status=='COMPLETED')sha(path)else NA),marker,auto_unbox=TRUE,pretty=TRUE)
 ans$status
}
if(mode=='prepare'){
 capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'));summary<-list()
 bg<-unique(scope[,c('analysis','context')])
 for(i in seq_len(nrow(bg))){a<-bg$analysis[i];c<-bg$context[i];name<-paste(a,c,sep='__');cat('PREPARE',name,'\n');flush.console()
  b<-baseline(a,c);ef<-full_expression_features(b);f<-ef$features;d<-file.path(out,'corrected',name);dir.create(d,recursive=TRUE,showWarnings=FALSE)
  tsv(f,file.path(d,'expression_features.tsv'));tsv(ef$libraries,file.path(d,'library_denominator_audit.tsv'))
  features<-list()
  for(seed in cfg$seeds){g<-readRDS(file.path(origin(a),c,paste0('seed_',seed),'WT_network.rds'))
   stopifnot(identical(rownames(g$WT),rownames(b$counts)));ff<-f[match(rownames(g$WT),f$gene),]
   ff$out_degree<-rowSums(g$WT!=0);ff$out_abs_weight_sum<-rowSums(abs(g$WT));features[[as.character(seed)]]<-ff
   tsv(ff,file.path(d,paste0('WT_features_seed_',seed,'.tsv')));rm(g);gc()
  }
  pool<-Reduce(intersect,lapply(features,function(q)q$gene[q$out_degree>0&!q$gene%in%excluded]))
  tt<-scope$target[scope$analysis==a&scope$context==c];sel<-select_references(features,tt,pool,cfg$matched_controls)
  tsv(sel$selected,file.path(d,'matched_target_controls.tsv'));tsv(sel$ranking,file.path(d,'complete_reference_ranking.tsv'))
  old<-read.delim(file.path(cfg$source_dir,'matched',name,'matched_target_controls.tsv'))
  for(t in tt){oldc<-unique(old$control[old$target==t]);newc<-unique(sel$selected$control[sel$selected$target==t]);ii<-match(t,f$gene)
   summary[[length(summary)+1L]]<-data.frame(analysis=a,context=c,target=t,old_CP10K=f$old_mean_CP10K[ii],corrected_CP10K=f$mean_CP10K[ii],
    old_controls=length(oldc),new_controls=length(newc),retained_controls=length(intersect(oldc,newc)),added_controls=paste(setdiff(newc,oldc),collapse=';'),removed_controls=paste(setdiff(oldc,newc),collapse=';'))
  }
  if(a==cfg$expanded_analysis&&c==cfg$expanded_context){
   low<-f[f$detected_cells>=rcfg$target_retention_min_cells&f$detection_fraction<ncfg$gene_min_fraction&!f$gene%in%excluded,]
   stopifnot(nrow(low)>0);pre<-list()
   for(t in cfg$expanded_targets){q<-rbind(f[f$gene==t,],low);x<-apply(as.matrix(q[,c('detection_fraction','mean_CP10K')]),2,function(v)rank(v,ties.method='average')/(length(v)+1))
    dist<-sqrt(rowMeans(sweep(x,2,x[1,],'-')^2));q<-q[-1,];q$expression_distance<-dist[-1];q<-q[order(q$expression_distance,q$gene),]
    q$preselected<-seq_len(nrow(q))<=cfg$preselection_per_target;q$target<-t;pre[[t]]<-q}
   pre<-do.call(rbind,pre);tsv(pre,file.path(out,'low_expression_preselection_ranking.tsv'))
   refs<-sort(unique(pre$gene[pre$preselected]));tsv(data.frame(gene=refs),file.path(out,'frozen_expansion_reference_genes.tsv'))
   genes<-sort(unique(c(rownames(b$counts),refs)));z<-b$all_gene_counts[genes,,drop=FALSE]
   atomic_rds(list(counts=z,metadata=b$metadata,original_genes=rownames(b$counts),reference_genes=refs),file.path(interim,'expanded_input.rds'))
   write_json(list(status='FROZEN_BEFORE_NETWORK_AND_KO',reference_genes=length(refs),original_genes=nrow(b$counts),expanded_genes=nrow(z),
    input_sha256=sha(file.path(interim,'expanded_input.rds')),reference_sha256=sha(file.path(out,'frozen_expansion_reference_genes.tsv'))),file.path(out,'expansion_freeze.json'),auto_unbox=TRUE,pretty=TRUE)
  }
  rm(b,ef,f,features);gc()
 }
 tsv(do.call(rbind,summary),file.path(out,'CP10K_and_reference_change.tsv'));cat('PREPARE_COMPLETE\n');flush.console()
}
if(mode=='original_KO'){
 bg<-unique(scope[,c('analysis','context')]);statuses<-list()
 for(i in seq_len(nrow(bg))){a<-bg$analysis[i];c<-bg$context[i];name<-paste(a,c,sep='__');d<-file.path(out,'corrected',name)
  m<-read.delim(file.path(d,'matched_target_controls.tsv'));old<-read.delim(file.path(cfg$source_dir,'matched',name,'KO_status.tsv'))
  for(seed in cfg$seeds){graph<-NULL
   for(t in sort(unique(m$control))){hit<-old[old$seed==seed&old$target==t&old$status=='COMPLETED',]
    if(nrow(hit)==1){p<-hit$path;stopifnot(file.exists(p));s<-'COMPLETED';reuse<-TRUE}else{
     if(is.null(graph)){gp<-file.path(origin(a),c,paste0('seed_',seed),'WT_network.rds');graph<-readRDS(gp);gh<-sha(gp)}
     sd<-file.path(d,paste0('seed_',seed));dir.create(sd,showWarnings=FALSE);p<-file.path(sd,paste0(t,'_diffRegulation.tsv'))
     cat('NEW_ORIGINAL_REFERENCE_KO',name,seed,t,'\n');flush.console();s<-ko(graph$WT,t,seed,p,gh);reuse<-FALSE}
    statuses[[length(statuses)+1L]]<-data.frame(analysis=a,context=c,seed=seed,target=t,status=s,path=p,reused=reuse,network_kind='original_WT')
   };rm(graph);gc()
  }
 }
 tsv(do.call(rbind,statuses),file.path(out,'corrected_reference_KO_status.tsv'));cat('ORIGINAL_REFERENCE_KO_COMPLETE\n');flush.console()
}
if(mode=='expanded_WT'){
 seed<-as.integer(args[3]);stopifnot(seed%in%cfg$seeds);x<-readRDS(file.path(interim,'expanded_input.rds'));freeze<-fromJSON(file.path(out,'expansion_freeze.json'))
 stopifnot(freeze$input_sha256==sha(file.path(interim,'expanded_input.rds')));p<-file.path(interim,paste0('expanded_WT_seed_',seed,'.rds'))
 if(file.exists(p)){m<-fromJSON(paste0(p,'.json'));stopifnot(m$input_sha256==freeze$input_sha256,m$sha256==sha(p))}
 else{cat('BUILD_EXPANDED_WT',seed,nrow(x$counts),'\n');flush.console();wt<-build_knk_wt(x$counts,ncfg,seed);atomic_rds(wt,p)}
 wt<-readRDS(p);stopifnot(identical(rownames(wt$WT),rownames(x$counts)),all(is.finite(wt$WT)),all(diag(wt$WT)==0))
 write_json(list(status='PASS',seed=seed,sha256=sha(p),input_sha256=freeze$input_sha256,reference_sha256=freeze$reference_sha256),paste0(p,'.json'),auto_unbox=TRUE,pretty=TRUE)
 cat('EXPANDED_WT_COMPLETE',seed,'\n');flush.console()
}
if(mode=='expanded_select'){
 b<-baseline(cfg$expanded_analysis,cfg$expanded_context);ef<-full_expression_features(b);features<-list()
 for(seed in cfg$seeds){p<-file.path(interim,paste0('expanded_WT_seed_',seed,'.rds'));g<-readRDS(p);stopifnot(sha(p)==fromJSON(paste0(p,'.json'))$sha256)
  f<-ef$features[match(rownames(g$WT),ef$features$gene),];f$out_degree<-rowSums(g$WT!=0);f$out_abs_weight_sum<-rowSums(abs(g$WT))
  features[[as.character(seed)]]<-f;tsv(f,file.path(out,paste0('expanded_WT_features_seed_',seed,'.tsv')));rm(g);gc()}
 refs<-read.delim(file.path(out,'frozen_expansion_reference_genes.tsv'))$gene
 pool<-Reduce(intersect,lapply(features,function(f)f$gene[f$out_degree>0&f$gene%in%refs]))
 stopifnot(length(pool)>0);sel<-select_references(features,cfg$expanded_targets,pool,cfg$matched_controls)
 tsv(sel$selected,file.path(out,'expanded_matched_target_controls.tsv'));tsv(sel$ranking,file.path(out,'expanded_complete_reference_ranking.tsv'))
 write_json(list(status='FROZEN_BEFORE_REFERENCE_KO',both_seed_estimable_reference_pool=length(pool),selected_reference_genes=sort(unique(sel$selected$control)),
   selection_sha256=sha(file.path(out,'expanded_matched_target_controls.tsv'))),file.path(out,'expanded_selection_freeze.json'),auto_unbox=TRUE,pretty=TRUE)
 cat('EXPANDED_SELECTION_COMPLETE',length(pool),'\n');flush.console()
}
if(mode=='expanded_KO'){
 seed<-as.integer(args[3]);p<-file.path(interim,paste0('expanded_WT_seed_',seed,'.rds'));g<-readRDS(p)
 gh<-sha(p);stopifnot(gh==fromJSON(paste0(p,'.json'))$sha256);m<-read.delim(file.path(out,'expanded_matched_target_controls.tsv'))
 stopifnot(sha(file.path(out,'expanded_matched_target_controls.tsv'))==fromJSON(file.path(out,'expanded_selection_freeze.json'))$selection_sha256)
 sd<-file.path(out,'expanded',paste0('seed_',seed));dir.create(sd,recursive=TRUE,showWarnings=FALSE);rows<-list()
 for(t in sort(unique(c(cfg$expanded_targets,m$control)))){cat('EXPANDED_KO',seed,t,'\n');flush.console();path<-file.path(sd,paste0(t,'_diffRegulation.tsv'));s<-ko(g$WT,t,seed,path,gh)
  rows[[t]]<-data.frame(analysis=cfg$expanded_analysis,context=cfg$expanded_context,seed=seed,target=t,status=s,path=path,reused=FALSE,network_kind='expanded_WT')}
 tsv(do.call(rbind,rows),file.path(out,paste0('expanded_KO_status_seed_',seed,'.tsv')));cat('EXPANDED_KO_COMPLETE',seed,'\n');flush.console()
}
