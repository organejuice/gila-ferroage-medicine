# New task/output adapter for audited 033 functions. No original outputs are changed.
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(Matrix);library(yaml);library(jsonlite);library(digest)})
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)>=2)
cfg<-read_yaml(args[1]);mode<-args[2];out<-cfg$output_dir;interim<-cfg$interim_dir
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(interim,recursive=TRUE,showWarnings=FALSE)
ncfg<-read_yaml(cfg$network_config);rcfg<-read_yaml(cfg$retention_config)
scope<-read.delim(cfg$scope_path,check.names=FALSE,encoding='UTF-8')
stopifnot(grepl('analysis$',names(scope)[1]));names(scope)[1]<-'analysis' # BOM under Windows C locale.
stopifnot(nrow(scope)==29,length(unique(scope$target))==15,all(c('analysis','context','target')%in%names(scope)))
source('src/130_virtual_ko_cached_backend.R')
write_tsv<-function(x,p)write.table(x,p,sep='\t',quote=FALSE,row.names=FALSE,na='')
atomic_rds<-function(x,p){saveRDS(x,paste0(p,'.partial'));stopifnot(file.rename(paste0(p,'.partial'),p))}
sha<-function(p)digest(file=p,algo='sha256')
panels<-list(Astrocyte=c('AQP4','ALDH1L1','SLC1A2','SLC1A3','SOX9','GFAP'),
 Microglia=c('P2RY12','TMEM119','CX3CR1','CSF1R','AIF1'),
 Neuron=c('SNAP25','RBFOX3','SYT1','STMN2','GRIN1'),
 Oligodendrocyte=c('MBP','PLP1','MOG','MOBP','OLIG1','OLIG2'),
 OPC=c('PDGFRA','CSPG4','VCAN'),Endothelial=c('CLDN5','KDR','EMCN','FLT1'),
 Pericyte=c('RGS5','PDGFRB','CSPG4','MCAM'))
targets<-sort(unique(scope$target))
sample_instances<-function(z,md,donors,seed){
 set.seed(seed);columns<-integer();mm<-list()
 for(i in seq_along(donors)){
  indices<-which(md$donor_id==donors[i]);ss<-split(indices,as.character(md$region[indices]))
  kept<-sort(unlist(lapply(ss,function(ix)if(length(ix)>ncfg$max_cells_per_donor_region)sample(ix,ncfg$max_cells_per_donor_region)else ix)))
  m<-md[kept,,drop=FALSE];m$source_cell_id<-rownames(m);m$source_donor_id<-m$donor_id
  m$donor_instance_id<-paste0('instance_',i);m$instance_cell_id<-paste0('instance_',i,'__',m$source_cell_id)
  columns<-c(columns,kept);mm[[i]]<-m
 }
 m<-do.call(rbind,mm);rownames(m)<-m$instance_cell_id;x<-z[,columns,drop=FALSE];colnames(x)<-rownames(m)
 list(counts=x,metadata=m)
}
filter_genes<-function(z,analysis){
 keep<-Matrix::rowMeans(z>0)>=ncfg$gene_min_fraction
 if(analysis=='target_retention')keep<-keep|(rownames(z)%in%unique(c(ncfg$targets_microglia,ncfg$targets_astrocyte)) & Matrix::rowSums(z>0)>=rcfg$target_retention_min_cells)
 # Caller supplies cell-specific whitelist below, preserving original background.
 keep
}
cell_whitelist<-function(context)if(startsWith(context,'Microglia'))ncfg$targets_microglia else ncfg$targets_astrocyte
filtered_input<-function(x,analysis,context){
 z<-x$counts;keep<-Matrix::rowMeans(z>0)>=ncfg$gene_min_fraction
 if(analysis=='target_retention')keep<-keep|(rownames(z)%in%cell_whitelist(context)&Matrix::rowSums(z>0)>=rcfg$target_retention_min_cells)
 z[keep,,drop=FALSE]
}
load_qc<-function(context){
 x<-readRDS(file.path(cfg$input_dir,paste0(context,'_raw_counts.rds')))
 z<-x$counts;md<-x$metadata
 stopifnot(inherits(z,'sparseMatrix'),identical(colnames(z),rownames(md)),!anyDuplicated(colnames(z)),all(z@x>=0),all(z@x==round(z@x)),all(c('donor_id','region','sample_id')%in%names(md)))
 lib<-Matrix::colSums(z);mt<-grepl('^MT-',rownames(z),ignore.case=TRUE)
 mf<-Matrix::colSums(z[mt,,drop=FALSE])/lib;keep<-lib>=ncfg$cell_min_counts & is.finite(mf)&mf<=ncfg$cell_max_mt_fraction
 list(counts=z[,keep,drop=FALSE],metadata=md[keep,,drop=FALSE],before=x)
}
if(mode=='prepare'){
 expression<-list();markers<-list();cells<-list();availability<-list();tasks<-list();ordinal<-0L
 capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
 for(context in unique(scope$context)){
  cat('PREPARE',context,'\n');flush.console();x<-load_qc(context);stopifnot(length(unique(x$metadata$donor_id))==6L)
  baseline<-sample_instances(x$counts,x$metadata,sort(unique(x$metadata$donor_id)),42)
  for(stage in c('before_QC','after_QC','balanced')){
   xx<-switch(stage,before_QC=x$before,after_QC=x,balanced=baseline);z<-xx$counts;md<-xx$metadata
   lib<-Matrix::colSums(z);group<-split(seq_len(ncol(z)),paste(md$donor_id,md$region,sep='__'))
   for(g in targets){
    present<-g%in%rownames(z);v<-if(present)as.numeric(z[g,])else rep(NA_real_,ncol(z));cp<-10000*v/lib
    for(ix in group){k<-sum(v[ix]>0,na.rm=TRUE);N<-length(ix)
     expression[[length(expression)+1L]]<-data.frame(stage=stage,context=context,donor_id=md$donor_id[ix[1]],region=md$region[ix[1]],target=g,present=present,N=N,k=if(present)k else NA,fraction=if(present)k/N else NA,target_UMI=if(present)sum(v[ix])else NA,mean_CP10K=if(present)mean(cp[ix])else NA,mean_positive_CP10K=if(present&&k>0)mean(cp[ix][v[ix]>0])else NA)
    }
   }
  }
  z<-x$counts;md<-x$metadata;lib<-Matrix::colSums(z)
  scores<-sapply(panels,function(p){p<-intersect(p,rownames(z));if(!length(p))rep(NA_real_,ncol(z))else Matrix::colMeans(log1p(t(t(z[p,,drop=FALSE])/lib)*10000))})
  detection<-sapply(panels,function(p){p<-intersect(p,rownames(z));if(!length(p))rep(NA_real_,ncol(z))else Matrix::colSums(z[p,,drop=FALSE]>0)})
  cells[[context]]<-data.frame(source_cell_id=colnames(z),context=context,md[,intersect(c('donor_id','region','sample_id','seurat_clusters','cell_type'),names(md)),drop=FALSE],scores,setNames(as.data.frame(detection),paste0(names(panels),'_n_detected_markers')),identity_status='DESCRIPTIVE_MULTIMARKER_REVIEW_NO_AUTOMATIC_RELABEL',check.names=FALSE)
  gg<-split(seq_len(ncol(z)),paste(md$donor_id,md$region,md$seurat_clusters,sep='__'))
  for(g in unique(c(unlist(panels),'INPP5D'))){if(!g%in%rownames(z))next;v<-as.numeric(z[g,]);for(ix in gg)markers[[length(markers)+1L]]<-data.frame(context=context,donor_id=md$donor_id[ix[1]],region=md$region[ix[1]],cluster=md$seurat_clusters[ix[1]],gene=g,N=length(ix),k=sum(v[ix]>0),fraction=mean(v[ix]>0),mean_logCP10K=mean(log1p(10000*v[ix]/lib[ix])))}
  available_tool<-requireNamespace('scDblFinder',quietly=TRUE)
  samples<-sort(unique(md$sample_id));archives<-basename(list.files('data/raw/D3_GSE160936',pattern='filtered_feature_bc_matrix',full.names=TRUE))
  availability[[context]]<-data.frame(context=context,sample_id=samples,full_cell_filtered_archive=vapply(samples,function(s)any(startsWith(archives,paste0(s,'_'))),logical(1)),unfiltered_droplets_available=FALSE,scDblFinder_installed=available_tool,doublet_status=if(available_tool)'PENDING_FULL_LIBRARY_QC'else'NOT_ASSESSABLE_IN_LOCKED_ENV_TOOL_UNAVAILABLE',ambient_status='NO_UNFILTERED_DROPLET_INPUT_IDENTIFIED')
  bg<-unique(scope[scope$context==context,c('analysis','context')])
  for(a in bg$analysis){
   # Baseline cells saved once for matching expressions; retention uses complete original whitelist.
   bpath<-file.path(interim,paste0(a,'__',context,'__baseline.rds'));if(!file.exists(bpath))atomic_rds(list(counts=filtered_input(baseline,a,context),metadata=baseline$metadata,all_gene_counts=baseline$counts),bpath)
   donors<-sort(unique(x$metadata$donor_id))
   for(stage in c('LODO','BOOTSTRAP'))for(r in seq_len(if(stage=='LODO')length(donors)else cfg$bootstrap_replicates)){
    ordinal<-ordinal+1L;sid<-paste0(a,'__',context,'__',stage,'_',sprintf('%03d',r));ss<-cfg$sampling_seed_base+ordinal
    set.seed(ss);dd<-if(stage=='LODO')donors[-r]else sample(donors,length(donors),replace=TRUE)
    tasks[[length(tasks)+1L]]<-data.frame(task_id=sid,analysis=a,context=context,stage=stage,replicate=r,sampling_seed=ss,donors_json=toJSON(dd),omitted_donor=if(stage=='LODO')donors[r]else '',targets=paste(scope$target[scope$analysis==a&scope$context==context],collapse=';'))
   }
  }
  rm(x,baseline,z,md,xx);gc()
 }
 write_tsv(do.call(rbind,expression),file.path(out,'target_expression_by_donor_region.tsv'))
 write_tsv(do.call(rbind,markers),file.path(out,'cell_identity_marker_summary.tsv'))
 write_tsv(do.call(rbind,cells),file.path(out,'cell_identity_audit.tsv'))
 write_tsv(do.call(rbind,availability),file.path(out,'library_QC_availability.tsv'))
 write_tsv(do.call(rbind,tasks),file.path(out,'network_task_registry.tsv'))
 write_json(list(status='PREPARED',network_tasks=length(tasks),bootstrap_mode='exploratory_20',matched_mode='descriptive_20',scientific_validation=FALSE),file.path(out,'preparation.json'),pretty=TRUE,auto_unbox=TRUE)
}
if(mode=='test'){
 # Test actual adapter instance semantics and full retention whitelist on deterministic fixtures.
 set.seed(10);z<-Matrix(matrix(rpois(40*650,2),40,650),sparse=TRUE);rownames(z)<-paste0('G',1:40);colnames(z)<-paste0('c',1:650)
 md<-data.frame(donor_id=rep(c('d1','d2'),each=325),region='EC',row.names=colnames(z))
 t1<-sample_instances(z,md,c('d1','d1'),77);t2<-sample_instances(z,md,c('d1','d1'),77)
 stopifnot(identical(t1,t2),!anyDuplicated(colnames(t1$counts)),length(unique(t1$metadata$source_cell_id))==325,ncol(t1$counts)==650,length(unique(t1$metadata$source_donor_id))==1)
 zz<-rbind(z,Matrix(matrix(c(rep(1,30),rep(0,620)),nrow=1),sparse=TRUE));rownames(zz)[nrow(zz)]<-'MPO'
 stopifnot(!'MPO'%in%rownames(filtered_input(list(counts=zz),'primary','Microglia_AD')),'MPO'%in%rownames(filtered_input(list(counts=zz),'target_retention','Microglia_AD')))
 fp<-file.path(interim,'adapter_JSON_fixture.tsv');write_tsv(data.frame(donors_json=toJSON(c('d1','d1'))),fp)
 stopifnot(identical(fromJSON(read.delim(fp,quote='')$donors_json),c('d1','d1')))
 # Compare implicit and dense deletion through the new interface on the same small WT.
 set.seed(9);wt<-matrix(runif(70^2),70,70);wt[wt<.9]<-0;diag(wt)<-0;dimnames(wt)<-list(paste0('g',1:70),paste0('g',1:70));saved<-wt
 tc<-ncfg;tc$implicit_manifold<-FALSE;d<-delete_knk_gene(wt,'g1',tc,42,FALSE);tc$implicit_manifold<-TRUE;i<-delete_knk_gene(wt,'g1',tc,42,FALSE)
 stopifnot(identical(wt,saved),d$status=='COMPLETED',i$status=='COMPLETED')
 q<-merge(d$diffRegulation,i$diffRegulation,by='gene',suffixes=c('_dense','_implicit'));error<-max(abs(q$distance_dense-q$distance_implicit));stopifnot(error<1e-5)
 write_json(list(status='PASS',duplicated_donor_instance_test=TRUE,deterministic_sampling=TRUE,full_retention_whitelist_test=TRUE,task_JSON_roundtrip=TRUE,WT_unchanged=TRUE,dense_implicit_max_distance_error=error,adapter_sha256=sha('src/verify_virtual_ko_035.R')),file.path(out,'adapter_deterministic_tests.json'),pretty=TRUE,auto_unbox=TRUE)
 # Actual original graph replay; mandatory before new networks are scheduled.
 graph<-readRDS('results/2026-10-04-stratified-virtual-ko-033/Microglia_AD/seed_42/WT_network.rds')
 orig<-read.delim('results/2026-10-04-stratified-virtual-ko-033/Microglia_AD/seed_42/AKR1B1_diffRegulation.tsv')
 ans<-delete_knk_gene(graph$WT,'AKR1B1',ncfg,42,FALSE);j<-merge(orig,ans$diffRegulation,by='gene',suffixes=c('_old','_new'))
 err<-max(abs(j$distance_old-j$distance_new));stopifnot(nrow(j)==nrow(orig),err<1e-5)
 write_json(list(status='PASS',original_real_WT_replay=TRUE,n_genes=nrow(j),max_distance_error=err,adapter_sha256=sha('src/verify_virtual_ko_035.R')),file.path(out,'adapter_real_replay.json'),pretty=TRUE,auto_unbox=TRUE)
}
if(mode=='network'){
 stopifnot(length(args)>=3);table<-read.delim(file.path(out,'network_task_registry.tsv'),check.names=FALSE,quote='');task<-table[table$task_id==args[3],];stopifnot(nrow(task)==1)
 td<-file.path(out,'networks',task$task_id);dir.create(td,recursive=TRUE,showWarnings=FALSE)
 x<-load_qc(task$context);xx<-sample_instances(x$counts,x$metadata,fromJSON(task$donors_json),task$sampling_seed);z<-filtered_input(xx,task$analysis,task$context)
 write_tsv(xx$metadata,file.path(td,'cell_instances.tsv'))
 sig<-digest(list(input=sha(file.path(cfg$input_dir,paste0(task$context,'_raw_counts.rds'))),cells=colnames(z),genes=rownames(z),sampling_seed=task$sampling_seed,ncfg=ncfg,rcfg=rcfg,code=sapply(c('src/verify_virtual_ko_035.R','src/130_virtual_ko_cached_backend.R','src/135_sparse_cp_equivalent.R','src/138_implicit_manifold.R','external/virtual_ko_033/pcNet_upstream.R'),sha),versions=sapply(c('scTenifoldNet','scTenifoldKnk'),function(p)as.character(packageVersion(p)))),algo='sha256')
 rm(x);gc();status<-list()
 for(seed in ncfg$seeds){
  sd<-file.path(td,paste0('seed_',seed));dir.create(sd,showWarnings=FALSE);gp<-file.path(sd,'WT_network.rds')
  if(ncol(z)<ncfg$network_n_cells){graph<-NULL}else if(file.exists(gp)){graph<-readRDS(gp);stopifnot(graph$signature==sig,graph$seed==seed)}else{
   cat('BUILD',task$task_id,seed,nrow(z),ncol(z),'\n');flush.console();graph<-build_knk_wt(z,ncfg,seed);graph$signature<-sig;graph$seed<-seed;atomic_rds(graph,gp)
  }
  for(target in strsplit(task$targets,';',fixed=TRUE)[[1]]){
   p<-file.path(sd,paste0(target,'_diffRegulation.tsv'));done<-paste0(p,'.json')
   if(file.exists(done)){a<-fromJSON(done);stopifnot(a$signature==sig,a$seed==seed);status[[length(status)+1L]]<-data.frame(task_id=task$task_id,seed=seed,target=target,status=a$status,path=if(a$status=='COMPLETED')p else '',error='');next}
   a<-if(is.null(graph))list(status='INSUFFICIENT_CELLS')else if(!target%in%rownames(graph$WT))list(status='LOW_EXPRESSION_NOT_ESTIMABLE')else tryCatch(delete_knk_gene(graph$WT,target,ncfg,seed,FALSE),error=function(e)list(status='ERROR',error=conditionMessage(e)))
   if(a$status=='ERROR')stop('KO failed: ',target,' ',a$error)
   if(a$status=='COMPLETED'){dr<-a$diffRegulation;stopifnot(!anyDuplicated(dr$gene),all(is.finite(dr$distance)),all(dr$distance>=0));write_tsv(dr,paste0(p,'.partial'));stopifnot(file.rename(paste0(p,'.partial'),p))}
   write_json(list(status=a$status,signature=sig,seed=seed,target=target),done,pretty=TRUE,auto_unbox=TRUE)
   status[[length(status)+1L]]<-data.frame(task_id=task$task_id,seed=seed,target=target,status=a$status,path=if(a$status=='COMPLETED')p else '',error='')
   cat('KO',task$task_id,seed,target,a$status,'\n');flush.console()
  }
  rm(graph);gc()
 }
 write_tsv(do.call(rbind,status),file.path(td,'KO_status.tsv'));write_json(list(status='PASS',task_id=task$task_id,signature=sig),file.path(td,'completion.json'),pretty=TRUE,auto_unbox=TRUE)
}
if(mode=='match'){
 stopifnot(length(args)>=4);analysis<-args[3];context<-args[4];tt<-scope$target[scope$analysis==analysis&scope$context==context]
 base<-readRDS(file.path(interim,paste0(analysis,'__',context,'__baseline.rds')));z<-base$counts;lib<-Matrix::colSums(z)
 expr<-data.frame(gene=rownames(z),detection_fraction=as.numeric(Matrix::rowMeans(z>0)),mean_CP10K=as.numeric(Matrix::rowMeans(t(t(z)/lib)*10000)))
 origin<-if(analysis=='primary')'results/2026-10-04-stratified-virtual-ko-033'else'results/2026-10-05-virtual-ko-target-retention-sensitivity-033'
 features<-list();directories<-character()
 for(seed in ncfg$seeds){gp<-file.path(origin,context,paste0('seed_',seed),'WT_network.rds');g<-readRDS(gp)
  stopifnot(identical(rownames(g$WT),rownames(z)))
  f<-expr;f$out_degree<-rowSums(g$WT!=0);f$out_abs_weight_sum<-rowSums(abs(g$WT));f$seed<-seed;features[[as.character(seed)]]<-f;rm(g);gc()
 }
 excluded<-unique(c(ncfg$targets_microglia,ncfg$targets_astrocyte))
 pool<-Reduce(intersect,lapply(features,function(f)f$gene[f$out_degree>0&!f$gene%in%excluded]))
 rows<-list()
 for(target in tt){
  dd<-list()
  for(seed in ncfg$seeds){f<-features[[as.character(seed)]];f<-f[f$gene%in%c(pool,target),];cols<-c('detection_fraction','mean_CP10K','out_degree','out_abs_weight_sum');x<-as.matrix(f[,cols]);pr<-apply(x,2,function(v)rank(v,ties.method='average')/(length(v)+1));tidx<-match(target,f$gene)
   dist<-sqrt(rowMeans(sweep(pr,2,pr[tidx,],'-')^2));ref<-which(f$gene%in%pool)
   q<-data.frame(target=target,control=f$gene[ref],seed=seed,distance=dist[ref],matching_limited=any(vapply(cols,function(n)f[tidx,n]<min(f[ref,n])||f[tidx,n]>max(f[ref,n]),logical(1))))
   for(n in cols){q[[paste0('target_',n)]]<-f[tidx,n];q[[paste0('control_',n)]]<-f[ref,n]};dd[[as.character(seed)]]<-q
  }
  d<-merge(dd[['42']][,c('control','distance')],dd[['43']][,c('control','distance')],by='control');d$distance<-(d$distance.x+d$distance.y)/2;chosen<-head(d$control[order(d$distance,d$control)],cfg$matched_controls)
  q<-do.call(rbind,lapply(dd,function(q)q[q$control%in%chosen,]));q$analysis<-analysis;q$context<-context;q$interpretation<-'DESCRIPTIVE_NEAREST_REFERENCE_NOT_CAUSAL_NULL';rows[[target]]<-q
 }
 md<-file.path(out,'matched',paste0(analysis,'__',context));dir.create(md,recursive=TRUE,showWarnings=FALSE)
 write_tsv(do.call(rbind,rows),file.path(md,'matched_target_controls.tsv'));controls<-sort(unique(do.call(rbind,rows)$control));statuses<-list()
 for(seed in ncfg$seeds){gp<-file.path(origin,context,paste0('seed_',seed),'WT_network.rds');graph<-readRDS(gp);sd<-file.path(md,paste0('seed_',seed));dir.create(sd,showWarnings=FALSE)
  sig<-digest(list(graph_hash=sha(gp),config_hash=sha(args[1]),adapter_hash=sha('src/verify_virtual_ko_035.R')),algo='sha256')
  for(target in controls){p<-file.path(sd,paste0(target,'_diffRegulation.tsv'));marker<-paste0(p,'.json')
   if(file.exists(marker)){a<-fromJSON(marker);stopifnot(a$signature==sig,a$status=='COMPLETED')}else{
    cat('MATCHED KO',analysis,context,seed,target,'\n');flush.console();ans<-delete_knk_gene(graph$WT,target,ncfg,seed,FALSE);stopifnot(ans$status=='COMPLETED');write_tsv(ans$diffRegulation,paste0(p,'.partial'));stopifnot(file.rename(paste0(p,'.partial'),p));write_json(list(status='COMPLETED',signature=sig),marker,pretty=TRUE,auto_unbox=TRUE)
   }
   statuses[[length(statuses)+1L]]<-data.frame(analysis=analysis,context=context,seed=seed,target=target,status='COMPLETED',path=p)
  };rm(graph);gc()
 }
 write_tsv(do.call(rbind,statuses),file.path(md,'KO_status.tsv'));write_json(list(status='PASS',n_unique_controls=length(controls)),file.path(md,'completion.json'),pretty=TRUE,auto_unbox=TRUE)
}
