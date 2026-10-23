# Adapted from audited 034 enrichment; dynamic config/output/registry and stable per-job RNG.
# Adapt the audited stage-A hypergeometric/fgsea implementation to unsigned KO distances.
options(stringsAsFactors=FALSE,width=180)
suppressPackageStartupMessages({library(yaml);library(fgsea);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)>=1);config_path<-args[1];cfg<-read_yaml(config_path);out<-cfg$output_dir
jobs<-read.delim(file.path(out,'job_registry.tsv'),check.names=FALSE)
worker<-if(length(args)>1)as.integer(args[2])else 1L;workers<-if(length(args)>2)as.integer(args[3])else 1L
stopifnot(worker>=1L,worker<=workers)
signature<-digest::digest(list(config=digest::digest(file=config_path,algo='sha256'),code=digest::digest(file='src/enrich_virtual_ko_verification_035.R',algo='sha256'),libraries=sapply(cfg$libraries,function(n)digest::digest(file=file.path(cfg$library_dir,paste0(n,'.gmt')),algo='sha256'))),algo='sha256')
read_gmt<-function(p){x<-strsplit(readLines(p,warn=FALSE),'\t');list(sets=setNames(lapply(x,function(r)sort(unique(r[-c(1,2)]))),sapply(x,`[`,1)),names=setNames(sapply(x,`[`,2),sapply(x,`[`,1)))}
libs<-setNames(lapply(cfg$libraries,function(n)read_gmt(file.path(cfg$library_dir,paste0(n,'.gmt')))),cfg$libraries)
write_gz<-function(x,p){con<-gzfile(p,'wt');on.exit(close(con));write.table(x,con,sep='\t',row.names=FALSE,quote=FALSE,na='')}
annotate<-function(x,r){for(n in c('job_id','analysis','context','cell_type','state','target','seed'))x[[n]]<-r[[n]];x}
if(worker==1L){
 capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
 write.table(data.frame(package=c('R','fgsea','yaml','jsonlite'),version=c(as.character(getRversion()),sapply(c('fgsea','yaml','jsonlite'),function(p)as.character(packageVersion(p))))),file.path(out,'environment_versions.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
}
dir.create(file.path(out,'jobs'),showWarnings=FALSE)
if(nrow(jobs)>0)for(i in seq.int(worker,nrow(jobs),by=workers)){
 r<-jobs[i,];d<-file.path(out,'jobs',r$job_id);dir.create(d,recursive=TRUE,showWarnings=FALSE)
 marker<-file.path(d,'completion.json')
 if(file.exists(marker)){a<-fromJSON(marker);stopifnot(a$source_sha256==r$sha256,a$status=='PASS',a$signature==signature);cat('CACHED',r$job_id,'\n');next}
 dr<-read.delim(r$path,check.names=FALSE);dr<-dr[dr$gene!=r$target,];stopifnot(!anyDuplicated(dr$gene),all(is.finite(dr$distance)),all(dr$distance>=0),all(is.finite(dr$p.adj)))
 dr<-dr[order(-dr$distance,dr$gene),];rank<-setNames(dr$distance,dr$gene);stopifnot(max(rank)>0)
 sig<-dr$gene[dr$p.adj<cfg$fdr]
 ora<-list();gsea<-list();qc<-list();mapping<-list()
 for(lib in names(libs)){
  full<-libs[[lib]]$sets;desc<-libs[[lib]]$names
  annotated<-unique(unlist(full,use.names=FALSE));universe<-intersect(dr$gene,annotated);query<-intersect(sig,universe)
  bgsets<-lapply(full,intersect,y=universe);sizes<-lengths(bgsets);psets<-bgsets[sizes>=cfg$ora$min_size&sizes<=cfg$ora$max_size]
  N<-length(universe);Q<-length(query);K<-lengths(psets);hits<-lapply(psets,intersect,y=query);k<-lengths(hits)
  p<-if(Q>0)phyper(k-1,K,N-K,Q,lower.tail=FALSE)else rep(1,length(K))
  z<-data.frame(term=names(psets),description=unname(desc[names(psets)]),library=lib,background_size=N,query_size=Q,term_size=K,overlap=k,gene_ratio=if(Q>0)k/Q else 0,background_ratio=K/N,fold_enrichment=if(Q>0)(k/Q)/(K/N)else NA_real_,P=p,FDR_library=p.adjust(p,'BH'),genes=sapply(hits,paste,collapse=';'),row.names=NULL)
  ora[[lib]]<-annotate(z,r)
  mapping[[lib]]<-annotate(data.frame(gene=dr$gene,library=lib,annotated=dr$gene%in%annotated,significant=dr$gene%in%sig),r)
  for(mode in names(cfg$gsea$modes)){
   set.seed(r$enrichment_seed)
   warnings<-character()
   z<-withCallingHandlers(as.data.frame(fgseaSimple(pathways=full,stats=rank,minSize=cfg$gsea$min_size,maxSize=cfg$gsea$max_size,nperm=cfg$gsea$nperm,nproc=cfg$gsea$workers,scoreType=cfg$gsea$scoreType,gseaParam=cfg$gsea$modes[[mode]])),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
   z$leadingEdge<-sapply(z$leadingEdge,paste,collapse=';');z$library<-lib;z$mode<-mode;z$description<-unname(desc[z$pathway]);z$FDR_library<-p.adjust(z$pval,'BH');z<-annotate(z,r)
   gsea[[paste(lib,mode)]]<-z
   qc[[paste(lib,mode)]]<-annotate(data.frame(library=lib,mode=mode,tested_genes=nrow(dr),annotated_genes=length(universe),significant_total=length(sig),annotated_significant=Q,ora_tests=length(psets),gsea_tests=nrow(z),missing_P=sum(is.na(z$pval)),tied_rank_values=sum(duplicated(rank)),zero_distances=sum(rank==0),p_at_resolution=sum(z$nMoreExtreme==0,na.rm=TRUE),warnings=paste(unique(warnings),collapse=' | ')),r)
   cat('GSEA',i,'/',nrow(jobs),r$job_id,lib,mode,nrow(z),'\n');flush.console()
  }
 }
 o<-do.call(rbind,ora);o$FDR_job_all_libraries<-p.adjust(o$P,'BH')
 g<-do.call(rbind,gsea);g$FDR_job_all_libraries<-ave(g$pval,g$mode,FUN=function(x)p.adjust(x,'BH'))
 write_gz(o,file.path(d,'ORA_all.tsv.gz'));write_gz(g,file.path(d,'GSEA_all.tsv.gz'));write_gz(do.call(rbind,mapping),file.path(d,'gene_mapping.tsv.gz'))
 write.table(do.call(rbind,qc),file.path(d,'QC.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='')
 write_json(list(status='PASS',source_sha256=r$sha256,signature=signature,ora_tests=nrow(o),gsea_tests=nrow(g),job_id=r$job_id),marker,auto_unbox=TRUE,pretty=TRUE)
}
write_json(list(status='WORKER_COMPLETE_PENDING_AUDIT',worker=worker,workers=workers),file.path(out,paste0('worker_',worker,'_state.json')),auto_unbox=TRUE,pretty=TRUE)
