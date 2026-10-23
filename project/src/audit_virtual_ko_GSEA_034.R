# Replay enrichment scores/leading edges from raw distances for every job/library/mode.
suppressPackageStartupMessages({library(yaml);library(fgsea);library(jsonlite)})
cfg<-read_yaml('configs/virtual_ko_pathways_034.yaml');out<-cfg$output_dir
jobs<-read.delim(file.path(out,'job_registry.tsv'));libs<-list()
for(lib in cfg$libraries){x<-strsplit(readLines(file.path(cfg$library_dir,paste0(lib,'.gmt')),warn=FALSE),'\t');libs[[lib]]<-setNames(lapply(x,function(s)unique(s[-c(1,2)])),sapply(x,`[`,1))}
checks<-list()
for(i in seq_len(nrow(jobs))){
 r<-jobs[i,];dr<-read.delim(r$path);dr<-dr[dr$gene!=r$target,];dr<-dr[order(-dr$distance,dr$gene),];st<-setNames(dr$distance,dr$gene)
 tab<-read.delim(gzfile(file.path(out,'jobs',r$job_id,'GSEA_all.tsv.gz')),check.names=FALSE)
 for(lib in cfg$libraries)for(mode in names(cfg$gsea$modes)){
  x<-tab[tab$library==lib&tab$mode==mode,];stopifnot(nrow(x)>0)
  # Deterministically inspect strongest, weakest, middle and extreme-ES terms.
  candidates<-unique(c(order(x$pval)[1],order(-x$pval)[1],ceiling(nrow(x)/2),which.max(x$ES),which.min(x$ES)))
  for(j in candidates){
   row<-x[j,];selected<-sort(which(names(st)%in%libs[[lib]][[row$pathway]]))
   actual<-calcGseaStat(st,selected,gseaParam=cfg$gsea$modes[[mode]],scoreType='pos',returnLeadingEdge=TRUE)
   expected<-if(is.na(row$leadingEdge)||row$leadingEdge=='')character()else strsplit(row$leadingEdge,';',fixed=TRUE)[[1]]
   err<-abs(actual$res-row$ES);stopifnot(err<1e-10,setequal(names(st)[actual$leadingEdge],expected))
   checks[[length(checks)+1]]<-data.frame(job_id=r$job_id,library=lib,mode=mode,term=row$pathway,ES_error=err,leading_edge_exact=TRUE)
  }
 }
}
z<-do.call(rbind,checks);write.table(z,file.path(out,'GSEA_score_replay_audit.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write_json(list(status='PASS',jobs=nrow(jobs),library_mode_groups=nrow(jobs)*length(cfg$libraries)*length(cfg$gsea$modes),selected_terms_replayed=nrow(z),max_ES_error=max(z$ES_error),leading_edge_exact=all(z$leading_edge_exact),scope='Selected representative score/leading-edge replay for every job/library/mode; full scope/ORA replay separately audited'),file.path(out,'GSEA_score_replay_audit.json'),pretty=TRUE,auto_unbox=TRUE)
cat('PASS GSEA score replay',nrow(z),'terms in',nrow(jobs),'jobs; maxESerror',max(z$ES_error),'\n')
