options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(jsonlite);library(yaml)})
args<-commandArgs(trailingOnly=TRUE)
cfg<-yaml::read_yaml(if(length(args))args[1]else'configs/virtual_ko_033.yaml');out<-cfg$output_dir
checks<-list();jobs<-list();issues<-character()
provenance<-jsonlite::read_json(file.path(out,'provenance.json'),simplifyVector=TRUE)
inputs<-provenance$files[grepl('^data/',provenance$files$path),,drop=FALSE]
for(i in seq_len(nrow(inputs))) {
  p<-inputs$path[i]
  if(!file.exists(p)||digest::digest(file=p,algo='sha256')!=inputs$sha256[i])issues<-c(issues,paste('Input identity mismatch:',p))
}
states<-if(is.null(cfg$states))c('AD','Control')else cfg$states
for(ct in c('Microglia','Astrocyte'))for(state in states)for(seed in cfg$seeds) {
  ctx<-paste(ct,state,sep='_');dir<-file.path(out,ctx,paste0('seed_',seed))
  if(!file.exists(file.path(dir,'WT_network.rds'))){issues<-c(issues,paste('Missing graph:',ctx,seed));next}
  graph<-readRDS(file.path(dir,'WT_network.rds'));wt<-graph$WT
  stopifnot(nrow(wt)==ncol(wt),identical(rownames(wt),colnames(wt)),all(diag(wt)==0),all(is.finite(wt)))
  if(!isTRUE(graph$decomposition$converged))issues<-c(issues,paste('CP failed convergence:',ctx,seed))
  graph_keys<-setdiff(names(cfg),c('implicit_manifold','ma_n_dim','fdr'))
  stopifnot(identical(graph$provenance$params[graph_keys],cfg[graph_keys]),graph$provenance$context==ctx,graph$provenance$seed==seed)
  qc<-read.delim(file.path(out,ctx,'target_QC.tsv'))
  summary<-read.delim(file.path(dir,'target_summary.tsv'),check.names=FALSE)
  stopifnot(!anyDuplicated(summary$target))
  targets<-if(ct=='Microglia')cfg$targets_microglia else cfg$targets_astrocyte
  stopifnot(setequal(summary$target,targets),all(summary$context==ctx),all(summary$seed==seed))
  for(target in targets) {
    p<-file.path(dir,paste0(target,'_KO.rds'))
    if(!file.exists(p)){issues<-c(issues,paste('Missing KO:',ctx,seed,target));next}
    result<-readRDS(p);dr<-result$diffRegulation
    stopifnot(result$provenance$context==ctx,result$provenance$seed==seed,result$provenance$target==target,
      identical(result$provenance$params,cfg),qc$eligible[match(target,qc$target)]==(target%in%rownames(wt)))
    if(result$status=='COMPLETED') {
      stopifnot(target%in%rownames(wt),any(wt[target,]!=0),!is.null(dr),setequal(dr$gene,rownames(wt)),!anyDuplicated(dr$gene),
        all(is.finite(dr$distance)),all(is.finite(dr$p.value)),all(is.finite(dr$p.adj)),all(dr$distance>=0),all(dr$p.adj>=0&dr$p.adj<=1))
      ma<-result$manifoldAlignment
      stopifnot(nrow(ma)==2*nrow(wt),ncol(ma)==cfg$ma_n_dim,all(is.finite(ma)))
      replay<-getFromNamespace('dRegulation','scTenifoldKnk')(ma,target)
      ix<-match(dr$gene,replay$gene)
      stopifnot(max(abs(dr$distance-replay$distance[ix]))<1e-12,max(abs(dr$p.adj-replay$p.adj[ix]))<1e-8)
      sig<-sum(dr$p.adj<cfg$fdr&dr$gene!=target)
    } else {
      sig<-NA_integer_
      if(result$status=='LOW_EXPRESSION_NOT_ESTIMABLE')stopifnot(!target%in%rownames(wt),!qc$eligible[match(target,qc$target)])
      else if(result$status=='PRIMARY_MODEL_ALREADY_EVALUATED')stopifnot(isTRUE(cfg$only_newly_eligible_targets),qc$primary_eligible[match(target,qc$target)])
      else if(result$status=='NO_OUTGOING_EDGES')stopifnot(target%in%rownames(wt),!any(wt[target,]!=0))
      else issues<-c(issues,paste('KO error:',ctx,seed,target,result$error))
    }
    record<-summary[summary$target==target,,drop=FALSE]
    stopifnot(nrow(record)==1,record$status==result$status,
      if(is.na(sig))is.na(record$n_significant)else record$n_significant==sig,
      record$n_tested==if(is.null(dr))0 else sum(dr$gene!=target))
    jobs[[paste(ctx,seed,target)]]<-data.frame(context=ctx,seed=seed,target=target,status=result$status,n_significant=sig)
  }
  rm(wt,graph);gc()
}
expected<-length(states)*length(cfg$seeds)*(length(cfg$targets_microglia)+length(cfg$targets_astrocyte))
if(length(jobs)!=expected)issues<-c(issues,paste('Job coverage:',length(jobs),'of',expected))
write.table(do.call(rbind,jobs),file.path(out,'artifact_audit_jobs.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
audit<-list(status=if(!length(issues))'PASS'else'FAIL',expected_jobs=expected,audited_jobs=length(jobs),input_hashes_checked=nrow(inputs),issues=issues,
  validation=c('raw input identity','all target/state/seed coverage','WT dimensions and finite values','CP convergence','eligibility and outgoing-edge gates','KO artifact readability','saved alignment and DR replay','self-target excluded from significant counts'))
jsonlite::write_json(audit,file.path(out,'artifact_audit.json'),pretty=TRUE,auto_unbox=TRUE)
print(audit);if(length(issues))stop('Artifact audit failed')
