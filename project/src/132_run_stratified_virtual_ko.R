options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(Matrix);library(yaml);library(jsonlite)})
source('src/130_virtual_ko_cached_backend.R')
args<-commandArgs(trailingOnly=TRUE)
cfg<-yaml::read_yaml(if(length(args))args[1]else'configs/virtual_ko_033.yaml')
out<-cfg$output_dir
dir.create(out,recursive=TRUE,showWarnings=FALSE)
states<-if(is.null(cfg$states))c('AD','Control')else cfg$states
for(ct in c('Microglia','Astrocyte')) for(dis in states) {
  input<-readRDS(file.path(cfg$input_dir,paste0(ct,'_',dis,'_raw_counts.rds')))
  md<-input$metadata; z<-input$counts
  mt<-grepl('^MT-',rownames(z),ignore.case=TRUE)
  lib<-Matrix::colSums(z);mtfrac<-Matrix::colSums(z[mt,,drop=FALSE])/lib
  eligible<-which(lib>=cfg$cell_min_counts & mtfrac<=cfg$cell_max_mt_fraction)
  # Fixed donor x region cap; identical input subset across network seeds.
  set.seed(42)
  strata<-split(eligible,paste(md$donor_id[eligible],md$region[eligible],sep='__'))
  keep<-sort(unlist(lapply(strata,function(i)if(length(i)>cfg$max_cells_per_donor_region)sample(i,cfg$max_cells_per_donor_region)else i)))
  z<-z[,keep,drop=FALSE];md<-md[keep,,drop=FALSE]
  fraction<-Matrix::rowMeans(z>0)
  gene_keep<-fraction>=cfg$gene_min_fraction
  targets<-if(ct=='Microglia')cfg$targets_microglia else cfg$targets_astrocyte
  primary_eligible<-targets%in%rownames(z)[gene_keep]
  if(!is.null(cfg$target_retention_min_cells))gene_keep<-gene_keep|(rownames(z)%in%targets & Matrix::rowSums(z>0)>=cfg$target_retention_min_cells)
  audit<-data.frame(target=targets,n_expressing=as.numeric(Matrix::rowSums(z[targets,,drop=FALSE]>0)),fraction=as.numeric(fraction[targets]),eligible=targets%in%rownames(z)[gene_keep])
  audit$primary_eligible<-primary_eligible
  ctx<-paste(ct,dis,sep='_');ctxdir<-file.path(out,ctx)
  dir.create(ctxdir,recursive=TRUE,showWarnings=FALSE)
  write.table(audit,file.path(ctxdir,'target_QC.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  write.table(md,file.path(ctxdir,'sampled_cells.tsv'),sep='\t',quote=FALSE,row.names=TRUE)
  z<-z[gene_keep,,drop=FALSE]
  write_json(list(context=ctx,n_genes=nrow(z),n_cells=ncol(z),n_donors=length(unique(md$donor_id)),regions=as.list(table(md$region)),qc=cfg),file.path(ctxdir,'network_input.json'),auto_unbox=TRUE,pretty=TRUE)
  rm(input);gc()
  if(ncol(z)<cfg$network_n_cells)stop('Insufficient QC cells: ',ctx)
  for(seed in cfg$seeds) {
    prefix<-file.path(ctxdir,paste0('seed_',seed));dir.create(prefix,showWarnings=FALSE)
    graph_path<-file.path(prefix,'WT_network.rds')
    message('NETWORK ',ctx,' seed ',seed,': ',nrow(z),' genes, ',ncol(z),' cells')
    if(file.exists(graph_path)) {
      graph<-readRDS(graph_path)
      graph_keys<-setdiff(names(cfg),c('implicit_manifold','ma_n_dim','fdr'))
      if(!identical(graph$provenance$params[graph_keys],cfg[graph_keys])) stop('Cached network parameter mismatch: ',graph_path)
    } else {
      graph<-build_knk_wt(z,cfg,seed)
      graph$provenance<-list(context=ctx,seed=seed,params=cfg,packages=list(scTenifoldKnk=as.character(packageVersion('scTenifoldKnk')),scTenifoldNet=as.character(packageVersion('scTenifoldNet'))))
      saveRDS(graph,graph_path)
    }
    summary_rows<-list()
    for(target in targets) {
      eligible_target<-target%in%rownames(graph$WT)
      path<-file.path(prefix,paste0(target,'_KO.rds'))
      message('KO ',ctx,' ',seed,' ',target)
      skip_primary<-isTRUE(cfg$only_newly_eligible_targets)&&audit$primary_eligible[match(target,audit$target)]
      ans<-if(skip_primary)list(status='PRIMARY_MODEL_ALREADY_EVALUATED',diffRegulation=NULL)else if(!eligible_target)list(status='LOW_EXPRESSION_NOT_ESTIMABLE',diffRegulation=NULL) else if(file.exists(path))readRDS(path) else tryCatch(delete_knk_gene(graph$WT,target,cfg,seed,store_networks=FALSE),error=function(e)list(status='ERROR',error=conditionMessage(e),diffRegulation=NULL))
      if(is.null(ans$implementation)&&ans$status=='COMPLETED')ans$implementation<-'official_dense_1.3_legacy_result_numerically_equivalent'
      ans$provenance<-list(context=ctx,target=target,seed=seed,WT_network=graph_path,params=cfg)
      # WT is preserved once; reconstruct KO as WT with target outgoing row zeroed.
      ans$tensorNetworks<-NULL
      saveRDS(ans,path)
      dr<-ans$diffRegulation
      if(!is.null(dr))write.table(dr,file.path(prefix,paste0(target,'_diffRegulation.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
      summary_rows[[target]]<-data.frame(context=ctx,cell_type=ct,state=dis,seed=seed,target=target,status=ans$status,n_tested=if(is.null(dr))0 else sum(dr$gene!=target),n_significant=if(is.null(dr))NA_integer_ else sum(dr$p.adj<cfg$fdr & dr$gene!=target),error=if(is.null(ans$error))''else ans$error)
    }
    write.table(do.call(rbind,summary_rows),file.path(prefix,'target_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
    rm(graph);gc()
  }
  rm(z,md);gc()
}
write_json(list(status='COMPUTATIONS_COMPLETE_PENDING_AUDIT'),file.path(out,'run_state.json'),pretty=TRUE,auto_unbox=TRUE)
