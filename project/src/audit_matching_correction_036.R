# Post-computation checks independent of the execution adapter.
suppressPackageStartupMessages({library(Matrix);library(yaml);library(jsonlite);library(digest)})
cfg<-read_yaml('configs/virtual_ko_matching_correction_036.yaml');out<-cfg$output_dir
sha<-function(p){con<-file(p,'rb');on.exit(close(con));digest(readBin(con,'raw',n=file.info(p)$size),algo='sha256',serialize=FALSE)}
input<-readRDS(file.path(cfg$interim_dir,'expanded_input.rds'));freeze<-fromJSON(file.path(out,'expansion_freeze.json'))
stopifnot(sha(file.path(cfg$interim_dir,'expanded_input.rds'))==freeze$input_sha256,
 all(input$original_genes%in%rownames(input$counts)),all(input$reference_genes%in%rownames(input$counts)),all(cfg$expanded_targets%in%rownames(input$counts)))
rows<-list()
for(seed in cfg$seeds){p<-file.path(cfg$interim_dir,paste0('expanded_WT_seed_',seed,'.rds'));m<-fromJSON(paste0(p,'.json'));w<-readRDS(p)
 stopifnot(m$seed==seed,m$input_sha256==freeze$input_sha256,m$reference_sha256==freeze$reference_sha256,sha(p)==m$sha256,
 identical(rownames(w$WT),rownames(input$counts)),identical(rownames(w$WT),colnames(w$WT)),all(is.finite(w$WT)),all(diag(w$WT)==0))
 rows[[as.character(seed)]]<-data.frame(seed=seed,genes=nrow(w$WT),cells=ncol(input$counts),CP_converged=w$decomposition$converged,
 CP_iterations=w$decomposition$iterations,CP_final_residual=tail(w$decomposition$all_residuals,1),norm_explained=w$decomposition$norm_explained,sha256=m$sha256)
 rm(w);gc()
}
write.table(do.call(rbind,rows),file.path(out,'expanded_network_audit.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('PASS: network hashes, both seeds, identical cells/genes, original genes and whitelist preserved, finite weights and zero diagonals\n')
