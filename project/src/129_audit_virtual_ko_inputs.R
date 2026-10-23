options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(jsonlite)})
out <- 'results/2026-10-04-stratified-virtual-ko-033'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
writeLines(capture.output(sessionInfo()), file.path(out,'audit_sessionInfo.txt'))
writeLines(capture.output({print(formals(scTenifoldKnk::scTenifoldKnk)); print(scTenifoldKnk::scTenifoldKnk)}),file.path(out,'installed_knk_function.txt'))
message('Loading existing D3 Seurat object')
x <- readRDS('data/processed/D3_GSE160936_glia_seurat.rds')
md <- x[[]]
writeLines(capture.output({print(dim(x)); print(names(md)); print(table(md$cell_type,md$disease)); print(SeuratObject::Layers(x[['RNA']]))}), file.path(out,'input_object_audit.txt'))
counts <- tryCatch(SeuratObject::GetAssayData(x, assay='RNA', layer='counts'), error=function(e) NULL)
if(is.null(counts)) { x <- SeuratObject::JoinLayers(x,assay='RNA'); counts <- SeuratObject::GetAssayData(x,assay='RNA',layer='counts') }
stopifnot(inherits(counts,'sparseMatrix'), identical(colnames(counts),rownames(md)), all(counts@x>=0),all(counts@x==round(counts@x)))
targets <- c('AKR1B1','LTA4H','ALOX15B','CA2','SNCA','PKM','GAPDH','SLC2A1','MMP9','RELA','DPP4','NLRP3','PTGS2','DYRK1A','SETD7','ALOX5','ALOX12','ALOX15','MPO','ABCC1')
rows <- list()
dir.create('data/interim/2026-10-04-stratified-virtual-ko-033',recursive=TRUE,showWarnings=FALSE)
for(ct in c('Microglia','Astrocyte')) for(dis in c('AD','Control')) {
  idx <- which(md$cell_type==ct & md$disease==dis)
  z <- counts[,idx,drop=FALSE]
  rows[[paste(ct,dis)]] <- do.call(rbind,lapply(targets,function(g) {
    present <- g %in% rownames(z)
    data.frame(cell_type=ct,disease=dis,target=g,n_cells=ncol(z),n_donors=length(unique(md$donor_id[idx])),present=present,n_expressing=if(present) sum(z[g,]>0) else 0, fraction=if(present)mean(z[g,]>0)else 0)
  }))
  saveRDS(list(counts=z,metadata=md[idx,,drop=FALSE]),file.path('data/interim/2026-10-04-stratified-virtual-ko-033',paste0(ct,'_',dis,'_raw_counts.rds')))
}
write.table(do.call(rbind,rows),file.path(out,'target_input_coverage.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write_json(list(status='INPUT_AUDIT_COMPLETE',raw_integer_counts=TRUE,source='D3_GSE160936',target_count=length(targets),packages=lapply(c('Seurat','scTenifoldKnk','scTenifoldNet'),function(p)list(package=p,version=as.character(packageVersion(p))))),file.path(out,'input_audit.json'),pretty=TRUE,auto_unbox=TRUE)
message('INPUT_AUDIT_COMPLETE')
