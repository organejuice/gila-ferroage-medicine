suppressPackageStartupMessages(library(Matrix))
source('src/138_implicit_manifold.R')
args<-commandArgs(trailingOnly=TRUE)
dir<-if(length(args)>=2)args[2]else'results/2026-10-04-stratified-virtual-ko-033/Microglia_AD/seed_42'
target<-if(length(args)>=3)args[3]else'AKR1B1'
seed<-if(length(args)>=4)as.integer(args[4])else 42L
wt<-readRDS(file.path(dir,'WT_network.rds'))$WT
ref<-readRDS(file.path(dir,paste0(target,'_KO.rds')))$diffRegulation
stopifnot(!is.null(ref),target%in%rownames(wt))
opts<-if(length(args))list(ncv=as.integer(args[1]))else NULL
elapsed<-system.time({set.seed(seed);ma<-implicit_ko_alignment(wt,target,opts=opts)})['elapsed']
dr<-getFromNamespace('dRegulation','scTenifoldKnk')(ma,target)
x<-merge(ref,dr,by='gene',suffixes=c('_dense','_implicit'))
relative_error<-sqrt(sum((x$distance_dense-x$distance_implicit)^2)/sum(x$distance_dense^2))
fdr_disagreement<-sum((x$p.adj_dense<.05)!=(x$p.adj_implicit<.05))
result<-data.frame(target=target,seed=seed,context_directory=dir,n_genes=nrow(x),elapsed_seconds=unname(elapsed),solver_iterations=attr(ma,'solver_iterations'),relative_distance_error=relative_error,FDR_disagreement=fdr_disagreement,pass=relative_error<1e-4 && fdr_disagreement==0)
file<-if(length(args))paste0('real_implicit_alignment_ncv',args[1],'_tests.tsv')else'real_implicit_alignment_tests.tsv'
output_dir<-if(length(args)>=2)dir else'results/2026-10-04-stratified-virtual-ko-033'
if(length(args)>=2)file<-paste0(target,'_',file)
write.table(result,file.path(output_dir,file),sep='\t',quote=FALSE,row.names=FALSE)
print(result);stopifnot(result$pass)
