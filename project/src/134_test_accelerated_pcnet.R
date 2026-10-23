suppressPackageStartupMessages({library(Matrix);library(jsonlite)})
fast<-new.env(parent=globalenv());sys.source('external/virtual_ko_033/pcNet_upstream.R',envir=fast)
set.seed(1033)
toy<-matrix(rpois(100*200,5),100,dimnames=list(paste0('G',1:100),NULL))
real<-readRDS('data/interim/2026-10-04-stratified-virtual-ko-033/Microglia_AD_raw_counts.rds')$counts
keep<-which(Matrix::rowMeans(real>0)>=.05)
real<-as.matrix(real[sample(keep,300),sample(ncol(real),200)])
tests<-list(toy=toy,real_subset=real);rows<-list()
for(name in names(tests)) {
  z<-tests[[name]]
  a<-as.matrix(scTenifoldNet::pcNet(z,nComp=3,q=.9,nCores=1,verbose=FALSE))
  b<-as.matrix(fast$pcNet(z,nComp=3,q=.9,nCores=1,verbose=FALSE))
  max_error<-max(abs(a-b));relative_error<-sqrt(sum((a-b)^2)/sum(a^2))
  rows[[name]]<-data.frame(test=name,max_abs_error=max_error,relative_Frobenius_error=relative_error,edge_disagreement=sum((a!=0)!=(b!=0)),pass=max_error<1e-5 && relative_error<1e-5)
}
result<-do.call(rbind,rows)
write.table(result,'results/2026-10-04-stratified-virtual-ko-033/accelerated_pcnet_tests.tsv',sep='\t',quote=FALSE,row.names=FALSE)
print(result);stopifnot(all(result$pass))
