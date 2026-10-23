suppressPackageStartupMessages(library(Matrix))
source('src/135_sparse_cp_equivalent.R')
set.seed(1033)
nets<-lapply(1:3,function(i){z<-matrix(rnorm(60^2),60);z[abs(z)<1]<-0;diag(z)<-0;dimnames(z)<-list(paste0('G',1:60),paste0('G',1:60));as(z,'dgCMatrix')})
tests<-list()
for(it in c(10L,100L)) {
  a<-scTenifoldNet::tensorDecomposition(nets,K=3,maxIter=it,maxError=1e-5,nDecimal=3)$X
  b<-tensor_decomposition_sparse(nets,3,it,1e-5,3)
  error<-max(abs(as.matrix(a)-as.matrix(b$X)))
  tests[[as.character(it)]]<-data.frame(max_iter=it,max_abs_error=error,pass=error<=.001,converged=b$converged,iterations=b$iterations)
}
result<-do.call(rbind,tests)
write.table(result,'results/2026-10-04-stratified-virtual-ko-033/sparse_cp_tests.tsv',sep='\t',quote=FALSE,row.names=FALSE)
print(result);stopifnot(all(result$pass))
