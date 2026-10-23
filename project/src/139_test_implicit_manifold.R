suppressPackageStartupMessages(library(Matrix))
source('src/138_implicit_manifold.R')
set.seed(1033);wt<-matrix(rnorm(60^2,sd=.05),60);diag(wt)<-0
dimnames(wt)<-list(paste0('G',1:60),paste0('G',1:60));target<-'G3'
ko<-wt;ko[target,]<-0
p<-nrow(wt);wx<-wt+1;wy<-ko+1;xy<-.9*(sum(wx)+sum(wy))/(2*p)*diag(p)
w<--rbind(cbind(wx,xy),cbind(t(xy),wy));diag(w)<-0;diag(w)<--colSums(w)
op<-make_implicit_ko_operator(wt,target)$operator
v<-rnorm(2*p)
product_error<-max(abs(drop(w%*%v)-op(v)))
set.seed(42);a<-scTenifoldNet::manifoldAlignment(wt,ko,d=2,nCores=1)
args<-commandArgs(trailingOnly=TRUE)
opts<-if(length(args))list(ncv=as.integer(args[1]))else NULL
set.seed(42);b<-implicit_ko_alignment(wt,target,opts=opts)
distance<-function(z)sqrt(rowSums((z[seq_len(p),]-z[p+seq_len(p),])^2))
dist_error<-max(abs(distance(a)-distance(b)))
result<-data.frame(matrix_product_max_error=product_error,KO_distance_max_error=dist_error,pass=product_error<1e-10 && dist_error<1e-7)
file<-if(length(args))paste0('implicit_manifold_ncv',args[1],'_tests.tsv')else'implicit_manifold_tests.tsv'
write.table(result,file.path('results/2026-10-04-stratified-virtual-ko-033',file),sep='\t',quote=FALSE,row.names=FALSE)
print(result);stopifnot(result$pass)
