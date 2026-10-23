# GPL >=2 adaptation of scTenifoldNet 1.3 manifoldAlignment.
# Exact matrix-vector product for its dense Laplacian, including nonsymmetric
# column-degree convention, +1 connections and 0.9 cross-network coupling.
# KO is represented as one removed WT row, without allocating duplicate graphs.
make_implicit_ko_operator <- function(wt,target) {
  stopifnot(all(diag(wt)==0))
  p<-nrow(wt);i<-match(target,rownames(wt));stopifnot(!is.na(i))
  row<-wt[i,];column_sums<-colSums(wt)
  mu<-.9*(2*sum(wt)-sum(row)+2*p*p)/(2*p)
  degree_x<-column_sums+p+mu
  degree_y<-column_sums-row+p+mu
  op<-function(v,args=NULL) {
    x<-v[seq_len(p)];y<-v[p+seq_len(p)]
    wx<-drop(wt%*%x);wy<-drop(wt%*%y);wy[i]<-0
    c(degree_x*x-wx-sum(x)-mu*y,degree_y*y-wy-sum(y)-mu*x)
  }
  list(operator=op,n=2L*p,mu=mu)
}
implicit_ko_alignment <- function(wt,target,d=2,ncores=1,opts=NULL) {
  gc()
  op<-make_implicit_ko_operator(wt,target)
  RhpcBLASctl::omp_set_num_threads(ncores);RhpcBLASctl::blas_set_num_threads(ncores)
  if(is.null(opts))opts<-list()
  eig<-suppressWarnings(RSpectra::eigs(op$operator,k=2*d,n=op$n,which='SR',opts=opts))
  if(eig$nconv<2*d)stop('Manifold eigensolver did not converge')
  values<-Re(eig$values);vectors<-Re(eig$vectors);ix<-order(values)
  vectors<-vectors[,ix,drop=FALSE];values<-values[ix]
  vectors<-vectors[,values>1e-8,drop=FALSE]
  if(ncol(vectors)<d)stop('Insufficient positive manifold eigenvalues')
  ma<-vectors[,seq_len(d),drop=FALSE]
  dimnames(ma)<-list(c(paste0('X_',rownames(wt)),paste0('Y_',rownames(wt))),paste0('NLMA ',seq_len(d)))
  attr(ma,'eigenvalues')<-values;attr(ma,'nconv')<-eig$nconv
  attr(ma,'solver_opts')<-opts;attr(ma,'solver_iterations')<-eig$niter
  ma
}
