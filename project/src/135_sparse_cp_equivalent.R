# GPL >=2 adaptation of scTenifoldNet 1.3 cpDecomposition/tensorDecomposition.
# Exact four-mode CP-ALS updates, initial RNG, L1 normalization and stop rule.
# Computes tensor contractions slice-wise; avoids dense tensors/unfoldings.
tensor_decomposition_sparse <- function(nets,K,maxIter,maxError,nDecimal) {
  p<-nrow(nets[[1]]);n<-length(nets);modes<-c(p,p,1L,n)
  set.seed(1)
  U<-lapply(modes,function(m)matrix(rnorm(m*K),nrow=m,ncol=K))
  xnorm2<-sum(vapply(nets,function(z)sum(z@x^2),numeric(1)))
  prev<-NA_real_;resids<-numeric();converged<-FALSE
  for(iter in seq_len(maxIter-1L)) {
    for(m in 1:4) {
      V<-Reduce('*',lapply(U[-m],crossprod))
      contraction<-matrix(0,modes[m],K)
      if(m==1)for(i in seq_len(n))contraction<-contraction+sweep(as.matrix(nets[[i]]%*%U[[2]]),2,U[[3]][1,]*U[[4]][i,],'*')
      if(m==2)for(i in seq_len(n))contraction<-contraction+sweep(as.matrix(Matrix::t(nets[[i]])%*%U[[1]]),2,U[[3]][1,]*U[[4]][i,],'*')
      if(m==3)for(i in seq_len(n))contraction[1,]<-contraction[1,]+colSums(U[[1]]*as.matrix(nets[[i]]%*%U[[2]]))*U[[4]][i,]
      if(m==4)for(i in seq_len(n))contraction[i,]<-colSums(U[[1]]*as.matrix(nets[[i]]%*%U[[2]]))*U[[3]][1,]
      tmp<-contraction%*%solve(V)
      lambdas<-colSums(abs(tmp)) # original norm(as.matrix(vec)) uses norm type O
      U[[m]]<-sweep(tmp,2,lambdas,'/')
    }
    gram<-Reduce('*',lapply(U,crossprod))
    estnorm2<-sum(outer(lambdas,lambdas)*gram)
    inner<-0
    for(i in seq_len(n))inner<-inner+sum(colSums(U[[1]]*as.matrix(nets[[i]]%*%U[[2]]))*lambdas*U[[3]][1,]*U[[4]][i,])
    resid<-sqrt(max(0,xnorm2+estnorm2-2*inner));resids<-c(resids,resid)
    if(iter>1 && abs(resid-prev)/sqrt(xnorm2)<maxError){converged<-TRUE;break}
    prev<-resid
    if(iter%%10==0)message('Sparse CP iteration ',iter,' residual ',resid)
  }
  avg<-sweep(U[[1]],2,lambdas*U[[3]][1,]*colMeans(U[[4]]),'*')%*%t(U[[2]])
  avg<-round(avg/max(abs(avg)),nDecimal)
  dimnames(avg)<-dimnames(nets[[1]])
  list(X=as(avg,'dgCMatrix'),converged=converged,iterations=length(resids),all_residuals=resids,norm_explained=100*(1-tail(resids,1)/sqrt(xnorm2)))
}
