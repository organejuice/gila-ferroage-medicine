# Adapted from installed scTenifoldKnk 1.0.3 (GPL >=2), cailab-tamu/scTenifoldKnk.
# Shares the identical WT graph across single-gene deletions to avoid rebuilding it.
# Explicitly fixes the installed wrapper's unassigned qc_minCells filter externally.
build_knk_wt <- function(counts,cfg,seed) {
  set.seed(seed)
  make_nets <- scTenifoldNet::makeNetworks
  if(isTRUE(cfg$accelerated_pcnet)) {
    pc_env<-new.env(parent=globalenv())
    sys.source('external/virtual_ko_033/pcNet_upstream.R',envir=pc_env)
    net_env<-new.env(parent=environment(make_nets));net_env$pcNet<-pc_env$pcNet
    environment(make_nets)<-net_env
  }
  nets <- make_nets(X=counts,q=cfg$network_q,nNet=cfg$network_n_net,
    nCells=cfg$network_n_cells,scaleScores=TRUE,symmetric=FALSE,nComp=cfg$network_n_comp,nCores=cfg$ncores)
  if(isTRUE(cfg$sparse_cp)) {
    source('src/135_sparse_cp_equivalent.R',local=TRUE)
    td<-tensor_decomposition_sparse(nets,cfg$td_K,cfg$td_max_iter,cfg$td_max_error,cfg$td_n_decimal)
  } else td <- scTenifoldNet::tensorDecomposition(xList=nets,K=cfg$td_K,
      maxError=cfg$td_max_error,maxIter=cfg$td_max_iter,nDecimal=cfg$td_n_decimal)
  wt <- as.matrix(getFromNamespace('strictDirection','scTenifoldKnk')(td$X,lambda=0))
  diag(wt) <- 0
  list(WT=t(wt),decomposition=td[setdiff(names(td),'X')])
}
delete_knk_gene <- function(wt,target,cfg,seed,store_networks=TRUE) {
  gc()
  stopifnot(target %in% rownames(wt))
  if(!any(wt[target,]!=0)) return(list(status='NO_OUTGOING_EDGES',diffRegulation=NULL))
  set.seed(seed)
  if(isTRUE(cfg$implicit_manifold)) {
    source('src/138_implicit_manifold.R',local=TRUE)
    ma<-implicit_ko_alignment(wt,target,d=cfg$ma_n_dim,ncores=cfg$ncores)
  } else {
    ko <- wt; ko[target,] <- 0
    ma <- scTenifoldNet::manifoldAlignment(wt,ko,d=cfg$ma_n_dim,nCores=cfg$ncores)
  }
  networks<-NULL
  if(store_networks){ko<-wt;ko[target,]<-0;networks<-list(WT=Matrix::Matrix(wt,sparse=TRUE),KO=Matrix::Matrix(ko,sparse=TRUE))}
  list(status='COMPLETED',tensorNetworks=networks,implementation=if(isTRUE(cfg$implicit_manifold))'implicit_matrix_operator_equivalent_to_1.3'else'official_dense_1.3',
    manifoldAlignment=ma,diffRegulation=getFromNamespace('dRegulation','scTenifoldKnk')(ma,target))
}
