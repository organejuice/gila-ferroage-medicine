full_expression_features <- function(base) {
 z<-base$counts;all<-base$all_gene_counts
 stopifnot(inherits(z,'sparseMatrix'),inherits(all,'sparseMatrix'),identical(colnames(z),colnames(all)),
   !anyDuplicated(colnames(z)),all(rownames(z)%in%rownames(all)))
 lib<-Matrix::colSums(all);filtered<-Matrix::colSums(z)
 stopifnot(all(is.finite(lib)),all(lib>0),all(lib>=filtered),
   identical(as.numeric(z@x),as.numeric(all[rownames(z),,drop=FALSE]@x)))
 f<-data.frame(gene=rownames(all),detected_cells=as.numeric(Matrix::rowSums(all>0)),
   detection_fraction=as.numeric(Matrix::rowMeans(all>0)),
   mean_CP10K=as.numeric(all %*% (10000/lib))/ncol(all))
 f$old_mean_CP10K<-as.numeric(all %*% (10000/filtered))/ncol(all)
 list(features=f,libraries=data.frame(cell_id=colnames(all),full_library_UMI=lib,
   filtered_library_UMI=filtered,retained_UMI_fraction=filtered/lib))
}
select_references <- function(features,targets,pool,M=20) {
 cols<-c('detection_fraction','mean_CP10K','out_degree','out_abs_weight_sum');rows<-list();orders<-list()
 for(target in targets){dd<-list()
  for(seed in names(features)){
   f<-features[[seed]];f<-f[f$gene%in%c(pool,target),]
   stopifnot(sum(f$gene==target)==1,nrow(f)>1)
   pr<-apply(as.matrix(f[,cols]),2,function(v)rank(v,ties.method='average')/(length(v)+1))
   ix<-match(target,f$gene);dist<-sqrt(rowMeans(sweep(pr,2,pr[ix,],'-')^2));ref<-which(f$gene%in%pool)
   q<-data.frame(target=target,control=f$gene[ref],seed=as.integer(seed),distance=dist[ref])
   for(n in cols){q[[paste0('target_',n)]]<-f[ix,n];q[[paste0('control_',n)]]<-f[ref,n]}
   q$matching_limited<-any(vapply(cols,function(n)f[ix,n]<min(f[ref,n])||f[ix,n]>max(f[ref,n]),logical(1)))
   dd[[seed]]<-q
  }
  d<-merge(dd[['42']][,c('control','distance')],dd[['43']][,c('control','distance')],by='control')
  d$mean_distance<-(d$distance.x+d$distance.y)/2;d<-d[order(d$mean_distance,d$control),]
  chosen<-head(d$control,M);d$target<-target;d$selected<-d$control%in%chosen;orders[[target]]<-d
  q<-do.call(rbind,lapply(dd,function(q)q[q$control%in%chosen,]));q$insufficient_controls<-length(chosen)<M
  rows[[target]]<-q
 }
 list(selected=do.call(rbind,rows),ranking=do.call(rbind,orders))
}
