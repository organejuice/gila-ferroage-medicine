# Describe target-expressing cells using the already frozen marker panel.
# No relabeling, removal, doublet substitution, corrected counts or new network.
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(Matrix);library(yaml);library(jsonlite)})
cfg<-read_yaml('configs/virtual_ko_scientific_review_035.yaml');out<-cfg$output_dir
old<-read_yaml(cfg$source_config);scope<-read.delim(cfg$scope_path,check.names=FALSE,encoding='UTF-8');names(scope)[1]<-'analysis'
identity<-read.delim(file.path(old$output_dir,'cell_identity_audit.tsv'),check.names=FALSE)
targets<-sort(unique(scope$target));results<-list();clusters<-list()
capture.output(sessionInfo(),file=file.path(out,'identity_review_R_sessionInfo.txt'))
for(context in unique(scope$context)){
  x<-readRDS(file.path(old$input_dir,paste0(context,'_raw_counts.rds')))
  stopifnot(identical(colnames(x$counts),rownames(x$metadata)))
  id<-identity[identity$context==context,];ix<-match(id$source_cell_id,colnames(x$counts));stopifnot(!anyNA(ix))
  z<-x$counts[,ix,drop=FALSE];stopifnot(identical(colnames(z),id$source_cell_id))
  marker_names<-c('CSF1R','P2RY12','INPP5D','AQP4','ALDH1L1','SLC1A2','SLC1A3','SOX9','GFAP')
  marker_detected<-as.matrix(z[intersect(marker_names,rownames(z)),,drop=FALSE]>0)
  target_detected<-as.matrix(z[intersect(targets,rownames(z)),,drop=FALSE]>0)
  original<-strsplit(context,'_',fixed=TRUE)[[1]][1];alternative<-if(original=='Astrocyte')'Microglia'else'Astrocyte'
  for(target in targets){
    present<-target%in%rownames(target_detected);detected<-if(present)target_detected[target,] else rep(FALSE,ncol(z))
    pos<-which(detected);N<-ncol(z);k<-length(pos)
    for(level in c('context','donor_region_cluster')){
      groups<-if(level=='context')list(all=seq_len(N))else split(seq_len(N),paste(id$donor_id,id$region,id$seurat_clusters,sep='__'))
      for(g in names(groups)){
        jj<-groups[[g]];pp<-jj[detected[jj]]
        q<-data.frame(context=context,target=target,level=level,group=g,N=length(jj),k=length(pp),
          donor_id=if(level=='context')''else id$donor_id[jj[1]],region=if(level=='context')''else id$region[jj[1]],
          cluster=if(level=='context')''else as.character(id$seurat_clusters[jj[1]]),
          positive_cells_original_panel_mean=if(length(pp))mean(id[pp,original])else NA_real_,
          positive_cells_alternative_panel_mean=if(length(pp))mean(id[pp,alternative])else NA_real_,
          positive_cells_alternative_higher=if(length(pp))sum(id[pp,alternative]>id[pp,original])else NA_integer_,
          positive_cells_both_panels_two_markers=if(length(pp))sum(id$Astrocyte_n_detected_markers[pp]>=2&id$Microglia_n_detected_markers[pp]>=2)else NA_integer_,
          interpretation='DESCRIPTIVE_TARGET_POSITIVE_CELLS_NOT_IDENTITY_PROBABILITY_OR_DOUBLET_RATE')
        for(marker in marker_names)
          q[[paste0('positive_cells_detecting_',marker)]]<-if(marker%in%rownames(marker_detected)&&length(pp))sum(marker_detected[marker,pp])else NA_integer_
        results[[length(results)+1L]]<-q
      }
    }
  }
  cat('IDENTITY_TARGET_REVIEW',context,'cells',N,'\n');flush.console();rm(x,z,id);gc()
}
write.table(do.call(rbind,results),file.path(out,'target_positive_cell_identity_review.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='')
write_json(list(status='DESCRIPTIVE_REVIEW_COMPLETE',new_labels=FALSE,doublet_detection=FALSE,ambient_correction=FALSE,new_networks=FALSE),file.path(out,'identity_target_review.json'),auto_unbox=TRUE,pretty=TRUE)
