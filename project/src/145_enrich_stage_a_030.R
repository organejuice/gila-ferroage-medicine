options(stringsAsFactors=FALSE, width=220)
suppressPackageStartupMessages({library(jsonlite);library(AnnotationDbi);library(org.Hs.eg.db);library(GO.db);library(fgsea);library(Matrix)})
cfg<-fromJSON('configs/stage_a_030.json'); run<-cfg$run_id
int<-file.path('data/interim',run); raw<-file.path('data/raw',run); out<-file.path('results',run)
read<-function(n) read.delim(file.path(int,n),check.names=FALSE)
write<-function(x,n) {p<-file.path(out,n);if(grepl('gz$',n)){con<-gzfile(p,'wt');on.exit(close(con));write.table(x,con,sep='\t',row.names=FALSE,quote=FALSE,na='')}else write.table(x,p,sep='\t',row.names=FALSE,quote=FALSE,na='')}
symbols<-keys(org.Hs.eg.db,keytype='SYMBOL');mp<-AnnotationDbi::select(org.Hs.eg.db,keys=symbols,keytype='SYMBOL',columns='ENTREZID');mp<-unique(mp[complete.cases(mp),])
write(mp,'symbol_entrez_mapping.tsv')
go<-AnnotationDbi::select(org.Hs.eg.db,keys=symbols,keytype='SYMBOL',columns=c('GOALL','ONTOLOGYALL'))
go<-unique(go[!is.na(go$GOALL)&go$ONTOLOGYALL=='BP',c('SYMBOL','GOALL')]);names(go)<-c('gene','term')
got<-AnnotationDbi::select(GO.db,keys=unique(go$term),columns='TERM',keytype='GOID');names(got)<-c('term','name')
paths<-list(GO_BP=split(go$gene,go$term)); namesets<-list(GO_BP=setNames(got$name,got$term))
kl<-read.delim(file.path(raw,'KEGG/link_hsa_pathway.tsv'),header=FALSE)
names(kl)<-c('term','entrez');kl$entrez<-sub('hsa:','',kl$entrez);kl$term<-sub('path:','',kl$term)
km<-merge(kl,mp,by.x='entrez',by.y='ENTREZID');paths$KEGG<-lapply(split(km$SYMBOL,km$term),unique)
kn<-read.delim(file.path(raw,'KEGG/list_pathway_hsa.tsv'),header=FALSE);namesets$KEGG<-setNames(kn$V2,sub('path:','',kn$V1))
rg<-strsplit(readLines(file.path(int,'Reactome_source.gmt'),warn=FALSE),'\t');paths$Reactome<-setNames(lapply(rg,function(x)unique(x[-c(1,2)])),sapply(rg,function(x)x[2]));namesets$Reactome<-setNames(sapply(rg,function(x)x[1]),sapply(rg,function(x)x[2]))
for(lib in names(paths)){paths[[lib]]<-lapply(paths[[lib]],function(x)sort(unique(x)));con<-file(file.path(int,paste0(lib,'.gmt')),'wt');for(t in names(paths[[lib]]))writeLines(paste(c(t,namesets[[lib]][t],paths[[lib]][[t]]),collapse='\t'),con);close(con)}
write(data.frame(package=c('R','fgsea','org.Hs.eg.db','GO.db','AnnotationDbi','Matrix','clusterProfiler'),version=c(as.character(getRversion()),sapply(c('fgsea','org.Hs.eg.db','GO.db','AnnotationDbi','Matrix','clusterProfiler'),function(p)as.character(packageVersion(p))))),'R_environment_versions.tsv')
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
meta<-read('G1_full_meta.tsv');meta<-meta[is.finite(meta$meta_z)&is.finite(meta$meta_logFC),];stopifnot(!anyDuplicated(meta$gene))
registry<-read('D0_interface_registry.tsv');members<-read('D0_interface_genes.tsv');links<-read('connections.tsv.gz');drugs<-read('all_drug_targets.tsv');bg<-read('disease_backgrounds.tsv');g1<-read('G1_244.tsv')
inputs<-list();details<-list();add<-function(id,genes,gate,kind,scope){inputs[[id]]<<-unique(as.character(genes));details[[id]]<<-data.frame(input=id,gate=gate,kind=kind,scope=scope)}
for(d in c('up','down')){x<-g1[g1$meta_logFC*(if(d=='up')1 else -1)>0,];add(paste0('G1_244_',d),x$gene,'G1','disease_ORA','robust human RNA association');x<-meta[meta$meta_FDR<cfg$fdr&meta$meta_logFC*(if(d=='up')1 else -1)>0,];add(paste0('G1_413_',d),x$gene,'G1','disease_ORA','significant sensitivity')}
for(i in seq_len(nrow(registry))){r<-registry[i,];if(r$gate!='G1')add(paste0(r$interface,'__membership'),members$gene[members$module==r$interface],r$gate,'program_annotation','PRESPECIFIED membership, not independent discovery')}
for(s in sort(unique(links$disease_set)))for(t in c('E1','E2','P','H')){x<-links[links$disease_set==s&links$tier==t,];if(nrow(x))add(paste0(s,'__',t),x$gene,substr(s,1,2),'drug_intersection',if(s=='G1_robust244'&&t=='E1')'G1 source reviewed; not formula efficacy' else 'exploratory source/program connection')}
det<-do.call(rbind,details);write(det,'ORA_input_registry.tsv');write(do.call(rbind,lapply(names(inputs),function(id)data.frame(input=id,gene=inputs[[id]]))),'ORA_input_genes.tsv')
ora<-list();oraqc<-list();mappingqc<-list();n<-0; cache<-new.env(parent=emptyenv())
for(lib in names(paths)){
  annotated<-sort(unique(unlist(paths[[lib]],use.names=FALSE)))
  for(gate in unique(det$gate)){
    background<-intersect(bg$gene[bg$gate==gate],annotated)
    for(mode in c('disease_measured','drug_record_selection_sensitivity')){
      ids<-det$input[det$gate==gate]
      if(mode!='disease_measured')ids<-ids[det$kind[match(ids,det$input)]=='drug_intersection']
      for(id in ids){
        universe<-background
        if(mode!='disease_measured'){tier<-sub('.*__','',id);universe<-intersect(universe,drugs$gene[drugs$tier==tier])}
        query<-intersect(inputs[[id]],universe);N<-length(universe);Q<-length(query)
        mappingqc[[length(mappingqc)+1]]<-data.frame(input=id,library=lib,background=mode,gene=inputs[[id]],in_background=inputs[[id]]%in%universe,annotated=inputs[[id]]%in%annotated)
        oraqc[[length(oraqc)+1]]<-data.frame(input=id,library=lib,background=mode,original_genes=length(inputs[[id]]),mapped_query=Q,universe=N,context=if(gate=='G2')'RNA proxy sensitivity; occupancy ORA unavailable' else det$scope[match(id,det$input)])
        if(N==0||Q==0)next
        cache_key<-paste(lib,gate,mode,if(mode=='disease_measured')'all' else tier,sep='::')
        if(exists(cache_key,envir=cache,inherits=FALSE)){psets<-get(cache_key,envir=cache)}else{psets<-lapply(paths[[lib]],intersect,y=universe);size<-lengths(psets);psets<-psets[size>=cfg$ora$min_size&size<=cfg$ora$max_size];assign(cache_key,psets,envir=cache)}
        if(!length(psets))next
        K<-lengths(psets);hits<-lapply(psets,intersect,y=query);k<-lengths(hits)
        p<-phyper(k-1,K,N-K,Q,lower.tail=FALSE)
        n<-n+1;ora[[n]]<-data.frame(input=id,library=lib,background=mode,term=names(psets),description=unname(namesets[[lib]][names(psets)]),query_size=Q,background_size=N,term_size=K,overlap=k,gene_ratio=k/Q,background_ratio=K/N,fold_enrichment=(k/Q)/(K/N),P=p,FDR_family=p.adjust(p,'BH'),genes=sapply(hits,paste,collapse=';'),interpretation=if(gate=='G2')'RNA_PROXY_ANNOTATION_ONLY' else if(det$kind[match(id,det$input)]=='program_annotation')'PRESELECTED_PROGRAM_ANNOTATION' else 'ASSOCIATION_ONLY',row.names=NULL)
      }
    }
  }
  cat('ORA complete',lib,'\n');flush.console()
}
ora<-do.call(rbind,ora);ora$FDR_global<-p.adjust(ora$P,'BH');write(ora,'ORA_all_tests.tsv.gz');write(ora[ora$FDR_family<cfg$fdr&ora$overlap>=cfg$ora$min_overlap_for_interpretation,],'ORA_significant_family.tsv');write(do.call(rbind,oraqc),'ORA_background_QC.tsv');write(do.call(rbind,mappingqc),'ORA_gene_mapping_QC.tsv')
gsea<-list();gq<-list()
for(rank in c('meta_z','meta_logFC'))for(size_mode in c('main','sensitivity')){
  if(rank=='meta_logFC'&&size_mode=='sensitivity')next
  bounds<-if(size_mode=='main')cfg$gsea$main_sizes else cfg$gsea$sensitivity_sizes
  stats<-setNames(meta[[rank]],meta$gene);stats<-stats[order(-stats,names(stats))]
  for(lib in names(paths)){
    set.seed(cfg$seed)
    result<-as.data.frame(fgseaMultilevel(pathways=paths[[lib]],stats=stats,minSize=bounds[1],maxSize=bounds[2],eps=cfg$gsea$eps,nPermSimple=cfg$gsea$nPermSimple,nproc=1))
    result$leadingEdge<-sapply(result$leadingEdge,paste,collapse=';');result$ranking<-rank;result$size_mode<-size_mode;result$library<-lib;result$description<-unname(namesets[[lib]][result$pathway]);result$FDR_family<-p.adjust(result$pval,'BH')
    gsea[[length(gsea)+1]]<-result
    gq[[length(gq)+1]]<-data.frame(ranking=rank,size_mode=size_mode,library=lib,input_genes=length(stats),annotated_input=sum(names(stats)%in%unique(unlist(paths[[lib]]))),tests=nrow(result),missing_P=sum(is.na(result$pval)),tied_rank_values=sum(duplicated(stats)),zero_rank=sum(stats==0))
    cat('GSEA',rank,size_mode,lib,nrow(result),'\n');flush.console()
  }
}
gsea<-do.call(rbind,gsea);gsea$FDR_all_libraries<-ave(gsea$pval,interaction(gsea$ranking,gsea$size_mode),FUN=function(x)p.adjust(x,'BH'));write(gsea,'GSEA_all_tests.tsv.gz');write(gsea[gsea$FDR_family<cfg$fdr,],'GSEA_significant_family.tsv');write(do.call(rbind,gq),'GSEA_QC.tsv')
write(data.frame(library=names(paths),terms=lengths(paths),genes=sapply(paths,function(x)length(unique(unlist(x))))),'pathway_library_summary.tsv')
cat('enrichment_finished\n')
