suppressPackageStartupMessages({library(limma);library(edgeR);library(jsonlite)})
run <- '2026-10-03-public-mechanism-analysis-026'
indir <- file.path('data/interim',run);out <- file.path('results',run)
cfg <- fromJSON('configs/public_mechanism_analysis_026.json',simplifyVector=FALSE)
args <- commandArgs(trailingOnly=TRUE)
datasets <- fromJSON(file.path(indir,if(length(args)) args[1] else 'dataset_manifest.json'),simplifyVector=FALSE)
if(length(args)>1) {out <- file.path(out,args[2]);dir.create(out,recursive=TRUE,showWarnings=FALSE)}
modules <- read.delim(file.path(indir,'frozen_modules.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
dir.create(file.path(out,'gene_effects'),recursive=TRUE,showWarnings=FALSE)
write_tsv <- function(x,p) write.table(x,p,sep='\t',quote=TRUE,row.names=FALSE,na='NA')
ols <- function(y,design,contrast) {
  fit <- lm.fit(design,y);df <- nrow(design)-qr(design)$rank
  beta <- sum(contrast*fit$coefficients)
  se <- sqrt(sum(fit$residuals^2)/df*as.numeric(t(contrast)%*%solve(crossprod(design))%*%contrast))
  c(effect=beta,SE=se,df=df,CI_low=beta-qt(.975,df)*se,CI_high=beta+qt(.975,df)*se,P=if(se>0) 2*pt(-abs(beta/se),df) else NA)
}
all_module <- list();all_loo <- list();qc <- list();summaries <- list();sample_qc <- list();score_tables <- list();i <- 0L
for(ds in datasets) {
  id <- ds$id;cat('Fitting',id,'\n');flush.console()
  input <- read.delim(gzfile(file.path(indir,paste0(id,'_expression.tsv.gz'))),check.names=FALSE,stringsAsFactors=FALSE)
  meta <- read.delim(file.path(indir,paste0(id,'_samples.tsv')),check.names=FALSE,stringsAsFactors=FALSE)
  maps <- read.delim(file.path(indir,paste0(id,'_gene_map.tsv')),check.names=FALSE,stringsAsFactors=FALSE)
  x <- as.matrix(input[,-1,drop=FALSE]);storage.mode(x) <- 'double';rownames(x) <- input$feature_id
  stopifnot(identical(colnames(x),meta$sample),all(is.finite(x)))
  levels_group <- sort(unique(meta$group))
  design <- outer(meta$group,levels_group,function(a,b) as.numeric(a==b));colnames(design) <- levels_group
  if(!is.null(ds$block_covariate) && isTRUE(ds$block_covariate)) {
    blevels <- sort(unique(meta$batch));bdesign <- outer(meta$batch,blevels[-1],function(a,b) as.numeric(a==b))
    colnames(bdesign) <- paste0('batch_',blevels[-1]);design <- cbind(design,bdesign)
  }
  stopifnot(qr(design)$rank==ncol(design))
  PCs <- if(ds$kind=='abundance') c(1,.5,2) else 1
  for(pc in PCs) {
    scenario <- if(pc==1) 'main' else paste0('pc',pc)
    if(ds$kind=='counts') {
      stopifnot(all(x>=0),all(abs(x-round(x))<1e-7))
      dge <- DGEList(x)
      keep <- filterByExpr(dge,design=design,min.count=cfg$counts$min_count,min.total.count=cfg$counts$min_total_count,large.n=cfg$counts$large_n,min.prop=cfg$counts$min_prop)
      dge <- calcNormFactors(dge[keep,,keep.lib.sizes=FALSE],method='TMM')
      v <- voom(dge,design,plot=FALSE);E <- v$E;fit <- lmFit(v,design)
      fit <- eBayes(fit,trend=FALSE)
    } else {
      transformed <- ds$kind %in% c('microarray_log2','rlog_normalized','protein_log2')
      keep <- apply(x,1,function(z) sd(z)>0 && (transformed || any(z>0)))
      E <- if(transformed) x[keep,,drop=FALSE] else log2(x[keep,,drop=FALSE]+pc)
      fit <- eBayes(lmFit(E,design),trend=TRUE)
    }
    gene_symbols <- maps$human_symbol[match(rownames(E),maps$feature_id)];gene_symbols[is.na(gene_symbols)] <- ''
    valid <- gene_symbols!=''
    # Exact duplicate gene mappings collapsed for scores only, not DE measurements.
    G <- rowsum(E[valid,,drop=FALSE],gene_symbols[valid],reorder=TRUE)/as.numeric(table(gene_symbols[valid]))
    gsd <- apply(G,1,sd);G <- G[is.finite(gsd)&gsd>0,,drop=FALSE]
    Z <- sweep(sweep(G,1,rowMeans(G)),1,apply(G,1,sd),'/')
    scores <- list();coverage <- list()
    for(mod in unique(modules$module)) {
      definition <- modules[modules$module==mod,];genes <- intersect(definition$gene,rownames(Z))
      coverage[[mod]] <- c(available=length(genes),total=nrow(definition))
      if(length(genes)>=cfg$minimum_score_genes) scores[[mod]] <- colMeans(Z[genes,,drop=FALSE]*definition$weight[match(genes,definition$gene)])
    }
    S <- do.call(rbind,scores)
    if(scenario=='main') {
      qc[[id]] <- data.frame(dataset=id,n_samples=ncol(E),n_input_features=nrow(x),n_tested_features=nrow(E),n_mapped_symbols=nrow(G),model=ds$kind,inferential=ds$formal,limitation=ds$limitation)
      pca <- prcomp(t(E),center=TRUE,scale.=FALSE)$x
      sample_qc[[id]] <- data.frame(dataset=id,sample=colnames(E),group=meta$group,PC1=pca[,1],PC2=pca[,2],library_sum=colSums(x),zero_fraction=colMeans(x==0))
      score_tables[[id]] <- data.frame(dataset=id,sample=rep(colnames(S),each=nrow(S)),module=rep(rownames(S),ncol(S)),score=as.vector(S))
    }
    contrasts <- ds$contrasts
    if(id=='GSE325931') contrasts <- list('iron_RSL3-tri_vehicle','co_vehicle-mono_vehicle','tri_vehicle-mono_vehicle')
    for(ct in contrasts) {
      C <- makeContrasts(contrasts=ct,levels=design)[,1]
      if(isTRUE(ds$formal)) {
        cf <- contrasts.fit(fit,C);cf <- eBayes(cf,trend=ds$kind!='counts')
        tab <- topTable(cf,number=Inf,sort.by='none',confint=TRUE)
        tab$feature_id <- rownames(tab);tab$human_symbol <- gene_symbols;tab$dataset <- id;tab$contrast <- ct
        if(scenario=='main') {
          dest <- file.path(out,'gene_effects',paste0(id,'_',gsub('[^A-Za-z0-9]','_',ct),'.tsv.gz'))
          con <- gzfile(dest,'wt',encoding='UTF-8');write_tsv(tab,con);close(con)
          summaries[[paste(id,ct)]] <- data.frame(dataset=id,contrast=ct,features=nrow(tab),FDR_lt_005=sum(tab$adj.P.Val<cfg$fdr,na.rm=TRUE),up=sum(tab$adj.P.Val<cfg$fdr & tab$logFC>0,na.rm=TRUE),down=sum(tab$adj.P.Val<cfg$fdr & tab$logFC<0,na.rm=TRUE),scope='within_dataset_replicate_expression; not human AD therapeutic validation')
        }
      }
      for(mod in rownames(S)) {
        i <- i+1L;y <- S[mod,]
        if(isTRUE(ds$formal)) est <- ols(y,design,C) else {
          beta <- lm.fit(design,y)$coefficients;est <- c(effect=sum(C*beta),SE=NA,df=NA,CI_low=NA,CI_high=NA,P=NA)
        }
        all_module[[i]] <- data.frame(dataset=id,contrast=ct,module=mod,scenario=scenario,t(est),genes_available=coverage[[mod]]['available'],genes_total=coverage[[mod]]['total'],inferential=ds$formal,limitation=ds$limitation,row.names=NULL)
        if(scenario=='main' && isTRUE(ds$formal)) {
          loo <- sapply(seq_along(y),function(j) {
            dd <- design[-j,,drop=FALSE]
            if(qr(dd)$rank<ncol(dd) || nrow(dd)<=ncol(dd)) return(NA_real_)
            sum(C*lm.fit(dd,y[-j])$coefficients)
          })
          all_loo[[length(all_loo)+1]] <- data.frame(dataset=id,contrast=ct,module=mod,sample_left_out=meta$sample,effect=unname(loo),full_effect=unname(est['effect']),row.names=NULL)
        }
      }
    }
  }
}
modtab <- do.call(rbind,all_module);modtab$FDR <- NA_real_
for(sc in unique(modtab$scenario)) {
  use <- which(modtab$scenario==sc & is.finite(modtab$P));modtab$FDR[use] <- p.adjust(modtab$P[use],method='BH')
}
write_tsv(modtab,file.path(out,'module_effects.tsv'));write_tsv(do.call(rbind,all_loo),file.path(out,'module_leave_one_out.tsv'))
write_tsv(do.call(rbind,qc),file.path(out,'fit_qc.tsv'));write_tsv(do.call(rbind,summaries),file.path(out,'gene_contrast_summary.tsv'))
write_tsv(do.call(rbind,sample_qc),file.path(out,'sample_qc.tsv'));write_tsv(do.call(rbind,score_tables),file.path(out,'sample_module_scores.tsv'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_environment_lock.txt'),useBytes=TRUE)
writeLines(capture.output(warnings()),file.path(out,'R_warnings.txt'),useBytes=TRUE)
cat('Completed',length(datasets),'datasets',nrow(modtab),'module records\n')
