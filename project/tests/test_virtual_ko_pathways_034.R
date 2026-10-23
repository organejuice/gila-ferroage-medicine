suppressPackageStartupMessages(library(fgsea))
stats<-setNames(seq(500,1),paste0('G',seq_len(500)))
sets<-list(TOP=names(stats)[1:20],BOTTOM=names(stats)[481:500])
# The unweighted running sum reaches exactly 1 for an uninterrupted top set.
es<-calcGseaStat(stats,seq_len(20),gseaParam=0,scoreType='pos')
stopifnot(abs(es-1)<1e-12)
run<-function(){set.seed(34034);as.data.frame(fgseaSimple(sets,stats,minSize=15,maxSize=500,scoreType='pos',gseaParam=0,nperm=10000,nproc=1))}
a<-run();b<-run();stopifnot(isTRUE(all.equal(a,b)),a$pval[a$pathway=='TOP']<.05,a$pval[a$pathway=='BOTTOM']>.05)
cat('PASS: positive-tail enrichment, exact unweighted planted-set ES, deterministic fgsea replay.\n')
