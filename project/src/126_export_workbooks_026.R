options(warn=1)
suppressPackageStartupMessages(library(readxl))
run <- '2026-10-03-public-mechanism-analysis-026'
out <- file.path('data/interim',run,'workbooks')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
files <- c(list.files(out,pattern='[.]xls$',full.names=TRUE),list.files(file.path('data/raw',run),pattern='[.]xlsx$',full.names=TRUE))
for (p in files) {
  staging <- tempfile(fileext=paste0('.',tools::file_ext(p)))
  if (!file.copy(p,staging,overwrite=FALSE)) stop('Cannot stage workbook')
  con <- file(staging,'rb'); magic <- readBin(con,'raw',8); close(con)
  is_excel <- identical(magic[1:4],as.raw(c(208,207,17,224))) || identical(magic[1:2],as.raw(c(80,75)))
  sheets <- if (is_excel) excel_sheets(staging) else 'tab_delimited_despite_xls_suffix'
  for (i in seq_along(sheets)) {
    x <- if(is_excel) as.data.frame(read_excel(staging,sheet=i,.name_repair='minimal',col_types='text')) else read.delim(staging,check.names=FALSE,quote='',comment.char='',fileEncoding='UTF-8')
    dest <- file.path(out,paste0(basename(p),'_sheet',i,'.tsv'))
    write.table(x,dest,sep='\t',row.names=FALSE,quote=TRUE,na='NA',fileEncoding='UTF-8')
    cat(basename(p),'SHEET',i,sheets[i],nrow(x),ncol(x),'\n')
    print(names(x)[seq_len(min(ncol(x),22))]);print(utils::head(x[,seq_len(min(ncol(x),6)),drop=FALSE],2))
  }
  unlink(staging)
}
