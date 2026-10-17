source(file.path(Sys.getenv('BCPN_REPOSITORY'), 'src/molecular/runtime.R'))
out <- file.path(root,'clinical_context')
read_csv <- function(p)read.csv(p,check.names=FALSE,fileEncoding='UTF-8-BOM')
write_csv <- function(d,n)write.csv(d,file.path(out,n),row.names=FALSE,na='')
save_json <- function(d,n)write_json(d,file.path(out,n),pretty=TRUE,auto_unbox=TRUE,digits=NA)

dir.create(file.path(out,'work'),recursive=TRUE,showWarnings=FALSE)
