options(warn=1)
root <- normalizePath(Sys.getenv('BCPN_WORKSPACE'),winslash='/')
ev <- file.path(root,'clinical')
if(nzchar(Sys.getenv("BCPN_R_LIBRARY"))) .libPaths(c(Sys.getenv("BCPN_R_LIBRARY"),.libPaths()))
library(singscore);library(jsonlite)
stopifnot(as.character(packageVersion('singscore'))=='1.28.1')
fg <- read.csv(file.path(root,'cohort_NODE_LIKE_SCORE_GENES.csv'))
up <- fg$clinical_symbol[fg$direction=='Node_up'];dn <- fg$clinical_symbol[fg$direction=='Node_down']
stopifnot(length(up)==315,length(dn)==554)
config <- fromJSON(file.path(ev,'SCORE_INPUT_CONFIG.json'),simplifyVector=FALSE)
all <- list();constants <- list();summaries <- list()
ledger<-read.csv(file.path(root,'clinical_TECHNICAL_SCORE_EXCLUSION_LEDGER.csv'),stringsAsFactors=FALSE)
for (item in config) {
 cohort <- item$cohort;ids<-unlist(item$patient_ids);bg<-unlist(item$background_genes)
 stopifnot(length(bg)==19349)
 conn<-file(file.path(root,item$binary_file),'rb');scores<-numeric(length(ids))
 for(start in seq(1,length(ids),by=128)) {
  end<-min(start+127,length(ids));n<-end-start+1
  expr<-matrix(readBin(conn,'double',n=n*19349,size=8,endian='little'),nrow=19349)
  stopifnot(ncol(expr)==n,all(is.finite(expr)))
  rownames(expr)<-bg;colnames(expr)<-ids[start:end]
  ranked<-rankGenes(expr,tiesMethod='min',stableGenes=NULL)
  ss<-simpleScore(ranked,upSet=up,downSet=dn,subSamples=NULL,centerScore=TRUE,dispersionFun=stats::mad,knownDirection=TRUE)
  stopifnot(identical(rownames(ss),colnames(expr)))
  scores[start:end]<-ss$TotalScore
  cat(cohort,'scored',end,'/',length(ids),'\n');flush.console()
 }
 close(conn);stopifnot(all(is.finite(scores)),sd(scores)>0)
 mu<-mean(scores);sigma<-sd(scores)
 lk<-ledger[ledger$cohort==cohort & tolower(as.character(ledger$score_eligible))=='true',];lk<-lk[match(ids,lk$expression_label),]
 stopifnot(!anyNA(lk$patient_id),length(ids)==if(cohort=='GSE96058')3273 else 1964)
 all[[cohort]]<-data.frame(cohort=cohort,patient_id=lk$patient_id,sample_id=lk$sample_id,raw_NodeLikeScore=scores,standardized_NodeLikeScore=(scores-mu)/sigma,score_eligible=TRUE,technical_QC_pass=TRUE)
 constants[[cohort]]<-data.frame(cohort=cohort,N=length(ids),mean=mu,SD=sigma,SD_definition='sample SD; n-1')
 summaries[[cohort]]<-data.frame(cohort=cohort,N=length(ids),missing_score_N=0,min=min(scores),max=max(scores),mean=mu,SD=sigma,median=median(scores),Q1=unname(quantile(scores,.25)),Q3=unname(quantile(scores,.75)),IQR=IQR(scores))
}
stopifnot(!file.exists(file.path(root,'clinical_NODE_LIKE_SCORE_VALUES_PREOUTCOME.csv')))
write.csv(do.call(rbind,all),file.path(root,'clinical_NODE_LIKE_SCORE_VALUES_PREOUTCOME.csv'),row.names=FALSE)
write.csv(do.call(rbind,constants),file.path(root,'clinical_SCORE_STANDARDIZATION_CONSTANTS.csv'),row.names=FALSE)
write.csv(do.call(rbind,summaries),file.path(root,'clinical_SCORE_TECHNICAL_QC.csv'),row.names=FALSE)
