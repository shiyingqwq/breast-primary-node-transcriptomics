options(warn=1)
root <- normalizePath(Sys.getenv('BCPN_WORKSPACE'), winslash='/')
model_inputs <- file.path(root, 'models/inputs')
out <- file.path(root, 'models')
if(nzchar(Sys.getenv("BCPN_R_LIBRARY"))) .libPaths(c(Sys.getenv("BCPN_R_LIBRARY"),.libPaths()))
library(glmnet)
library(pROC)
library(jsonlite)
auc_fixed <- function(y, p) as.numeric(auc(roc(y,p,levels=c(0,1),direction='<',quiet=TRUE)))
folds_for <- function(y, seed) {
  set.seed(seed)
  f <- integer(length(y))
  for(v in c(0,1)) {
    ix <- which(y==v)
    f[ix] <- sample(rep(seq_len(10),length.out=length(ix)))
  }
  f
}
load_development <- function() {
  e <- read.csv(file.path(model_inputs,'GSE96058_EXPRESSION_Z.csv'),check.names=FALSE)
  d <- read.csv(file.path(model_inputs,'GSE96058_CLINICAL_TRAINING.csv'),stringsAsFactors=FALSE)
  d <- d[match(e$patient_id,d$patient_id),]
  stopifnot(identical(e$patient_id,d$patient_id),!anyDuplicated(d$patient_id))
  ok <- !is.na(d$node_positive)
  x <- as.matrix(e[ok,-1]);d <- d[ok,]
  stopifnot(nrow(x)==3109,ncol(x)==869,all(is.finite(x)))
  list(x=x,d=d,y=d$node_positive)
}
clinical_coding <- function(d) {
  d$grade <- factor(d$grade,levels=c(1,2,3))
  d$ER <- factor(d$ER,levels=c(0,1))
  d
}
lp_for <- function(m,x) {
  stopifnot(all(names(m$coefficients)%in%colnames(x)))
  as.numeric(m$intercept + x[,names(m$coefficients),drop=FALSE] %*% m$coefficients)
}
metric <- function(label,y,p) {
  stopifnot(length(y)==length(p),all(is.finite(p)),all(y%in%c(0,1)))
  r <- roc(y,p,levels=c(0,1),direction='<',quiet=TRUE)
  set.seed(20261006)
  bounds <- as.numeric(ci.auc(r,method='bootstrap',boot.n=2000,boot.stratified=TRUE))
  data.frame(model=label,N=length(y),node_positive_N=sum(y),AUC=as.numeric(auc(r)),
    AUC_CI_low=bounds[1],AUC_CI_high=bounds[3],Brier=mean((y-p)^2),
    mean_probability=mean(p),observed_node_fraction=mean(y),
    CI_method='Stratified patient bootstrap 2000; percentile; fixed predictions')
}

dir.create(file.path(out,'work'),recursive=TRUE,showWarnings=FALSE)
