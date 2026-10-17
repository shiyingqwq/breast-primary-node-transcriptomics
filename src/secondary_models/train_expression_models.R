source(file.path(Sys.getenv('BCPN_REPOSITORY'), 'src/secondary_models/runtime.R'))
dir.create(file.path(out,'work'),showWarnings=FALSE)
dev <- load_development();x <- dev$x;y <- dev$y;d <- dev$d
ref <- readRDS(file.path(Sys.getenv('BCPN_REPOSITORY'),'data/reference/ridge_reference.rds'))
stopifnot(identical(colnames(x),ref$feature_universe),length(ref$coefficients)==869)
seeds <- 20261006:20261010
folds <- lapply(seeds,function(s)folds_for(y,s))
old_f <- read.csv(file.path(model_inputs,'DEVELOPMENT_FOLDS.csv'))
stopifnot(identical(d$patient_id,old_f$patient_id),all(folds[[1]]==old_f$fold))
write.csv(do.call(rbind,lapply(seq_along(seeds),function(i)
  data.frame(patient_id=d$patient_id,seed=seeds[i],fold=folds[[i]]))),
  file.path(out,'work/REPEATED_DEVELOPMENT_FOLDS.csv'),row.names=FALSE)
alphas <- c(.25,.5,.75,1)
worker <- function(a) {
  logpath <- file.path(out,'work',paste0('alpha_',a,'.log'))
  results <- vector('list',length(seeds));summaries <- vector('list',length(seeds))
  for(i in seq_along(seeds)) {
    cat(format(Sys.time()),'START seed',seeds[i],'alpha',a,'\n',file=logpath,append=TRUE)
    cv <- cv.glmnet(x=x,y=y,family='binomial',alpha=a,type.measure='auc',
      nfolds=10,foldid=folds[[i]],keep=FALSE,standardize=TRUE,intercept=TRUE,parallel=FALSE)
    results[[i]] <- data.frame(alpha=a,seed=seeds[i],lambda=cv$lambda,
      CV_AUC=cv$cvm,CV_fold_SE=cv$cvsd,nonzero=cv$nzero)
    mi <- match(cv$lambda.min,cv$lambda);oi <- match(cv$lambda.1se,cv$lambda)
    summaries[[i]] <- data.frame(alpha=a,seed=seeds[i],lambda_min=cv$lambda.min,
      lambda_1se=cv$lambda.1se,CV_AUC_lambda_min=cv$cvm[mi],CV_AUC_lambda_1se=cv$cvm[oi],
      nonzero_lambda_min=cv$nzero[mi],nonzero_lambda_1se=cv$nzero[oi])
    cat(format(Sys.time()),'DONE seed',seeds[i],'best AUC',cv$cvm[mi],'\n',file=logpath,append=TRUE)
  }
  for(i in 2:length(results)) stopifnot(identical(results[[1]]$lambda,results[[i]]$lambda),
                                       identical(results[[1]]$nonzero,results[[i]]$nonzero))
  saveRDS(list(curves=results,summaries=summaries),
          file.path(out,'work',paste0('alpha_',a,'.rds')))
  TRUE
}
cat('Repeated EN CV: 20 fits, four independent alpha workers\n');flush.console()
cl <- parallel::makePSOCKcluster(4,outfile=file.path(out,'work/WORKERS.log'))
parallel::clusterExport(cl,c('x','y','folds','seeds','out','worker'))
libs <- .libPaths()
parallel::clusterCall(cl,function(paths){.libPaths(paths);library(glmnet);TRUE},libs)
status <- tryCatch(parallel::parLapply(cl,alphas,worker),finally=parallel::stopCluster(cl))
stopifnot(all(unlist(status)))
pieces <- lapply(alphas,function(a)readRDS(file.path(out,'work',paste0('alpha_',a,'.rds'))))
curves <- do.call(rbind,lapply(pieces,function(z)do.call(rbind,z$curves)))
repeat_summary <- do.call(rbind,lapply(pieces,function(z)do.call(rbind,z$summaries)))
write.csv(curves,file.path(out,'model_FULL_CV_RESULTS.csv'),row.names=FALSE)
write.csv(repeat_summary,file.path(out,'work/LAMBDA_MIN_1SE_BY_REPEAT.csv'),row.names=FALSE)
en <- do.call(rbind,lapply(seq_along(alphas),function(i) {
  z <- pieces[[i]]$curves
  means <- rowMeans(sapply(z,function(t)t$CV_AUC))
  sds <- apply(sapply(z,function(t)t$CV_AUC),1,sd)
  data.frame(alpha=alphas[i],lambda=z[[1]]$lambda,mean_CV_AUC=means,
    SD_CV_AUC=sds,nonzero=z[[1]]$nonzero,repeats=5)
}))
best <- max(en$mean_CV_AUC,na.rm=TRUE)
en$within_0_01_best <- is.finite(en$mean_CV_AUC)&en$mean_CV_AUC>=best-.01
eligible <- en[en$within_0_01_best&en$nonzero>0,]
eligible <- eligible[order(eligible$nonzero,-eligible$mean_CV_AUC,-eligible$alpha,-eligible$lambda),]
a <- eligible[1,]
en$selected <- en$alpha==a$alpha&en$lambda==a$lambda
write.csv(en,file.path(out,'model_ELASTIC_NET_CANDIDATES.csv'),row.names=FALSE)
final_a <- glmnet(x,y,family='binomial',alpha=a$alpha,standardize=TRUE,intercept=TRUE)
stopifnot(any(final_a$lambda==a$lambda))
ca <- as.matrix(coef(final_a,s=a$lambda))[,1];ba <- ca[-1];ba <- ba[ba!=0]
stopifnot(length(ba)==a$nonzero)
model_a <- list(name='Sparse_EN',method='Elastic Net',alpha=a$alpha,lambda=a$lambda,
  intercept=unname(ca[1]),coefficients=ba,genes=names(ba),gene_N=length(ba),
  development_CV_AUC=a$mean_CV_AUC,development_CV_SD=a$SD_CV_AUC,fit=final_a)
cat('Selected EN: alpha',a$alpha,'lambda',a$lambda,'genes',length(ba),'CV',a$mean_CV_AUC,'\n');flush.console()

ranked <- names(ref$coefficients)[order(-abs(ref$coefficients),names(ref$coefficients))]
write.csv(data.frame(gene=ranked,ridge_coefficient=ref$coefficients[ranked],rank=seq_along(ranked)),
          file.path(out,'work/REFERENCE_RIDGE_GENE_RANKING.csv'),row.names=FALSE)
safe_glm <- function(xx,yy) {
  warning_text <- character()
  f <- withCallingHandlers(glm.fit(cbind('(Intercept)'=1,xx),yy,family=binomial()),
    warning=function(w){warning_text <<- c(warning_text,conditionMessage(w));invokeRestart('muffleWarning')})
  ok <- isTRUE(f$converged)&&all(is.finite(f$coefficients))&&
    !any(grepl('did not converge|numerically 0 or 1',warning_text))
  list(fit=f,ok=ok,warnings=paste(warning_text,collapse='; '))
}
top_rows <- list();top_models <- list();kvalues <- c(10,20,30,50,75,100)
for(k in kvalues) {
  genes <- ranked[seq_len(k)];xx <- x[,genes,drop=FALSE]
  full <- safe_glm(xx,y);fold_auc <- numeric(10);oof <- rep(NA_real_,length(y))
  ok <- full$ok;issues <- full$warnings
  for(j in 1:10) {
    te <- folds[[1]]==j;g <- safe_glm(xx[!te,,drop=FALSE],y[!te]);ok <- ok&&g$ok
    issues <- paste(issues,g$warnings,sep='; ')
    if(g$ok) {
      oof[te] <- plogis(cbind(1,xx[te,,drop=FALSE])%*%g$fit$coefficients)
      fold_auc[j] <- auc_fixed(y[te],oof[te])
    } else fold_auc[j] <- NA_real_
  }
  if(ok) {
    beta <- full$fit$coefficients[-1];names(beta) <- genes
    m <- list(name='TopK',method='Ordinary logistic',alpha=NA_real_,lambda=NA_real_,
      intercept=unname(full$fit$coefficients[1]),coefficients=beta,genes=genes,gene_N=k,
      development_CV_AUC=mean(fold_auc),development_CV_SD=sd(fold_auc),fit=full$fit)
    fallback <- 'None'
  } else {
    cv <- cv.glmnet(xx,y,family='binomial',alpha=0,type.measure='auc',
                   foldid=folds[[1]],nfolds=10,standardize=TRUE,keep=TRUE)
    lam <- cv$lambda.min;ix <- match(lam,cv$lambda)
    cf <- as.matrix(coef(cv,s=lam))[,1];beta <- cf[-1];beta <- beta[beta!=0]
    m <- list(name='TopK',method='Panel-restricted ridge',alpha=0,lambda=lam,
      intercept=unname(cf[1]),coefficients=beta,genes=names(beta),gene_N=length(beta),
      development_CV_AUC=cv$cvm[ix],development_CV_SD=NA_real_,fit=cv$glmnet.fit)
    oof <- plogis(cv$fit.preval[,ix]);fallback <- paste('Ridge after logistic instability:',issues)
  }
  top_models[[as.character(k)]] <- m
  top_rows[[as.character(k)]] <- data.frame(K=k,nonzero=m$gene_N,method=m$method,
    CV_AUC=m$development_CV_AUC,CV_pooled_OOF_AUC=auc_fixed(y,oof),
    CV_fold_SD=m$development_CV_SD,alpha=m$alpha,lambda=m$lambda,fallback=fallback,
    selected=FALSE,CV_scope='Conditional on full-development ridge ranking; not nested')
  cat('TopK',k,'CV AUC',m$development_CV_AUC,'method',m$method,'\n');flush.console()
}
top <- do.call(rbind,top_rows)
choices <- top[top$CV_AUC>=max(top$CV_AUC)-.005,]
choices <- choices[order(choices$K,-choices$CV_AUC),]
b <- choices[1,];top$selected <- top$K==b$K
model_b <- top_models[[as.character(b$K)]]
write.csv(top,file.path(out,'model_TOPK_COMPRESSION_RESULTS.csv'),row.names=FALSE)
model_ref <- list(name='Ridge_reference',method='Full-development ridge reference',alpha=ref$alpha,
  lambda=ref$lambda,intercept=ref$intercept,coefficients=ref$coefficients,
  genes=ref$selected_genes,gene_N=length(ref$selected_genes),
  development_CV_AUC=ref$CV_fold_mean_AUC,development_CV_SD=NA_real_,fit=ref$molecular_model)
models <- list(Sparse_EN=model_a,TopK=model_b,Ridge_reference=model_ref)
selected <- list(models=models,feature_universe=colnames(x),seeds=seeds,
  EN_best_mean_CV_AUC=best,clinical_reference=ref$clinical_model,
  development_N=nrow(x),created_UTC=format(Sys.time(),tz='UTC',usetz=TRUE),
  METABRIC_used_for_new_candidate_fit_or_selection=FALSE,
  preprocessing='Use cohort-local gene z-scores with the specified gene universe and patient mask')
saveRDS(selected,file.path(out,'work/DEVELOPMENT_SELECTED_MODELS.rds'))
coefrows <- do.call(rbind,lapply(models,function(m)data.frame(model=m$name,
  gene=c('(Intercept)',names(m$coefficients)),coefficient=c(m$intercept,m$coefficients),
  alpha=m$alpha,lambda=m$lambda,type=c('intercept',rep('gene',length(m$coefficients))))))
write.csv(coefrows,file.path(out,'work/ALL_THREE_SELECTED_COEFFICIENTS.csv'),row.names=FALSE,na='')
write_json(list(development_N=nrow(x),EN=list(alpha=a$alpha,lambda=a$lambda,
  gene_N=model_a$gene_N,CV_AUC=model_a$development_CV_AUC,SD=model_a$development_CV_SD),
  TopK=list(K=b$K,method=model_b$method,CV_AUC=model_b$development_CV_AUC),
  reference_gene_N=model_ref$gene_N,external_used_for_new_candidate_selection=FALSE),
  file.path(out,'work/TRAINING_COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
capture.output(sessionInfo(),file=file.path(out,'work/TRAINING_SESSION.txt'))
cat('TRAINING COMPLETE: two compressed candidates selected before new external predictions\n')
