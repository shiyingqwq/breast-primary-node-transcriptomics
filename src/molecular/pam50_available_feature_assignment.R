source(file.path(Sys.getenv('BCPN_REPOSITORY'), 'src/molecular/runtime.R'))
stopifnot(!file.exists(file.path(out,'work/COMPLETION.json')))
rescale <- get('rescale',asNamespace('genefu'));intrinsic.cluster.predict <- genefu::intrinsic.cluster.predict
d <- readRDS(file.path(out,'work/ANALYSIS_RESULTS.rds'));input <- readRDS(file.path(root,'enrichment/work/INPUTS.rds'))
p <- readRDS(file.path(out,'work/PAM50_INPUT_MAPPING.rds'));use <- p$details$available
stopifnot(sum(use)==43)
xx <- t(input$x[p$index[use],,drop=FALSE]);colnames(xx) <- p$details$centroid_probe[use]
ann <- data.frame(probe=colnames(xx),EntrezGene.ID=p$details$EntrezGene.ID[use],row.names=colnames(xx))
if(file.exists(file.path(out,'work/PAM50_PREDICTION.rds')))pred <- readRDS(file.path(out,'work/PAM50_PREDICTION.rds')) else
 pred <- intrinsic.cluster.predict(sbt.model=p$model,data=xx,annot=ann,do.mapping=FALSE,do.prediction.strength=FALSE,verbose=TRUE)
labels <- c(Basal='Basal-like',Her2='HER2-enriched',LumA='Luminal A',LumB='Luminal B',Normal='Normal-like')
pam <- data.frame(d$samples,PAM50_subtype=unname(labels[as.character(pred$subtype[d$samples$GSM])]),
 status='CLASSIFIED_RESEARCH_USE_43_OF_50_AVAILABLE_FEATURE_RESEARCH_IMPLEMENTATION',available_genes=43,
 implementation='Released genefu 2.40.0 pam50.robust source; incomplete-gene research implementation, not clinical assay')
pam <- cbind(pam,as.data.frame(pred$cor[d$samples$GSM,,drop=FALSE]));pam$correlation_margin <- apply(pred$cor[d$samples$GSM,,drop=FALSE],1,function(v)diff(sort(v,decreasing=TRUE)[2:1]))
stopifnot(!anyNA(pam$PAM50_subtype),all(is.finite(as.matrix(pred$cor))),all(pam$correlation_margin>=0))
write_csv(pam,'molecular_PAM50_SAMPLE_RESULTS.csv');saveRDS(pred,file.path(out,'work/PAM50_PREDICTION.rds'))
transformed <- apply(xx,2,function(v)(rescale(v,q=p$model$rescale.q,na.rm=TRUE)-.5)*2)
quantiles <- t(apply(xx,2,quantile,probs=c(.025,.975)))
write_csv(data.frame(centroid_gene=rownames(quantiles),lower_2_5_pct=quantiles[,1],upper_97_5_pct=quantiles[,2]),'work/PAM50_POOLED_TRANSFORMATION.csv')
manual_cor <- cor(t(transformed),p$model$centroids[use,,drop=FALSE],method=p$model$method.cor)
stopifnot(identical(dim(manual_cor),dim(pred$cor)),max(abs(manual_cor-pred$cor))<1e-12)
pd <- d$discordance;pd$Primary_PAM50 <- pam$PAM50_subtype[match(pd$Primary_GSM,pam$GSM)];pd$Node_PAM50 <- pam$PAM50_subtype[match(pd$Node_GSM,pam$GSM)]
pd$evaluable <- TRUE;pd$PAM50_discordant <- pd$Primary_PAM50!=pd$Node_PAM50;pd$status <- 'EVALUABLE_INCOMPLETE_GENE_RESEARCH_CLASSIFICATION'
write_csv(pd,'molecular_PAM50_DISCORDANCE.csv')
lv <- c('Luminal A','Luminal B','HER2-enriched','Basal-like','Normal-like')
tm <- table(factor(pd$Primary_PAM50,levels=lv),factor(pd$Node_PAM50,levels=lv));transition <- as.data.frame(tm)
names(transition) <- c('Primary_subtype','Node_subtype','N');transition$N <- as.integer(transition$N);transition$evaluable_pairs <- 36
transition$status <- 'EVALUABLE_INCOMPLETE_GENE_RESEARCH_CLASSIFICATION';write_csv(transition,'molecular_PAM50_TRANSITION_MATRIX.csv')
n <- nrow(pd);discord <- sum(pd$PAM50_discordant);pe <- sum(rowSums(tm)*colSums(tm))/n^2;po <- sum(diag(tm))/n
kappa <- if(pe<1)(po-pe)/(1-pe) else NA_real_;pr <- discord/n;z <- 1.96;center <- (pr+z^2/(2*n))/(1+z^2/n);half <- z*sqrt(pr*(1-pr)/n+z^2/(4*n^2))/(1+z^2/n)
save_json(list(coverage=43,total_pairs=36,classified_pairs=n,concordant=n-discord,discordant=discord,discordance_pct=100*pr,
 Cohen_kappa=kappa,discordance_Wilson95_low=center-half,discordance_Wilson95_high=center+half,
 status='43_OF_50_GENE_RESEARCH_FALLBACK'), 'work/PAM50_SUMMARY.json')
drivers <- c('StromalScore','B_lineage','Fibroblast','Endothelial','EMT','E2F','G2M')
driver <- do.call(rbind,lapply(drivers,function(s){
 a <- abs(d$delta[s,pd$PAM50_discordant]);b <- abs(d$delta[s,!pd$PAM50_discordant]);a <- a[is.finite(a)];b <- b[is.finite(b)]
 test <- if(length(a)&&length(b))wilcox.test(a,b,exact=FALSE,correct=TRUE) else list(p.value=NA_real_,statistic=NA_real_)
 data.frame(score=s,discordant_N=length(a),concordant_N=length(b),discordant_median_abs_delta=if(length(a))median(a) else NA_real_,
  concordant_median_abs_delta=if(length(b))median(b) else NA_real_,raw_P=test$p.value,status=if(length(a)&&length(b))'EVALUABLE' else 'NOT_RUN_ONE_GROUP_EMPTY',
  test='Two-sided rank-sum, exact=FALSE, continuity correction')
}))
driver$BH_q <- p.adjust(driver$raw_P,'BH');driver$BH_family <- 'Seven prespecified discordance drivers';write_csv(driver,'molecular_DISCORDANCE_DRIVERS.csv')
master <- read_csv(file.path(out,'molecular_PATIENT_PAIRED_SIGNATURES.csv'))
master$Primary_PAM50 <- pd$Primary_PAM50;master$Node_PAM50 <- pd$Node_PAM50;master$evaluable <- pd$evaluable
master$PAM50_discordant <- pd$PAM50_discordant;master$status <- pd$status;write_csv(master,'molecular_PATIENT_PAIRED_SIGNATURES.csv')
d$pam <- pam;d$discordance <- pd;d$transition <- transition;d$drivers <- driver
saveRDS(d,file.path(out,'work/ANALYSIS_RESULTS.rds'))
q <- fromJSON(file.path(out,'work/ANALYSIS_COMPLETION.json'));q$PAM50_classified_pairs <- 36;q$discordance_driver_tests <- sum(is.finite(driver$raw_P))
q$PAM50_research_fallback_coverage <- 43;save_json(q,'work/ANALYSIS_COMPLETION.json')
cat('Research PAM50: discordant',discord,'/36; kappa',kappa,'\n')
