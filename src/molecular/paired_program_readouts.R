source(file.path(Sys.getenv('BCPN_REPOSITORY'), 'src/molecular/runtime.R'))
stopifnot(!file.exists(file.path(out,'work/ANALYSIS_RESULTS.rds')))
input <- readRDS(file.path(root,'enrichment/work/INPUTS.rds'));samples <- input$samples;x <- input$x
tme <- readRDS(file.path(root,'enrichment/work/TME_RESULTS.rds'));pm <- tme$pair_map
stopifnot(nrow(pm)==36,!anyDuplicated(pm$patient),identical(colnames(x),samples$GSM),
 all(samples$patient[match(pm$Primary,samples$GSM)]==pm$patient),all(samples$patient[match(pm$Node,samples$GSM)]==pm$patient))
qc <- fromJSON(file.path(out,'work/PAM50_QC.json'));prepared <- readRDS(file.path(out,'work/PAM50_INPUT_MAPPING.rds'))
pam <- data.frame(samples,PAM50_subtype=NA_character_,status=paste0('SKIP_COVERAGE_',qc$available,'_OF_50'),
 available_genes=qc$available,implementation='Research implementation; not clinical assay')
if(qc$available==43){
 rescale <- get('rescale',asNamespace('genefu'));intrinsic.cluster.predict <- genefu::intrinsic.cluster.predict
 use <- prepared$details$available;xx <- t(x[prepared$index[use],,drop=FALSE]);colnames(xx) <- prepared$details$centroid_probe[use]
 ann <- data.frame(probe=colnames(xx),EntrezGene.ID=prepared$details$EntrezGene.ID[use],row.names=colnames(xx))
 if(requireNamespace('genefu',quietly=TRUE))pred <- genefu::intrinsic.cluster.predict(prepared$model,xx,ann,do.mapping=FALSE,do.prediction.strength=FALSE) else
  pred <- intrinsic.cluster.predict(prepared$model,xx,ann,do.mapping=FALSE,do.prediction.strength=FALSE)
 labels <- c(Basal='Basal-like',Her2='HER2-enriched',LumA='Luminal A',LumB='Luminal B',Normal='Normal-like')
 pam$PAM50_subtype <- unname(labels[as.character(pred$subtype[samples$GSM])]);pam$status <- ifelse(is.na(pam$PAM50_subtype),'NOT_EVALUABLE','CLASSIFIED_RESEARCH_USE')
 pam <- cbind(pam,as.data.frame(pred$cor[samples$GSM,,drop=FALSE]));pam$correlation_margin <- apply(pred$cor[samples$GSM,,drop=FALSE],1,function(v)diff(sort(v,decreasing=TRUE)[2:1]))
 saveRDS(pred,file.path(out,'work/PAM50_PREDICTION.rds'))
}
write_csv(pam,'molecular_PAM50_SAMPLE_RESULTS.csv')
pd <- data.frame(patient=pm$patient,Primary_GSM=pm$Primary,Node_GSM=pm$Node,
 Primary_PAM50=pam$PAM50_subtype[match(pm$Primary,pam$GSM)],Node_PAM50=pam$PAM50_subtype[match(pm$Node,pam$GSM)])
pd$evaluable <- !is.na(pd$Primary_PAM50)&!is.na(pd$Node_PAM50);pd$PAM50_discordant <- ifelse(pd$evaluable,pd$Primary_PAM50!=pd$Node_PAM50,NA)
pd$status <- ifelse(pd$evaluable,'EVALUABLE',paste0('PAM50_SKIPPED_',qc$available,'_OF_50'))
write_csv(pd,'molecular_PAM50_DISCORDANCE.csv')
levels <- c('Luminal A','Luminal B','HER2-enriched','Basal-like','Normal-like')
tm <- table(factor(pd$Primary_PAM50[pd$evaluable],levels=levels),factor(pd$Node_PAM50[pd$evaluable],levels=levels))
transition <- as.data.frame(tm);names(transition) <- c('Primary_subtype','Node_subtype','N')
transition$N <- if(sum(pd$evaluable))as.integer(transition$N) else NA_integer_
transition$evaluable_pairs <- sum(pd$evaluable);transition$status <- if(sum(pd$evaluable))'EVALUABLE' else 'NOT_EVALUABLE_NO_CLASSIFICATIONS'
write_csv(transition,'molecular_PAM50_TRANSITION_MATRIX.csv')
n <- sum(pd$evaluable);discord <- if(n)sum(pd$PAM50_discordant[pd$evaluable]) else NA_integer_;kappa <- NA_real_
wilson_low <- wilson_high <- NA_real_
if(n){pe <- sum(rowSums(tm)*colSums(tm))/n^2;po <- sum(diag(tm))/n;kappa <- if(pe<1)(po-pe)/(1-pe) else NA_real_
 p <- discord/n;z <- 1.96;center <- (p+z^2/(2*n))/(1+z^2/n);half <- z*sqrt(p*(1-p)/n+z^2/(4*n^2))/(1+z^2/n)
 wilson_low <- center-half;wilson_high <- center+half}
save_json(list(coverage=qc$available,total_pairs=36,classified_pairs=n,concordant=if(n)n-discord else NA_integer_,
 discordant=discord,discordance_pct=if(n)100*discord/n else NA_real_,Cohen_kappa=kappa,
 discordance_Wilson95_low=wilson_low,discordance_Wilson95_high=wilson_high,
 status=if(n)'RESEARCH_CLASSIFICATION' else 'PAM50_SKIPPED_INSUFFICIENT_SELECTED_COVERAGE'), 'work/PAM50_SUMMARY.json')

read_gmt <- function(p){a <- strsplit(readLines(p),'\t');setNames(lapply(a,function(v)unique(v[-c(1,2)])),vapply(a,`[`,'',1))}
h <- read_gmt(file.path(root,'resources/h.all.v2025.1.Hs.symbols.gmt'))
go <- read_gmt(file.path(root,'enrichment/sources/GO_BP.v2025.1.Hs.gmt'));re <- read_gmt(file.path(root,'enrichment/sources/REACTOME.v2025.1.Hs.gmt'))
sets <- list(ER_EARLY=h[['HALLMARK_ESTROGEN_RESPONSE_EARLY']],ER_LATE=h[['HALLMARK_ESTROGEN_RESPONSE_LATE']],
 E2F=h[['HALLMARK_E2F_TARGETS']],G2M=h[['HALLMARK_G2M_CHECKPOINT']],EMT=h[['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']],
 TNFA_NFKB=h[['HALLMARK_TNFA_SIGNALING_VIA_NFKB']],Hypoxia=h[['HALLMARK_HYPOXIA']],Adipogenesis=h[['HALLMARK_ADIPOGENESIS']],
 BCR=go[['GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY']],ECM=re[['REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION']],
 ERBB2_SIGNALING=re[['REACTOME_SIGNALING_BY_ERBB2']])
stopifnot(length(sets)==11,all(lengths(sets)>0))
map <- input$mapping;use <- map$annotation_usable
xt <- x[use,,drop=FALSE];rownames(xt) <- map$approved_symbol[use];stopifnot(!anyDuplicated(rownames(xt)))
means <- rowMeans(xt);sds <- apply(xt,1,sd);zz <- (xt-means)/sds
constant <- !is.finite(sds)|sds<=0;zz[constant,] <- NA_real_
scores <- matrix(NA_real_,11,72,dimnames=list(names(sets),samples$GSM));coverage <- genedetails <- list()
for(i in seq_along(sets)){
 score <- names(sets)[i];members <- sets[[i]];mapped <- intersect(members,rownames(xt));used <- mapped[!constant[match(mapped,rownames(xt))]]
 if(length(used))scores[i,] <- colMeans(zz[used,,drop=FALSE])
 coverage[[i]] <- data.frame(score=score,set_size=length(members),mapped_N=length(mapped),used_N=length(used),coverage_pct=100*length(used)/length(members),
  missing_genes=paste(setdiff(members,mapped),collapse=';'),constant_genes=paste(setdiff(mapped,used),collapse=';'),
  full_members=paste(members,collapse=';'),method='Mean across available gene z-scores; pooled 72 sample mean and sample SD; no DEG subset')
 genedetails[[i]] <- data.frame(score=score,gene=members,available=members%in%used,
  pooled_gene_mean=means[match(members,rownames(xt))],pooled_gene_SD=sds[match(members,rownames(xt))],
  selected_original_symbol=map$gene_symbol[match(members,map$approved_symbol)],row.names=NULL)
}
write_csv(do.call(rbind,coverage),'work/SIGNATURE_COVERAGE.csv');write_csv(do.call(rbind,genedetails),'work/SIGNATURE_GENE_STANDARDIZATION.csv')
write_csv(data.frame(samples,t(scores),check.names=FALSE),'molecular_SIGNATURE_SCORES.csv')
genes <- c('ESR1','PGR','ERBB2','MKI67');rx <- matrix(NA_real_,4,72,dimnames=list(genes,samples$GSM))
for(g in genes)if(g%in%rownames(x))rx[g,] <- x[g,]
write_csv(data.frame(samples,t(rx),check.names=FALSE),'molecular_RECEPTOR_EXPRESSION.csv')
ts <- rbind(StromalScore=tme$ESTIMATE$scores['StromalScore',samples$GSM],ImmuneScore=tme$ESTIMATE$scores['ImmuneScore',samples$GSM],
 B_lineage=tme$MCP_COUNTER$scores['B lineage',samples$GSM],T_cells=tme$MCP_COUNTER$scores['T cells',samples$GSM],
 Fibroblast=tme$MCP_COUNTER$scores['Fibroblasts',samples$GSM],Endothelial=tme$MCP_COUNTER$scores['Endothelial cells',samples$GSM],
 Monocytic=tme$MCP_COUNTER$scores['Monocytic lineage',samples$GSM],Neutrophils=tme$MCP_COUNTER$scores['Neutrophils',samples$GSM])
readouts <- rbind(rx,scores,ts);stopifnot(nrow(readouts)==23,ncol(readouts)==72)
primary <- readouts[,pm$Primary,drop=FALSE];node <- readouts[,pm$Node,drop=FALSE];delta <- node-primary
colnames(primary) <- colnames(node) <- colnames(delta) <- pm$patient
units <- c(rep('author-processed log2 expression',4),rep('mean pooled gene z-score',11),rep('original ESTIMATE score',2),rep('original MCP-counter score',6))
tests <- do.call(rbind,lapply(seq_len(nrow(readouts)),function(i){
 a <- primary[i,];b <- node[i,];d <- delta[i,];keep <- is.finite(a)&is.finite(b);a <- a[keep];b <- b[keep];d <- d[keep]
 test <- if(!length(d))list(p.value=NA_real_,statistic=NA_real_) else if(all(d==0))list(p.value=1,statistic=0) else
  wilcox.test(d,mu=0,alternative='two.sided',exact=FALSE,correct=TRUE)
 data.frame(score=rownames(readouts)[i],unit=units[i],n_pairs=length(d),Primary_median=if(length(d))median(a) else NA_real_,
  Node_median=if(length(d))median(b) else NA_real_,median_delta=if(length(d))median(d) else NA_real_,mean_delta=mean(d),
  Node_higher_pct=if(length(d))100*mean(d>0) else NA_real_,Node_lower_pct=if(length(d))100*mean(d<0) else NA_real_,
  raw_P=test$p.value,Wilcoxon_W=unname(test$statistic),zero_deltas=sum(d==0),abs_delta_ties=anyDuplicated(abs(d[d!=0]))>0,
  display_SD=sd(readouts[i,],na.rm=TRUE),median_delta_over_pooled_readout_SD=median(d)/sd(readouts[i,],na.rm=TRUE),
  direction=if(!length(d))'NOT_EVALUABLE' else if(median(d)>0)'Node_higher' else if(median(d)<0)'Node_lower' else 'Median_zero',
  test='Two-sided paired Wilcoxon, exact=FALSE, continuity correction',row.names=NULL)
}))
tests$BH_q <- p.adjust(tests$raw_P,'BH');tests$BH_family_N <- sum(is.finite(tests$raw_P));tests$BH_family <- 'All study evaluable continuous readouts'
tests$significant <- is.finite(tests$BH_q)&tests$BH_q<.05
write_csv(tests,'molecular_PAIRED_SIGNATURE_TESTS.csv')
master <- pd
for(score in rownames(readouts)){
 master[[paste0(score,'_Primary')]] <- primary[score,];master[[paste0(score,'_Node')]] <- node[score,];master[[paste0(score,'_delta')]] <- delta[score,]
}
write_csv(master,'molecular_PATIENT_PAIRED_SIGNATURES.csv')
drivers <- c('StromalScore','B_lineage','Fibroblast','Endothelial','EMT','E2F','G2M')
driver <- do.call(rbind,lapply(drivers,function(s){
 dis <- pd$evaluable & !is.na(pd$PAM50_discordant) & pd$PAM50_discordant;con <- pd$evaluable & !is.na(pd$PAM50_discordant) & !pd$PAM50_discordant
 a <- abs(delta[s,dis]);b <- abs(delta[s,con]);a <- a[is.finite(a)];b <- b[is.finite(b)]
 test <- if(length(a)&&length(b))wilcox.test(a,b,exact=FALSE,correct=TRUE) else list(p.value=NA_real_,statistic=NA_real_)
 data.frame(score=s,discordant_N=length(a),concordant_N=length(b),discordant_median_abs_delta=if(length(a))median(a) else NA_real_,
  concordant_median_abs_delta=if(length(b))median(b) else NA_real_,raw_P=test$p.value,
  status=if(length(a)&&length(b))'EVALUABLE' else 'NOT_RUN_PAM50_GROUPS_UNAVAILABLE',test='Two-sided rank-sum, exact=FALSE, continuity correction')
}))
driver$BH_q <- p.adjust(driver$raw_P,'BH');driver$BH_family <- 'Seven prespecified discordance drivers';write_csv(driver,'molecular_DISCORDANCE_DRIVERS.csv')
corrs <- data.frame(x=c('StromalScore','StromalScore','StromalScore','B_lineage','Fibroblast','Fibroblast','E2F','E2F'),
 y=c('EMT','ECM','Adipogenesis','BCR','ECM','EMT','MKI67','G2M'))
correlations <- do.call(rbind,lapply(seq_len(nrow(corrs)),function(i){
 a <- delta[corrs$x[i],];b <- delta[corrs$y[i],];keep <- is.finite(a)&is.finite(b)
 test <- if(sum(keep)>2&&sd(a[keep])>0&&sd(b[keep])>0)cor.test(a[keep],b[keep],method='spearman',exact=FALSE) else list(estimate=NA_real_,p.value=NA_real_)
 xg <- sets[[corrs$x[i]]];yg <- sets[[corrs$y[i]]];overlap <- if(!is.null(xg)&&!is.null(yg))length(intersect(xg,yg)) else NA_integer_
 data.frame(x_delta=corrs$x[i],y_delta=corrs$y[i],n_pairs=sum(keep),rho=unname(test$estimate),raw_P=test$p.value,
  shared_module_genes_N=overlap,test='Two-sided Spearman, exact=FALSE',interpretation='Association only; shared markers/data may contribute')
}))
correlations$BH_q <- p.adjust(correlations$raw_P,'BH');correlations$BH_family_N <- 8
write_csv(correlations,'molecular_DELTA_CORRELATIONS.csv')
saveRDS(list(samples=samples,pair_map=pm,pam=pam,discordance=pd,transition=transition,readouts=readouts,
 primary=primary,node=node,delta=delta,tests=tests,drivers=driver,correlations=correlations,sets=sets),file.path(out,'work/ANALYSIS_RESULTS.rds'))
save_json(list(n_pairs=36,continuous_readouts=23,BH_family_evaluable_N=sum(is.finite(tests$raw_P)),
 significant_continuous=sum(tests$significant),PAM50_classified_pairs=n,signature_modules=11,
 new_DEG_fits=0,new_ML_fits=0,survival=FALSE,MR=FALSE,new_TME_estimation=FALSE,ROR='SKIPPED',
 correlation_family_N=8,discordance_driver_tests=sum(is.finite(driver$raw_P))), 'work/ANALYSIS_COMPLETION.json')
capture.output(sessionInfo(),file=file.path(out,'work/ANALYSIS_R_SESSION.txt'))
cat('study analysis completed; continuous significant',sum(tests$significant),'/23; PAM50 pairs',n,'\n')
