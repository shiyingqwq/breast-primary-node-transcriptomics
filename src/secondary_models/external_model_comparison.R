source(file.path(Sys.getenv('BCPN_REPOSITORY'), 'src/secondary_models/runtime.R'))
selected <- readRDS(file.path(out,'work/DEVELOPMENT_SELECTED_MODELS.rds'))
models <- selected$models
e <- read.csv(file.path(model_inputs,'METABRIC_EXPRESSION_Z.csv'),check.names=FALSE)
clinical <- read.csv(file.path(root,'clinical/CLINICAL_JOINED.csv'),stringsAsFactors=FALSE)
d <- clinical[clinical$cohort=='METABRIC',
              c('patient_id','sample_id','node_positive','age10','size10','grade','ER')]
d <- d[match(e$patient_id,d$patient_id),]
stopifnot(identical(e$patient_id,d$patient_id),!anyDuplicated(d$patient_id),nrow(e)==1964)
ok <- !is.na(d$node_positive);x <- as.matrix(e[ok,selected$feature_universe]);d <- d[ok,]
stopifnot(nrow(x)==1890,all(is.finite(x)))
lps <- lapply(models,function(m)lp_for(m,x));probs <- lapply(lps,plogis)
for(n in c('Sparse_EN','Ridge_reference')) {
  m <- models[[n]]
  native_lp <- as.numeric(predict(m$fit,newx=x,s=m$lambda,type='link'))
  stopifnot(max(abs(native_lp-lps[[n]]))<1e-8)
}
rows <- lapply(names(models),function(n) {
  m <- models[[n]];z <- metric(n,d$node_positive,probs[[n]])
  z$genes <- m$gene_N;z$development_CV_AUC <- m$development_CV_AUC
  z$alpha <- m$alpha;z$lambda <- m$lambda;z
})
external <- do.call(rbind,rows)
refauc <- external$AUC[external$model=='Ridge_reference']
old <- read.csv(file.path(root,'model_MODEL_PERFORMANCE.csv'))
stopifnot(abs(refauc-old$AUC[old$model=='MOLECULAR_EXTERNAL_ALL'])<1e-10)
external$delta_AUC_vs_ridge <- external$AUC-refauc
retention <- function(n,a) {
  if(a<.60||n>200) return('Failed compression')
  if(n<=50&&refauc-a<=.02) return('Excellent compression')
  if(n<=100&&refauc-a<=.02) return('Good compression')
  if(n<=100&&a>=.60) return('Acceptable compression')
  'Limited compression'
}
external$compression_retention <- mapply(retention,external$genes,external$AUC)
external$compression_retention[external$model=='Ridge_reference'] <- 'Reference'
candidates <- external[external$model!='Ridge_reference',]
candidates <- candidates[order(-candidates$AUC,candidates$genes,-candidates$development_CV_AUC),]
winner_name <- candidates$model[1];winner <- models[[winner_name]]
external$selected_winner <- external$model==winner_name
external$inference_scope <- 'Fixed development coefficients; METABRIC used to select between A and B; exploratory'
write.csv(external,file.path(out,'model_EXTERNAL_VALIDATION.csv'),row.names=FALSE,na='')
cat('External candidates evaluated; winner',winner_name,'genes',winner$gene_N,'AUC',candidates$AUC[1],'\n')

dev <- load_development();gd <- clinical_coding(dev$d)
gd$compressed_score <- lp_for(winner,dev$x)
gcc <- complete.cases(gd[,c('age10','size10','grade','ER')])
combined <- glm(node_positive~age10+size10+grade+ER+compressed_score,
                data=gd[gcc,],family=binomial())
stopifnot(combined$converged,all(is.finite(coef(combined))))
d <- clinical_coding(d);d$compressed_score <- lps[[winner_name]]
cc <- complete.cases(d[,c('age10','size10','grade','ER')]);dc <- d[cc,]
cp <- predict(selected$clinical_reference,newdata=dc,type='response')
bp <- predict(combined,newdata=dc,type='response')
mp <- probs[[winner_name]][cc];rp <- probs$Ridge_reference[cc]
increment <- rbind(metric('Clinical',dc$node_positive,cp),
                   metric('Compressed_molecular',dc$node_positive,mp),
                   metric('Clinical_plus_compressed',dc$node_positive,bp))
stopifnot(sum(cc)==1780,sum(gcc)==2871,
  abs(increment$AUC[1]-old$AUC[old$model=='CLINICAL_EXTERNAL_CC'])<1e-10)
rc <- roc(dc$node_positive,cp,levels=c(0,1),direction='<',quiet=TRUE)
rb <- roc(dc$node_positive,bp,levels=c(0,1),direction='<',quiet=TRUE)
dt <- roc.test(rc,rb,paired=TRUE,method='delong',alternative='two.sided')
delta <- as.numeric(auc(rb)-auc(rc))
increment$delta_AUC_vs_clinical <- increment$AUC-increment$AUC[1]
increment$DeLong_P_vs_clinical <- NA_real_
increment$DeLong_P_vs_clinical[3] <- dt$p.value
increment$winner <- winner_name
increment$inference_scope <- 'Same external complete cases; winner selected on METABRIC; DeLong P is exploratory, not selection-adjusted'
write.csv(increment,file.path(out,'model_CLINICAL_INCREMENT.csv'),row.names=FALSE,na='')
write.csv(data.frame(term=names(coef(combined)),coefficient=coef(combined),
  fit_cohort='GSE96058',development_N=sum(gcc),score='Selected compressed molecular LP'),
  file.path(out,'work/COMBINED_CLINICAL_COEFFICIENTS.csv'),row.names=FALSE)
prediction <- data.frame(patient_id=d$patient_id,sample_id=d$sample_id,node_positive=d$node_positive,
  Sparse_EN_LP=lps$Sparse_EN,Sparse_EN_probability=probs$Sparse_EN,
  TopK_LP=lps$TopK,TopK_probability=probs$TopK,
  Ridge_reference_LP=lps$Ridge_reference,Ridge_reference_probability=probs$Ridge_reference,
  clinical_complete_case=cc,clinical_probability=NA_real_,combined_probability=NA_real_)
prediction$clinical_probability[cc] <- cp;prediction$combined_probability[cc] <- bp
write.csv(prediction,file.path(out,'work/METABRIC_FIXED_PREDICTIONS.csv'),row.names=FALSE,na='')

mapping <- read.csv(file.path(root,'cohort_NODE_LIKE_SCORE_GENES.csv'),stringsAsFactors=FALSE)
ann <- read.csv(file.path(root,'bulk_INTEGRATED_CLASSES.csv'),stringsAsFactors=FALSE)
genes <- names(winner$coefficients)
ix <- match(genes,mapping$clinical_symbol);ix2 <- match(mapping$selected_symbol[ix],ann$gene_symbol)
signature <- data.frame(gene=c('(Intercept)',genes),coefficient=c(winner$intercept,winner$coefficients),
  type=c('intercept',rep('gene',length(genes))),
  sign=c(NA_character_,ifelse(winner$coefficients>0,'Positive','Negative')),
  absolute_coefficient_rank=c(NA_integer_,rank(-abs(winner$coefficients),ties.method='first')),
  discovery_selected_symbol=c(NA_character_,mapping$selected_symbol[ix]),
  discovery_bulk_direction=c(NA_character_,ann$bulk_direction[ix2]),
  discovery_integrated_class_A_D=c(NA_character_,ann$integrated_class_A_D[ix2]),
  model=winner_name,alpha=winner$alpha,lambda=winner$lambda,
  coefficient_scale='Cohort-local gene z-scores')
signature <- signature[order(is.na(signature$absolute_coefficient_rank)==FALSE,
                              signature$absolute_coefficient_rank),]
write.csv(signature,file.path(out,'model_FINAL_COMPRESSED_SIGNATURE.csv'),row.names=FALSE,na='')
final <- list(winner=winner,winner_name=winner_name,all_development_selected_models=selected,
  clinical_model=selected$clinical_reference,combined_model=combined,
  winner_selection='Highest METABRIC AUC, then fewer genes, then development CV AUC',
  METABRIC_used_for_winner_selection=TRUE,METABRIC_coefficient_refit=FALSE,
  inference_scope='Exploratory selection-conditioned METABRIC performance; not a new confirmatory independent test',
  development_scalers=fromJSON(file.path(model_inputs,'GSE96058_SCALERS.json')),
  external_scalers=fromJSON(file.path(model_inputs,'METABRIC_SCALERS.json')),
  clinical_coding='Age years/10; size mm/10; grade factor 1/2/3 ref1; ER factor 0/1 ref0',
  external_results=external,clinical_increment=increment,
  completed_UTC=format(Sys.time(),tz='UTC',usetz=TRUE))
saveRDS(final,file.path(out,'model_FINAL_MODEL.rds'))
a <- candidates$AUC[1];n <- winner$gene_N
write_json(list(winner=winner_name,winner_method=winner$method,
  winner_genes=n,gene_names=genes,external_N=nrow(d),external_clinical_N=sum(cc),
  development_clinical_N=sum(gcc),external_AUC=a,external_AUC_CI_low=candidates$AUC_CI_low[1],
  external_AUC_CI_high=candidates$AUC_CI_high[1],clinical_AUC=increment$AUC[1],
  combined_AUC=increment$AUC[3],delta_AUC=delta,DeLong_P=dt$p.value,
  new_external_molecular_models_evaluated=2,combined_models_fit=1,external_refits=0,
  external_winner_selection=TRUE,survival_analyses=0,subgroups=0,pathway_analyses=0,
  downloads=0),file.path(out,'work/VALIDATION_COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
capture.output(sessionInfo(),file=file.path(out,'work/VALIDATION_SESSION.txt'))
cat('Model comparisons completed; deltaAUC',delta,'DeLong P',dt$p.value,'\n')
