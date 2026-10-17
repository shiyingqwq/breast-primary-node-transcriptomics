source(file.path(Sys.getenv('BCPN_REPOSITORY'), 'src/clinical_context/runtime.R'))
p <- readRDS(file.path(out,'work/PREPARED.rds'))
input <- p$input;pm <- p$pm;primary <- p$primary;node <- p$node
samples <- input$samples;pam <- p$pam[match(samples$GSM,p$pam$GSM),]
stopifnot(identical(pam$GSM,samples$GSM),!anyNA(pam$PAM50_subtype))
if(!file.exists(file.path(out,'context_RISK_SIGNATURE_DISCORDANCE.csv'))){
risk <- data.frame(row.names=samples$GSM)
risk[['ROR-S-like']] <- .05*pam$Basal+.12*pam$Her2-.34*pam$LumA+.23*pam$LumB
for(key in c('GGI','RS')){
 d <- p$definitions[[key]]
 z <- sweep(sweep(input$x[d$index,,drop=FALSE],1,d$mean,'-'),1,d$sd,'/')
 stopifnot(all(is.finite(z)))
 if(key=='GGI')risk[['GGI-like']] <- colMeans(z*d$coefficients)
 if(key=='RS'){
  risk[['RS-like']] <- colSums(z*d$coefficients)
  groups <- split(seq_along(d$groups),d$groups)
  modules <- data.frame(samples,check.names=FALSE)
  for(g in names(groups)){
   ii <- groups[[g]]
   modules[[paste0(g,'_weighted_contribution')]] <- colSums(z[ii,,drop=FALSE]*d$coefficients[ii])
  }
  write_csv(modules,'work/RS_MODULE_COMPONENTS.csv')
 }
}
stopifnot(ncol(risk)==3,nrow(risk)==72,all(is.finite(as.matrix(risk))))
threshold <- do.call(rbind,lapply(names(risk),function(n){q<-quantile(risk[primary,n],c(1/3,2/3),type=7);
 data.frame(signature=n,lower_tertile=unname(q[1]),upper_tertile=unname(q[2]),reference_N=36,
 reference_tissue='Primary',type=7,low_rule='score<=lower',intermediate_rule='lower<score<=upper',high_rule='score>upper')}))
stopifnot(all(threshold$upper_tertile>threshold$lower_tertile))
write_csv(threshold,'work/PRIMARY_RISK_CATEGORY_CUTOFFS.csv')
categorize <- function(x,q)factor(ifelse(x<=q[1],'Low',ifelse(x<=q[2],'Intermediate','High')),
 levels=c('Low','Intermediate','High'),ordered=TRUE)
kappa <- function(a,b,weighted=TRUE){
 lev<-if(weighted)c('Low','Intermediate','High') else c('Luminal A','Luminal B','HER2-enriched','Basal-like','Normal-like')
 tab<-table(factor(a,levels=lev),factor(b,levels=lev));n<-sum(tab)
 w<-if(weighted)1-abs(outer(seq_along(lev),seq_along(lev),'-'))/(length(lev)-1) else diag(length(lev))
 po<-sum(tab*w)/n;pe<-sum(outer(rowSums(tab),colSums(tab))*w)/n^2
 c(kappa=(po-pe)/(1-pe),agreement=po,chance_agreement=pe)
}
wilcox <- function(a,b){delta<-b-a;nonzero<-delta[delta!=0];
 if(!length(nonzero))return(list(P=1,statistic=0,zeros=length(delta),ties=FALSE))
 v<-wilcox.test(b,a,paired=TRUE,alternative='two.sided',exact=FALSE,correct=TRUE)
 list(P=v$p.value,statistic=unname(v$statistic),zeros=sum(delta==0),ties=anyDuplicated(abs(nonzero))>0)}
long <- paired <- result <- transitions <- list()
wide <- data.frame(patient=pm$patient,Primary_GSM=pm$Primary,Node_GSM=pm$Node,
 Primary_PAM50=pam$PAM50_subtype[primary],Node_PAM50=pam$PAM50_subtype[node])
wide$PAM50_discordant <- wide$Primary_PAM50!=wide$Node_PAM50
for(n in names(risk)){
 q<-unlist(threshold[threshold$signature==n,c('lower_tertile','upper_tertile')],use.names=FALSE)
 catg<-categorize(risk[[n]],q);a<-risk[primary,n];b<-risk[node,n];delta<-b-a
 ca<-catg[primary];cb<-catg[node];order_delta<-as.integer(cb)-as.integer(ca)
 kk<-kappa(ca,cb);wt<-wilcox(a,b)
 long[[n]]<-data.frame(samples,signature=n,score=risk[[n]],risk_category=as.character(catg),
  lower_cutoff=q[1],upper_cutoff=q[2],implementation='Primary-tertile research approximation; not a validated clinical assay')
 paired[[n]]<-data.frame(patient=pm$patient,signature=n,Primary_GSM=pm$Primary,Node_GSM=pm$Node,
  Primary_score=a,Node_score=b,delta=delta,Primary_category=as.character(ca),Node_category=as.character(cb),
  category_discordant=ca!=cb,category_step=order_delta,transition=ifelse(order_delta>0,'Upward',ifelse(order_delta<0,'Downward','Stable')))
 result[[n]]<-data.frame(signature=n,status='EVALUABLE_RESEARCH_APPROXIMATION',evaluable_pairs=36,
  concordant_N=sum(ca==cb),discordant_N=sum(ca!=cb),discordance_pct=100*mean(ca!=cb),
  upward_shift_N=sum(order_delta>0),downward_shift_N=sum(order_delta<0),upward_shift_pct=100*mean(order_delta>0),
  downward_shift_pct=100*mean(order_delta<0),kappa=unname(kk[1]),kappa_type='Linear weighted Cohen kappa',
  Primary_median=median(a),Node_median=median(b),median_delta=median(delta),raw_P=wt$P,BH_q=NA_real_,
  Node_higher_pct=100*mean(delta>0),Node_lower_pct=100*mean(delta<0),Wilcoxon_V=wt$statistic,
  zero_deltas=wt$zeros,tied_absolute_nonzero_deltas=wt$ties,exact=FALSE,continuity_correction=TRUE,
  low_to_intermediate_or_high_N=sum(ca=='Low'&cb!='Low'),intermediate_to_high_N=sum(ca=='Intermediate'&cb=='High'),
  high_to_low_or_intermediate_N=sum(ca=='High'&cb!='High'),low_to_high_N=sum(ca=='Low'&cb=='High'),
  Primary_low_N=sum(ca=='Low'),Primary_intermediate_N=sum(ca=='Intermediate'),Primary_high_N=sum(ca=='High'))
 transitions[[n]]<-as.data.frame(table(Primary=factor(ca,levels=c('Low','Intermediate','High')),
  Node=factor(cb,levels=c('Low','Intermediate','High'))));transitions[[n]]$signature<-n
 wide[[paste0(n,'_Primary_score')]]<-a;wide[[paste0(n,'_Node_score')]]<-b;wide[[paste0(n,'_delta')]]<-delta
 wide[[paste0(n,'_Primary_category')]]<-as.character(ca);wide[[paste0(n,'_Node_category')]]<-as.character(cb)
 wide[[paste0(n,'_transition')]]<-paired[[n]]$transition
}
results<-do.call(rbind,result);results$BH_q<-p.adjust(results$raw_P,'BH');rownames(results)<-NULL
blank<-results[1,];blank[1,]<-NA;blank$signature<-'PAM50';blank$status<-'REUSED_molecular_43_OF_50_RESEARCH_CLASSIFICATION'
blank$evaluable_pairs<-36;blank$concordant_N<-28;blank$discordant_N<-8;blank$discordance_pct<-100*8/36
blank$kappa<-unname(kappa(wide$Primary_PAM50,wide$Node_PAM50,FALSE)[1]);blank$kappa_type<-'Ordinary Cohen kappa'
skip<-blank;skip[1,]<-NA;skip$signature<-'70-gene research score';skip$status<-'SKIP: 47/70 original source features; below 80%';skip$evaluable_pairs<-0
results<-rbind(blank,results[1:2,],skip,results[3,]);rownames(results)<-NULL
write_csv(results,'context_RISK_SIGNATURE_DISCORDANCE.csv')
pamlong<-long[[1]];pamlong$signature<-'PAM50';pamlong$score<-NA_real_;pamlong$risk_category<-pam$PAM50_subtype
pamlong$lower_cutoff<-NA_real_;pamlong$upper_cutoff<-NA_real_;pamlong$implementation<-'Reused study 43/50-gene research subtype; not an ordered risk category'
write_csv(rbind(pamlong,do.call(rbind,long)),'context_RISK_SIGNATURE_SCORES.csv')
write_csv(wide,'context_PATIENT_RISK_TRANSITIONS.csv')
pairedlong<-do.call(rbind,paired);rownames(pairedlong)<-NULL
write_csv(pairedlong,'work/RISK_PATIENT_TRANSITIONS_LONG.csv')
write_csv(do.call(rbind,transitions),'work/RISK_CATEGORY_TRANSITION_COUNTS.csv')
cat('RISK RESULTS\n');print(results[,c('signature','evaluable_pairs','discordance_pct','upward_shift_pct','kappa','median_delta','BH_q')])
}else{
 results<-read_csv(file.path(out,'context_RISK_SIGNATURE_DISCORDANCE.csv'))
 threshold<-read_csv(file.path(out,'work/PRIMARY_RISK_CATEGORY_CUTOFFS.csv'))
 wide<-read_csv(file.path(out,'context_PATIENT_RISK_TRANSITIONS.csv'))
 pairedlong<-read_csv(file.path(out,'work/RISK_PATIENT_TRANSITIONS_LONG.csv'))
 saved<-read_csv(file.path(out,'context_RISK_SIGNATURE_SCORES.csv'))
 risk<-data.frame(row.names=samples$GSM)
 for(n in c('ROR-S-like','GGI-like','RS-like')){
  v<-saved[saved$signature==n,];risk[[n]]<-v$score[match(samples$GSM,v$GSM)]
 }
 stopifnot(all(is.finite(as.matrix(risk))),nrow(wide)==36,nrow(pairedlong)==108)
 cat('RESUME: saved risk scores, thresholds, tests and transitions reused without recomputation\n')
}

clin<-read_csv(file.path(out,'work/METABRIC_CLINICAL_MAP.csv'))
if(!is.logical(clin$nodal_count_eligible))clin$nodal_count_eligible<-tolower(clin$nodal_count_eligible)=='true'
stopifnot(!anyNA(clin$nodal_count_eligible),sum(clin$nodal_count_eligible)==1890)
genes<-read_csv(file.path(out,'work/METABRIC_SIGNATURE_GENE_Z.csv'))
zz<-as.matrix(genes[,-1]);rownames(zz)<-genes$gene
stopifnot(identical(colnames(zz),clin$SAMPLE_ID),all(is.finite(zz)),nrow(clin)==1964)
nodal<-clin
for(n in names(p$nodal_sets)){
 use<-intersect(p$nodal_sets[[n]],rownames(zz));stopifnot(length(use)>0)
 nodal[[n]]<-colMeans(zz[use,,drop=FALSE])
}
write_csv(nodal,'context_METABRIC_NODAL_BURDEN_SCORES.csv')
evaluable<-nodal[nodal$nodal_count_eligible,];stopifnot(nrow(evaluable)==1890)
groups<-factor(evaluable$nodal_group,levels=c('N0','LOW','HIGH'),ordered=TRUE)
expected<-c(ECM=-1,BCR=1,E2F=1,G2M=1,Adipogenesis=-1)
stats<-descriptives<-ordinal<-list()
for(n in names(p$nodal_sets)){
 vals<-evaluable[[n]];rho<-cor.test(vals,evaluable$positive_node_count,method='spearman',exact=FALSE)
 kw<-kruskal.test(vals,groups)
 desc<-do.call(rbind,lapply(levels(groups),function(g){v<-vals[groups==g];q<-quantile(v,c(.25,.5,.75),type=7);
  data.frame(score=n,group=g,N=length(v),Q1=unname(q[1]),median=unname(q[2]),Q3=unname(q[3]),IQR=unname(q[3]-q[1]))}))
 descriptives[[n]]<-desc
 median_delta<-diff(desc$median);monotonic<-all(expected[n]*median_delta>0)
 stats[[n]]<-data.frame(score=n,evaluable_N=length(vals),N0_N=desc$N[1],LOW_N=desc$N[2],HIGH_N=desc$N[3],
  N0_median=desc$median[1],LOW_median=desc$median[2],HIGH_median=desc$median[3],
  N0_Q1=desc$Q1[1],N0_Q3=desc$Q3[1],LOW_Q1=desc$Q1[2],LOW_Q3=desc$Q3[2],HIGH_Q1=desc$Q1[3],HIGH_Q3=desc$Q3[3],
  expected_direction=ifelse(expected[n]>0,'Increasing','Decreasing'),rho=unname(rho$estimate),
  trend_direction=ifelse(rho$estimate>0,'Increasing','Decreasing'),trend_raw_P=rho$p.value,trend_BH_q=NA_real_,
  direction_matches_expected=sign(rho$estimate)==expected[n],medians_monotonic_in_expected_direction=monotonic,
  KW_chisq=unname(kw$statistic),KW_df=unname(kw$parameter),KW_raw_P=kw$p.value,KW_BH_q=NA_real_)
 if(requireNamespace('MASS',quietly=TRUE)){
  dat<-data.frame(nodal_group=groups,standardized_score=as.numeric(scale(vals)))
  fit<-tryCatch(MASS::polr(nodal_group~standardized_score,data=dat,Hess=TRUE,method='logistic',control=list(maxit=1000)),error=function(e)e)
  if(!inherits(fit,'error')&&fit$convergence==0){
   cc<-coef(summary(fit));b<-cc['standardized_score','Value'];se<-cc['standardized_score','Std. Error']
   ordinal[[n]]<-data.frame(score=n,N=nrow(dat),OR_per_SD=exp(b),CI95_low=exp(b-1.96*se),CI95_high=exp(b+1.96*se),
    raw_P=2*pnorm(-abs(b/se)),BH_q=NA_real_,status='CONVERGED_UNIVARIABLE_EXPLORATORY',score_mean=mean(vals),score_SD=sd(vals),
    proportional_odds_assumption='Not independently checked; supplementary association only')
  }else ordinal[[n]]<-data.frame(score=n,N=nrow(dat),OR_per_SD=NA_real_,CI95_low=NA_real_,CI95_high=NA_real_,
    raw_P=NA_real_,BH_q=NA_real_,status='SKIP_FIT_FAILURE',score_mean=mean(vals),score_SD=sd(vals),proportional_odds_assumption='Not checked')
 }
}
ntests<-do.call(rbind,stats);rownames(ntests)<-NULL
ntests$trend_BH_q<-p.adjust(ntests$trend_raw_P,'BH');ntests$KW_BH_q<-p.adjust(ntests$KW_raw_P,'BH')
ntests$consistent_significant_monotonic_trend<-ntests$trend_BH_q<.05 & ntests$direction_matches_expected & ntests$medians_monotonic_in_expected_direction
write_csv(ntests,'context_NODAL_BURDEN_TESTS.csv')
write_csv(do.call(rbind,descriptives),'work/NODAL_GROUP_DESCRIPTIVES.csv')
if(length(ordinal)){
 ord<-do.call(rbind,ordinal);ord$BH_q<-p.adjust(ord$raw_P,'BH');write_csv(ord,'work/NODAL_ORDINAL_MODELS.csv')
}
cat('NODAL RESULTS\n');print(ntests[,c('score','N0_median','LOW_median','HIGH_median','rho','trend_BH_q','consistent_significant_monotonic_trend')])
saveRDS(list(risk=risk,threshold=threshold,results=results,paired=pairedlong,wide=wide,nodal=nodal,ntests=ntests,descriptives=do.call(rbind,descriptives),ordinal=if(length(ordinal))ord else NULL),file.path(out,'work/RESULTS.rds'))
