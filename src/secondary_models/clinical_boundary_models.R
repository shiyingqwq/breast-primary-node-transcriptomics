options(warn=1)
root<-normalizePath(Sys.getenv('BCPN_WORKSPACE'),winslash='/');ev<-file.path(root,'clinical')
if(nzchar(Sys.getenv("BCPN_R_LIBRARY"))) .libPaths(c(Sys.getenv("BCPN_R_LIBRARY"),.libPaths()))
library(survival);library(pROC);library(jsonlite)
stopifnot(!file.exists(file.path(root,'clinical_PRIMARY_NODAL_ASSOCIATION.csv')))
started<-format(Sys.time(),tz='UTC',usetz=TRUE)
d<-read.csv(file.path(ev,'CLINICAL_JOINED.csv'),stringsAsFactors=FALSE)
names(d)[names(d)=='standardized_NodeLikeScore']<-'Z'
d$grade<-factor(d$grade,levels=c(1,2,3));d$ER<-factor(d$ER,levels=c(0,1));d$HER2<-factor(d$HER2,levels=c(0,1))
base<-c('age10','size10','grade','ER','HER2');details<-list();coeffs<-list();ph<-list();preds<-list();nodal<-list();support<-list();increment<-list();surv<-list();multiplicity<-list();models<-list();k<-0
write_table<-function(name,rows)write.csv(if(is.data.frame(rows))rows else do.call(rbind,rows),file.path(root,name),row.names=FALSE,na='')
run_fit<-function(formula,data,id,cohort,endpoint,kind){
 warnings<-character()
 model<-withCallingHandlers(if(kind=='logistic')glm(formula,data=data,family=binomial(link='logit'),control=glm.control(epsilon=1e-8,maxit=100)) else coxph(formula,data=data,ties='efron',x=TRUE,model=TRUE,control=coxph.control()),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 sm<-summary(model);cm<-sm$coefficients
 converged<-if(kind=='logistic')isTRUE(model$converged) else all(is.finite(coef(model))) && !length(warnings)
 stopifnot(all(is.finite(coef(model))),converged)
 if(kind=='logistic'){b<-cm[,1];se<-cm[,2];p<-cm[,4]}else{b<-cm[,1];se<-cm[,3];p<-cm[,5]}
 out<-data.frame(model_id=id,cohort=cohort,endpoint=endpoint,formula=paste(deparse(formula),collapse=''),analysis_N=nrow(data),event_or_positive_N=if(kind=='logistic')sum(data$node_positive) else sum(data[[paste0(endpoint,'_event')]]),
  term=rownames(cm),beta=b,SE=se,effect=exp(b),CI_low=exp(b-qnorm(.975)*se),CI_high=exp(b+qnorm(.975)*se),raw_P=p,effect_type=if(kind=='logistic')'OR' else 'HR',convergence=converged,warnings=paste(warnings,collapse=' | '),log_likelihood=as.numeric(logLik(model)),AIC=AIC(model),stringsAsFactors=FALSE)
 coeffs[[id]]<<-out;models[[id]]<<-model
 details[[id]]<<-data.frame(model_id=id,cohort=cohort,endpoint=endpoint,formula=paste(deparse(formula),collapse=''),N=nrow(data),positive_or_event_N=out$event_or_positive_N[1],iterations=paste(model$iter,collapse=';'),convergence=converged,warnings=paste(warnings,collapse=' | '))
 if(kind=='cox'){
  z<-cox.zph(model,transform='km',terms=TRUE,singledf=FALSE);tab<-z$table
  for(term in rownames(tab))ph[[paste(id,term)]]<<-data.frame(model_id=id,cohort=cohort,endpoint=endpoint,term=term,chisq=tab[term,'chisq'],df=tab[term,'df'],P=tab[term,'p'],PH_concern=tab[term,'p']<.05,diagnostic='Schoenfeld; transform=km')
 }
 model
}
score_est<-function(id){r<-coeffs[[id]];r[r$term=='Z',]}
lrt<-function(a,b){stat<-2*(as.numeric(logLik(b))-as.numeric(logLik(a)));df<-attr(logLik(b),'df')-attr(logLik(a),'df');c(stat=stat,df=df,P=pchisq(stat,df,lower.tail=FALSE))}
perf<-function(model,dd,id,cohort){
 lp<-predict(model,newdata=dd,type='link');prob<-plogis(lp);y<-dd$node_positive
 rocobj<-roc(y,prob,levels=c(0,1),direction='<',quiet=TRUE);ci<-as.numeric(ci.auc(rocobj,method='delong'))
 cal0<-glm(y~1+offset(lp),family=binomial());cal1<-glm(y~lp,family=binomial())
 preds[[id]]<<-data.frame(cohort=cohort,patient_id=dd$patient_id,model_id=id,node_positive=y,probability=prob,linear_predictor=lp)
 data.frame(cohort=cohort,model_id=id,N=nrow(dd),positive_N=sum(y),AUC=as.numeric(auc(rocobj)),AUC_CI_low=ci[1],AUC_CI_high=ci[3],Brier=mean((y-prob)^2),calibration_intercept=unname(coef(cal0)[1]),calibration_slope=unname(coef(cal1)[2]),calibration_free_intercept=unname(coef(cal1)[1]),performance_scope='APPARENT same-cohort fitted performance; no independent calibration claim',log_likelihood=as.numeric(logLik(model)),AIC=AIC(model),stringsAsFactors=FALSE)
}
for(cohort in c('GSE96058','METABRIC')){
 cc<-d[d$cohort==cohort,];nn<-cc[!is.na(cc$node_positive),];complete<-nn[complete.cases(nn[,base]),]
 id<-paste0(cohort,'_NODAL_PRIMARY');m<-run_fit(node_positive~Z+age10+size10+grade+ER+HER2,complete,id,cohort,'NODAL','logistic')
 nodal[[cohort]]<-score_est(id)
 b<-run_fit(node_positive~age10+size10+grade+ER+HER2,complete,paste0(cohort,'_NODAL_BASELINE'),cohort,'NODAL','logistic')
 s<-run_fit(node_positive~Z,complete,paste0(cohort,'_NODAL_SCORE_ONLY_CC'),cohort,'NODAL','logistic')
 u<-run_fit(node_positive~Z,nn,paste0(cohort,'_NODAL_SUPPORTIVE'),cohort,'NODAL','logistic');support[[cohort]]<-score_est(paste0(cohort,'_NODAL_SUPPORTIVE'))
 pb<-perf(b,complete,paste0(cohort,'_NODAL_BASELINE'),cohort);pm<-perf(m,complete,id,cohort);ps<-perf(s,complete,paste0(cohort,'_NODAL_SCORE_ONLY_CC'),cohort)
 lt<-lrt(b,m)
 for(obj in list(pb,pm,ps)){
  obj$delta_AUC_vs_baseline<-obj$AUC-pb$AUC;obj$delta_Brier_vs_baseline<-obj$Brier-pb$Brier;obj$incremental_LRT_stat<-if(obj$model_id==id)lt['stat'] else NA_real_;obj$incremental_LRT_df<-if(obj$model_id==id)lt['df'] else NA_real_;obj$incremental_LRT_P<-if(obj$model_id==id)lt['P'] else NA_real_
  increment[[obj$model_id]]<-obj
 }
 cat('Completed nodal models:',cohort,'N primary=',nrow(complete),'supportive=',nrow(nn),'\n');flush.console()
}
primary<-do.call(rbind,nodal);primary$Holm_adjusted_P<-p.adjust(primary$raw_P,method='holm',n=2);primary$positive_direction<-primary$beta>0;primary$primary_association_supported<-primary$positive_direction & primary$Holm_adjusted_P<.05
primary$primary_family<-'Two adjusted score Wald tests; Holm'
write_table('clinical_PRIMARY_NODAL_ASSOCIATION.csv',primary);write_table('clinical_SUPPORTIVE_UNADJUSTED_NODAL.csv',support)
np<-c(primary$raw_P,vapply(c('GSE96058','METABRIC'),function(c)increment[[paste0(c,'_NODAL_PRIMARY')]]$incremental_LRT_P,numeric(1)))
nq<-p.adjust(np,method='BH',n=4)
for(i in seq_len(2))increment[[paste0(c('GSE96058','METABRIC')[i],'_NODAL_PRIMARY')]]$incremental_LRT_BH_q<-nq[i+2]
for(id in names(increment))if(is.null(increment[[id]]$incremental_LRT_BH_q))increment[[id]]$incremental_LRT_BH_q<-NA_real_
for(i in seq_len(2)){
 multiplicity[[paste0('PRIMARY',i)]]<-data.frame(family='PRIMARY_NODAL',test=paste(primary$cohort[i],'adjusted score Wald'),raw_P=primary$raw_P[i],correction='Holm',family_N=2,adjusted_value=primary$Holm_adjusted_P[i],interpretation='Primary comparison')
}
tests<-c(paste(primary$cohort,'adjusted score Wald'),paste(c('GSE96058','METABRIC'),'incremental LRT'))
for(i in seq_len(4))multiplicity[[paste0('SECONDARY_NODAL',i)]]<-data.frame(family='SECONDARY_NODAL',test=tests[i],raw_P=np[i],correction='BH',family_N=4,adjusted_value=nq[i],interpretation='Secondary; same Wald tests also have primary Holm labels')
write_table('clinical_INCREMENTAL_VALUE.csv',increment)
nodalcoef<-do.call(rbind,coeffs);nodalcoef<-nodalcoef[nodalcoef$endpoint=='NODAL',]
nodalcoef$secondary_nodal_BH_q<-NA_real_;nodalcoef$primary_Holm_P<-NA_real_
for(i in seq_len(2)){
 ii<-nodalcoef$model_id==primary$model_id[i]&nodalcoef$term=='Z';nodalcoef$secondary_nodal_BH_q[ii]<-nq[i];nodalcoef$primary_Holm_P[ii]<-primary$Holm_adjusted_P[i]
}
write_table('clinical_ADJUSTED_NODAL_MODELS.csv',nodalcoef)
write_table('clinical_EXTERNAL_MODEL_TRANSPORT.csv',data.frame(status='NOT_EVALUABLE',reason='Insufficient cross-cohort clinical coding'))
train<-list(cohort='GSE96058',models=lapply(models[c('GSE96058_NODAL_BASELINE','GSE96058_NODAL_PRIMARY')],function(m)list(formula=paste(deparse(formula(m)),collapse=''),coefficients=as.list(coef(m)),contrasts=m$contrasts,xlevels=m$xlevels)),constants=read.csv(file.path(root,'clinical_SCORE_STANDARDIZATION_CONSTANTS.csv')))
write_json(train,file.path(root,'clinical_TRAINING_MODEL_PARAMETER_SETTINGS.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
for(spec in list(c('METABRIC','RFS'),c('GSE96058','OS'),c('METABRIC','OS'))){
 cohort<-spec[1];endpoint<-spec[2];cc<-d[d$cohort==cohort,];time<-paste0(endpoint,'_time');event<-paste0(endpoint,'_event')
 ss<-cc[!is.na(cc[[time]]) & cc[[time]]>0 & !is.na(cc[[event]]) & cc[[event]]%in%c(0,1),]
 complete<-ss[complete.cases(ss[,c(base,'node_positive')]),]
 formu<-as.formula(paste0('Surv(',time,',',event,')~Z'))
 forma<-as.formula(paste0('Surv(',time,',',event,')~Z+age10+size10+grade+ER+HER2+node_positive'))
 formb<-as.formula(paste0('Surv(',time,',',event,')~age10+size10+grade+ER+HER2+node_positive'))
 pre<-paste(cohort,endpoint,sep='_')
 u<-run_fit(formu,ss,paste0(pre,'_UNADJUSTED'),cohort,endpoint,'cox');a<-run_fit(forma,complete,paste0(pre,'_ADJUSTED'),cohort,endpoint,'cox');b<-run_fit(formb,complete,paste0(pre,'_BASELINE'),cohort,endpoint,'cox')
 lt<-lrt(b,a)
 for(ad in c(FALSE,TRUE)){
  id<-paste0(pre,if(ad)'_ADJUSTED' else '_UNADJUSTED');out<-score_est(id);out$adjusted<-ad;out$survival_BH_q<-NA_real_;out$incremental_LRT_stat<-if(ad)lt['stat'] else NA_real_;out$incremental_LRT_P<-if(ad)lt['P'] else NA_real_
  out$score_PH_P<-ph[[paste(id,'Z')]]$P;out$global_PH_P<-ph[[paste(id,'GLOBAL')]]$P;out$score_PH_concern<-out$score_PH_P<.05;out$global_PH_concern<-out$global_PH_P<.05
  surv[[id]]<-out
 }
 cat('Completed survival:',pre,'N=',nrow(ss),'events=',sum(ss[[event]]),'\n');flush.console()
}
unadj<-which(!vapply(surv,function(r)r$adjusted,logical(1)));sp<-vapply(surv[unadj],function(r)r$raw_P,numeric(1));sq<-p.adjust(sp,method='BH',n=3)
for(i in seq_along(unadj)){
 j<-unadj[i];surv[[j]]$survival_BH_q<-sq[i]
 multiplicity[[paste0('SURVIVAL',i)]]<-data.frame(family='SECONDARY_SURVIVAL',test=paste(surv[[j]]$cohort,surv[[j]]$endpoint,'unadjusted score Wald'),raw_P=sp[i],correction='BH',family_N=3,adjusted_value=sq[i],interpretation='Secondary survival; adjusted tests supporting and raw')
}
write_table('clinical_SURVIVAL_RESULTS.csv',surv);write_table('clinical_PH_ASSUMPTION_DETAILS.csv',ph);write_table('clinical_MULTIPLICITY_LEDGER.csv',multiplicity);write_table('clinical_ALL_MODEL_COEFFICIENTS.csv',coeffs);write_table('clinical_MODEL_DETAILS.csv',details)
write.csv(do.call(rbind,preds),file.path(ev,'NODAL_MODEL_PREDICTIONS.csv'),row.names=FALSE)
write_json(list(started_utc=started,completed_utc=format(Sys.time(),tz='UTC',usetz=TRUE),primary=primary,logistic_models=8,Cox_models=9,calibration_diagnostic_models=12,alternative_scores=0,subgroups=0,ML=0,drug=0,KM=0,model_warnings=lapply(details,function(x)x$warnings)),file.path(ev,'ANALYSIS_COMPLETION_SUMMARY.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
capture.output(sessionInfo(),file=file.path(ev,'ANALYSIS_SESSION.txt'))
cat('Specified clinical models completed.\n')
