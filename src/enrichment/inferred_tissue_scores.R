source(file.path(Sys.getenv('BCPN_REPOSITORY'), 'src/enrichment/runtime.R'))
library(tidyestimate)
input <- readRDS(file.path(out,'work/INPUTS.rds'))
x <- input$x;m <- input$mapping;samples <- input$samples
use <- m$annotation_usable
xt <- x[use,,drop=FALSE];rownames(xt)<-m$approved_symbol[use]
stopifnot(!anyDuplicated(rownames(xt)),ncol(xt)==72)
df <- data.frame(hgnc_symbol=rownames(xt),xt,check.names=FALSE)
est_filtered <- tidyestimate::filter_common_genes(df,id='hgnc_symbol',tidy=TRUE,tell_missing=TRUE,find_alias=FALSE)
stopifnot(ncol(est_filtered)==73,all(vapply(est_filtered[,-1],is.numeric,TRUE)))
est <- as.data.frame(tidyestimate::estimate_score(est_filtered,is_affymetrix=TRUE))
est <- est[match(samples$GSM,est$sample),]
stopifnot(identical(samples$GSM,est$sample))
est_scores <- t(as.matrix(est[,c('immune','stromal','estimate','purity')]))
rownames(est_scores)<-c('ImmuneScore','StromalScore','ESTIMATEScore','TumorPurity')
colnames(est_scores)<-samples$GSM
source(file.path(src,'MCPcounter.R'))
marker_genes <- read.delim(file.path(src,'MCP_genes.tsv'),check.names=FALSE,colClasses='character')
marker_probes <- read.delim(file.path(src,'MCP_probesets.tsv'),header=FALSE,colClasses='character')
mcp <- as.matrix(MCPcounter.estimate(xt,featuresType='HUGO_symbols',
                                   genes=marker_genes,probesets=marker_probes))
stopifnot(ncol(mcp)==72,all(is.finite(mcp)))
mcp <- mcp[,samples$GSM,drop=FALSE]
expected_populations <- unique(marker_genes[['Cell population']])
missing_populations <- setdiff(expected_populations,rownames(mcp))
if(length(missing_populations)) {
  extra <- matrix(NA_real_,nrow=length(missing_populations),ncol=72,
                  dimnames=list(missing_populations,samples$GSM))
  mcp <- rbind(mcp,extra)
  cat('MCP populations with no observed markers:',paste(missing_populations,collapse='; '),'\n')
}
patients <- unique(samples$patient)
pair_map <- data.frame(patient=patients,
  Primary=sapply(patients,function(p)samples$GSM[samples$patient==p&samples$tissue=='Primary']),
  Node=sapply(patients,function(p)samples$GSM[samples$patient==p&samples$tissue=='Node']))
analyze <- function(scores,method) {
  rows <- lapply(rownames(scores),function(score) {
    primary <- scores[score,pair_map$Primary];node <- scores[score,pair_map$Node]
    keep <- is.finite(primary)&is.finite(node);delta <- node[keep]-primary[keep]
    test <- if(length(delta))wilcox.test(delta,mu=0,alternative='two.sided',exact=FALSE,correct=TRUE) else
      list(p.value=NA_real_,statistic=NA_real_)
    data.frame(method=method,score=score,n_pairs=sum(keep),median_Primary=median(primary[keep]),
      median_Node=median(node[keep]),mean_delta=mean(delta),median_delta=median(delta),
      direction=if(!length(delta))'NOT_EVALUABLE' else if(median(delta)>0)'Node_higher' else if(median(delta)<0)'Node_lower' else 'Median_zero',
      raw_P=test$p.value,Wilcoxon_W=unname(test$statistic),zero_deltas=sum(delta==0),
      absolute_delta_ties=anyDuplicated(abs(delta[delta!=0]))>0,
      test='Two-sided paired Wilcoxon; asymptotic; continuity correction; zeros omitted',row.names=NULL)
  })
  result <- do.call(rbind,rows);result$BH_q <- p.adjust(result$raw_P,'BH',n=nrow(scores))
  result$significant <- is.finite(result$BH_q)&result$BH_q<.05;result$BH_family <- method
  write_csv(result,paste0('enrichment_',method,'_PAIRED_RESULTS.csv'))
  list(result=result,scores=scores,delta=scores[,pair_map$Node,drop=FALSE]-scores[,pair_map$Primary,drop=FALSE])
}
es <- analyze(est_scores,'ESTIMATE');mc <- analyze(mcp,'MCP_COUNTER')
colnames(es$delta)<-colnames(mc$delta)<-patients
marker_coverage <- do.call(rbind,lapply(unique(marker_genes[['Cell population']]),function(pop) {
  markers <- unique(marker_genes[['HUGO symbols']][marker_genes[['Cell population']]==pop])
  data.frame(method='MCP_COUNTER',score=pop,signature_genes_N=length(markers),
    observed_genes_N=sum(markers%in%rownames(xt)),missing_genes=paste(setdiff(markers,rownames(xt)),collapse=';'))
}))
write_csv(marker_coverage,'work/MCP_MARKER_COVERAGE.csv')
estimate_coverage <- do.call(rbind,lapply(names(tidyestimate::gene_sets),function(n){
  markers <- unique(na.omit(tidyestimate::gene_sets[[n]]))
  data.frame(method='ESTIMATE',score=n,signature_genes_N=length(markers),
    observed_genes_N=sum(markers%in%est_filtered[[1]]),missing_genes=paste(setdiff(markers,est_filtered[[1]]),collapse=';'))
}))
write_csv(estimate_coverage,'work/ESTIMATE_MARKER_COVERAGE.csv')
saveRDS(list(ESTIMATE=es,MCP_COUNTER=mc,pair_map=pair_map,samples=samples),file.path(out,'work/TME_RESULTS.rds'))
write_csv(data.frame(score=rownames(es$delta),es$delta,check.names=FALSE),'work/ESTIMATE_PAIRED_DELTAS.csv')
write_csv(data.frame(score=rownames(mc$delta),mc$delta,check.names=FALSE),'work/MCP_PAIRED_DELTAS.csv')
write_csv(data.frame(score=rownames(est_scores),est_scores,check.names=FALSE),'work/ESTIMATE_SAMPLE_SCORES.csv')
write_csv(data.frame(score=rownames(mcp),mcp,check.names=FALSE),'work/MCP_SAMPLE_SCORES.csv')
save_json(list(estimate_implementation='tidyestimate CRAN 1.1.1; original R-Forge download unavailable',
  input_gene_N=nrow(xt),ESTIMATE_common_background_N=nrow(est_filtered),ESTIMATE_family_N=4,
  MCP_family_N=nrow(mcp),MCP_evaluable_population_N=sum(mc$result$n_pairs>0),
  MCP_not_evaluable=missing_populations,xCell='SKIPPED: optional package and GSVA absent',
  clinical_analysis=FALSE,paired_patients=36), 'work/TME_COMPLETE.json')
capture.output(sessionInfo(),file=file.path(out,'work/TME_SESSION.txt'))
cat('TME COMPLETE: ESTIMATE significant',sum(es$result$significant),'MCP significant',sum(mc$result$significant),'\n')
