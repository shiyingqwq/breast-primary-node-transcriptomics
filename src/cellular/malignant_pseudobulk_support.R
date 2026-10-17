options(stringsAsFactors=FALSE)
root <- normalizePath(Sys.getenv("BCPN_WORKSPACE"),winslash="/",mustWork=TRUE)
if(nzchar(Sys.getenv("BCPN_R_LIBRARY"))) .libPaths(c(Sys.getenv("BCPN_R_LIBRARY"),.libPaths()))
library(edgeR)
ev <- file.path(root,"cellular/malignant")
read_csv <- function(path) {
 d <- read.csv(path,check.names=FALSE,encoding="UTF-8")
 names(d) <- sub("^\xEF\xBB\xBF","",names(d),useBytes=TRUE)
 d
}
save_csv <- function(value,name) write.csv(value,file.path(root,name),row.names=FALSE,na="NA")
assert <- function(ok,message) if(!isTRUE(ok)) stop(paste("STOP:",message))
candidates <- read_csv(file.path(root,"bulk/candidates.csv"))
assert(nrow(candidates)==890 && !anyDuplicated(candidates$gene_symbol),"candidate identity")
symbols <- readLines(file.path(ev,"genes.txt"))
samples <- read_csv(file.path(ev,"samples.csv"))
patients <- c("patient_1","patient_5","patient_6","patient_7","patient_8")
assert(identical(samples$patient_id,rep(patients,each=2)),"fixed patient order")
assert(identical(samples$tissue,rep(c("Primary","Node"),5)),"fixed tissue identity")
con <- file(file.path(ev,"counts.float64"),"rb")
raw <- readBin(con,what="double",n=253010,size=8,endian="little"); close(con)
assert(length(raw)==253010 && all(is.finite(raw)) && all(raw>=0) && all(raw==floor(raw)),"raw integer counts")
counts <- matrix(raw,nrow=25301,ncol=10,byrow=TRUE,dimnames=list(symbols,samples$sample_id))
sample_data <- data.frame(patient=factor(samples$patient_id,levels=patients),tissue=factor(samples$tissue,levels=c("Primary","Node")))
design <- model.matrix(~patient+tissue,sample_data)
rownames(design) <- samples$sample_id
rank <- qr(design)$rank
assert(rank==6 && rank==ncol(design) && nrow(design)==10,"full paired design rank")
assert("tissueNode" %in% colnames(design),"target Node minus Primary")
save_csv(data.frame(sample_id=samples$sample_id,patient=samples$patient_id,tissue=samples$tissue,design,check.names=FALSE),"malignant_DESIGN_MATRIX.csv")
y <- DGEList(counts=counts)
y <- calcNormFactors(y,method="TMM")
raw_lib <- y$samples$lib.size
factors <- y$samples$norm.factors
effective_lib <- raw_lib*factors
assert(all(is.finite(factors)) && all(factors>0),"TMM factors")
normalization <- data.frame(sample_id=samples$sample_id,patient=samples$patient_id,tissue=samples$tissue,
 raw_library_size=raw_lib,TMM_factor=factors,effective_library_size=effective_lib,log_offset=log(effective_lib))
write.csv(normalization,file.path(ev,"SELECTED_TMM_AND_OFFSETS.csv"),row.names=FALSE)
# The scientific settings explicitly specifies raw full-background library sizes for filtering.
keep <- filterByExpr(y,design=design,lib.size=raw_lib,min.count=10,min.total.count=15,large.n=10,min.prop=0.7)
assert(length(keep)==25301 && is.logical(keep) && !anyNA(keep),"selected abundance filter")
save_csv(data.frame(feature_row_one_based=seq_along(symbols),gene_symbol=symbols,filter_pass=keep,
 filter_reason=ifelse(keep,"SELECTED_ABUNDANCE_FILTER_PASS","SELECTED_ABUNDANCE_FILTER_FAIL")),"malignant_BACKGROUND_FILTER_DETAILS.csv")
feature_rows <- match(candidates$gene_symbol,symbols)
transferred <- !is.na(feature_rows)
passes <- rep(FALSE,890); passes[transferred] <- keep[feature_rows[transferred]]
tested_indices <- which(passes)
filtered_rows <- which(keep)
fit_indices <- match(feature_rows[tested_indices],filtered_rows)
assert(length(tested_indices)>0 && !anyNA(fit_indices),"candidate-only test indices")
write.csv(data.frame(gene_symbol=candidates$gene_symbol,transferred=transferred,filter_pass=passes,
 planned_test=passes),file.path(ev,"PRE_FIT_CANDIDATE_EVALUABILITY.csv"),row.names=FALSE)
yf <- y[keep,,keep.lib.sizes=TRUE]
assert(identical(yf$samples$lib.size,raw_lib) && identical(yf$samples$norm.factors,factors),"normalization held fixed")
cat("Filtered background:",sum(keep),"; candidate-only family:",length(tested_indices),"\n")
fit <- glmQLFit(yf,design=design,legacy=FALSE,abundance.trend=TRUE,robust=FALSE,
 winsor.tail.p=c(0.05,0.10),dispersion=NULL,top.proportion=NULL)
candidate_fit <- fit[fit_indices, ]
# Prove that the installed row subset preserves every fitted nuisance quantity used by QLF.
row_fields <- c("dispersion","prior.n","prior.df","s2.post","s2.prior","df.prior","df.residual",
 "df.residual.zeros","df.residual.adj","deviance","deviance.adj","iter","failed","AveLogCPM")
checks <- list()
for(n in row_fields) {
 original <- fit[[n]]; subvalue <- candidate_fit[[n]]
 if(is.null(original)) {
  ok <- is.null(subvalue)
 } else {
  expected <- if(length(original)==nrow(fit$counts)) original[fit_indices] else original
  ok <- identical(unname(subvalue),unname(expected))
 }
 checks[[n]] <- data.frame(quantity=n,exactly_preserved=ok)
}
for(n in c("counts","offset","coefficients","fitted.values","unshrunk.coefficients","leverage","unit.deviance.adj","unit.df.adj")) {
 original <- fit[[n]]; subvalue <- candidate_fit[[n]]
 ok <- if(is.null(original)) is.null(subvalue) else identical(unname(as.matrix(subvalue)),unname(as.matrix(original)[fit_indices,,drop=FALSE]))
 checks[[n]] <- data.frame(quantity=n,exactly_preserved=ok)
}
for(n in c("design","samples","prior.count")) {
 checks[[n]] <- data.frame(quantity=n,exactly_preserved=identical(candidate_fit[[n]],fit[[n]]))
}
preservation <- do.call(rbind,checks)
write.csv(preservation,file.path(ev,"FITTED_OBJECT_SUBSET_PRESERVATION.csv"),row.names=FALSE)
assert(all(preservation$exactly_preserved),"fitted nuisance quantities altered by subsetting")
assert(!is.null(fit$df.residual.adj),"nonlegacy adjusted residual df present")
# glmQLFTest caps total df at sum(df.residual); candidate-only subsetting must not change effective df.
full_df <- pmin(fit$df.prior+fit$df.residual.adj,sum(fit$df.residual))[fit_indices]
candidate_df <- pmin(candidate_fit$df.prior+candidate_fit$df.residual.adj,sum(candidate_fit$df.residual))
df_equal <- identical(unname(full_df),unname(candidate_df))
write.csv(data.frame(gene_symbol=candidates$gene_symbol[tested_indices],full_background_effective_df=full_df,
 candidate_subset_effective_df=candidate_df,equal=full_df==candidate_df),file.path(ev,"QLF_DF_BOUNDARY_CHECK.csv"),row.names=FALSE)
assert(df_equal,"candidate subset changes effective QLF degrees of freedom")
cat("Candidate-only QLF nuisance and df eligibility PASS\n")
ql <- glmQLFTest(candidate_fit,coef="tissueNode",poisson.bound=FALSE)
table <- ql$table
assert(nrow(table)==length(tested_indices) && identical(rownames(table),candidates$gene_symbol[tested_indices]),"candidate QLF output identity")
assert(isTRUE(all.equal(unname(ql$df.total),unname(candidate_df),tolerance=0)),"actual test df")
q <- p.adjust(table$PValue,method="BH",n=length(tested_indices))
failed <- if(is.null(candidate_fit$failed)) rep(FALSE,length(tested_indices)) else as.logical(candidate_fit$failed)
assert(length(failed)==length(tested_indices) && !anyNA(failed),"convergence flags")
evaluable_subset <- is.finite(table$logFC) & is.finite(table$PValue) & is.finite(q) & !failed
logcpm <- cpm(counts[feature_rows[tested_indices],,drop=FALSE],lib.size=effective_lib,log=TRUE,prior.count=2)
primary_cols <- match(paste0(patients,"_Primary"),samples$sample_id)
node_cols <- match(paste0(patients,"_Node"),samples$sample_id)
deltas <- logcpm[,node_cols,drop=FALSE]-logcpm[,primary_cols,drop=FALSE]
colnames(deltas) <- patients
assert(all(is.finite(deltas)),"paired descriptive effects finite")
grid <- as.matrix(expand.grid(rep(list(1:5),5)))
assert(nrow(grid)==3125,"all ordered bootstrap samples")
weights <- t(apply(grid,1,function(x) tabulate(x,nbins=5)/5))
bootstrap_means <- deltas %*% t(weights)
interval <- t(apply(bootstrap_means,1,quantile,probs=c(0.025,0.975),type=8,names=FALSE))
result <- data.frame(gene_symbol=candidates$gene_symbol,bulk_logFC=candidates$logFC,bulk_direction=candidates$bulk_direction,
 bulk_q=candidates$adj.P.Val,transferred=transferred,transferred_feature_row=feature_rows,filter_pass=passes,
 malignant_tested=passes,malignant_logFC=NA_real_,malignant_logCPM=NA_real_,malignant_F=NA_real_,malignant_P=NA_real_,malignant_q=NA_real_,
 malignant_direction=NA_character_,direction_compatible=NA,mean_patient_log2CPM_delta=NA_real_,
 descriptive_interval_low=NA_real_,descriptive_interval_high=NA_real_,support_class="S5",support_class_reason="NOT_TRANSFERRED_EXACT_SYMBOL")
result$support_class_reason[transferred & !passes] <- "SELECTED_ABUNDANCE_FILTER_FAIL"
result$malignant_logFC[tested_indices] <- table$logFC
result$malignant_logCPM[tested_indices] <- table$logCPM
result$malignant_F[tested_indices] <- table$F
result$malignant_P[tested_indices] <- table$PValue
result$malignant_q[tested_indices] <- q
result$mean_patient_log2CPM_delta[tested_indices] <- rowMeans(deltas)
result$descriptive_interval_low[tested_indices] <- interval[,1]
result$descriptive_interval_high[tested_indices] <- interval[,2]
result$malignant_direction[tested_indices] <- ifelse(table$logFC>0,"Node_up",ifelse(table$logFC<0,"Node_down","Zero"))
for(k in seq_along(tested_indices)) {
 i <- tested_indices[k]
 if(!evaluable_subset[k]) {
  result$support_class_reason[i] <- if(failed[k]) "NONCONVERGENT_FIT" else "NONFINITE_COEFFICIENT_OR_P_Q"
  next
 }
 sign_m <- sign(table$logFC[k]); sign_b <- sign(result$bulk_logFC[i])
 compatible <- sign_m!=0 && sign_m==sign_b
 result$direction_compatible[i] <- compatible
 if(sign_m!=0 && sign_m!=sign_b) {
  result$support_class[i] <- "S3"; result$support_class_reason[i] <- "OPPOSITE_NONZERO_EDGER_DIRECTION"
 } else if(compatible && q[k]<0.05) {
  result$support_class[i] <- "S1"; result$support_class_reason[i] <- "SAME_NONZERO_DIRECTION_AND_CANDIDATE_BH_Q_LT_0.05"
 } else if(sign_m==0) {
  result$support_class[i] <- "S4"; result$support_class_reason[i] <- "EXACT_ZERO_EDGER_COEFFICIENT"
 } else if(compatible && q[k]>=0.05 && interval[k,1]<=-1 && interval[k,2]>=1) {
  result$support_class[i] <- "S4"; result$support_class_reason[i] <- "SAME_DIRECTION_Q_GE_0.05_INTERVAL_SPANS_MINUS1_PLUS1"
 } else {
  assert(compatible && q[k]>=0.05,"exhaustive S2 precedence")
  result$support_class[i] <- "S2"; result$support_class_reason[i] <- "SAME_NONZERO_DIRECTION_Q_GE_0.05_NOT_S4_PRECISION_FLAG"
 }
}
assert(nrow(result)==890 && all(result$support_class %in% paste0("S",1:5)),"all selected candidates exhaustively classified")
save_csv(result,"bulk_malignant_ALL_CANDIDATES.csv")
all_deltas <- matrix(NA_real_,890,5,dimnames=list(candidates$gene_symbol,patients))
all_deltas[tested_indices,] <- deltas
effects <- data.frame(gene_symbol=candidates$gene_symbol,malignant_tested=passes,not_evaluable_reason=ifelse(result$support_class=="S5",result$support_class_reason,NA_character_),
 all_deltas,mean_delta=rowMeans(all_deltas),check.names=FALSE)
save_csv(effects,"bulk_malignant_PATIENT_EFFECTS.csv")
sets <- list(`30`=c("patient_1","patient_5","patient_7","patient_8"),`40`=c("patient_1","patient_5","patient_8"),`50`=c("patient_5","patient_8"))
sensitivities <- lapply(names(sets),function(threshold) {
 ids <- sets[[threshold]]; values <- all_deltas[,ids,drop=FALSE]
 mean <- rowMeans(values); median <- apply(values,1,median)
 data.frame(gene_symbol=candidates$gene_symbol,threshold=as.integer(threshold),n_patients=length(ids),
 retained_patients=paste(ids,collapse=";"),bulk_direction=candidates$bulk_direction,malignant_primary_tested=passes,
 not_evaluable_reason=ifelse(result$support_class=="S5",result$support_class_reason,NA_character_),
 patient_1=if("patient_1" %in% ids) all_deltas[,"patient_1"] else NA_real_,
 patient_5=if("patient_5" %in% ids) all_deltas[,"patient_5"] else NA_real_,
 patient_6=NA_real_,patient_7=if("patient_7" %in% ids) all_deltas[,"patient_7"] else NA_real_,
 patient_8=all_deltas[,"patient_8"],mean_delta=mean,median_delta=median,
 min_delta=apply(values,1,min),max_delta=apply(values,1,max),mean_direction=ifelse(is.na(mean),NA_character_,ifelse(mean>0,"Node_up",ifelse(mean<0,"Node_down","Zero"))),
 bulk_compatible=ifelse(is.na(mean),NA,sign(mean)!=0 & sign(mean)==sign(candidates$logFC)),
 inferential_refit=FALSE,normalization_reestimated=FALSE)
})
sens <- do.call(rbind,sensitivities)
assert(nrow(sens)==2670,"complete nested sensitivity rows")
save_csv(sens,"bulk_malignant_THRESHOLD_SENSITIVITY.csv")
lopo_means <- sapply(patients,function(omit) rowMeans(all_deltas[,setdiff(patients,omit),drop=FALSE]))
colnames(lopo_means) <- paste0("omit_",patients)
full_mean <- rowMeans(all_deltas)
lopo <- data.frame(gene_symbol=candidates$gene_symbol,malignant_primary_tested=passes,
 not_evaluable_reason=ifelse(result$support_class=="S5",result$support_class_reason,NA_character_),
 full_five_mean_log2CPM_delta=full_mean,full_five_EDGER_logFC=result$malignant_logFC,
 number_of_5_omissions_retaining_full_five_malignant_direction=rowSums(sign(lopo_means)==sign(result$malignant_logFC) & sign(result$malignant_logFC)!=0),
 number_of_5_omissions_retaining_full_five_descriptive_mean_direction=rowSums(sign(lopo_means)==sign(full_mean) & sign(full_mean)!=0),
 min_mean_effect=apply(lopo_means,1,min),median_mean_effect=apply(lopo_means,1,median),max_mean_effect=apply(lopo_means,1,max),
 lopo_means,inferential_refit=FALSE,normalization_reestimated=FALSE,check.names=FALSE)
save_csv(lopo,"bulk_malignant_LOPO_DESCRIPTIVE.csv")
class_counts <- setNames(as.integer(table(factor(result$support_class,levels=paste0("S",1:5)))),paste0("S",1:5))
summary_rows <- do.call(rbind,lapply(c("ALL","Node_up","Node_down"),function(group) {
 r <- if(group=="ALL") result else result[result$bulk_direction==group,]
 n <- sum(r$support_class!="S5"); compatible <- sum(r$direction_compatible %in% TRUE)
 data.frame(bulk_direction_group=group,support_class=paste0("S",1:5),
 N=as.integer(table(factor(r$support_class,levels=paste0("S",1:5)))),selected_candidates_N=nrow(r),
 transferred_N=sum(r$transferred),tested_N=sum(r$malignant_tested),evaluable_N=n,direction_compatible_N=compatible,
 support_fraction=as.integer(table(factor(r$support_class,levels=paste0("S",1:5))))/nrow(r))
}))
save_csv(summary_rows,"malignant_SUPPORT_CLASS_COUNTS.csv")
