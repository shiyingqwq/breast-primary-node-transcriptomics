options(stringsAsFactors=FALSE)
root <- normalizePath(Sys.getenv("BCPN_WORKSPACE"),winslash="/",mustWork=TRUE)
if(nzchar(Sys.getenv("BCPN_R_LIBRARY"))) .libPaths(c(Sys.getenv("BCPN_R_LIBRARY"),.libPaths()))
library(edgeR)
ev <- file.path(root,"cellular/nonmalignant")
read_csv <- function(path) {
 d <- read.csv(path,check.names=FALSE,encoding="UTF-8")
 names(d) <- sub("^\xEF\xBB\xBF","",names(d),useBytes=TRUE)
 d
}
assert <- function(ok,message) if(!isTRUE(ok)) stop(paste("STOP:",message))
labels <- c("B","CAFs","CD4+ T","CD8+ T","DCs","Macrophage","Mast cells","NK","TECs","plasma B")
patients <- c("patient_1","patient_5","patient_6","patient_7","patient_8")
bulk <- read_csv(file.path(root,"bulk/candidates.csv"))
assert(nrow(bulk)==890 && !anyDuplicated(bulk$gene_symbol),"full bulk candidate set")
symbols <- readLines(file.path(ev,"genes.txt"))
columns <- read_csv(file.path(ev,"samples.csv"))
stopifnot(!anyNA(columns$author_label_cell_count),all(columns$author_label_cell_count>=20))
release <- list(shape=c(length(symbols),nrow(columns)),released_population_names=unique(columns$population))
eligibility <- list(eligible_population_names=unique(columns$population))
assert(nrow(columns)==release$shape[2] && length(symbols)==25301,"complete raw background")
con <- file(file.path(ev,"counts.float64"),"rb")
raw <- readBin(con,what="double",n=prod(release$shape),size=8,endian="little"); close(con)
assert(length(raw)==prod(release$shape) && all(is.finite(raw)) && all(raw>=0) && all(raw==floor(raw)),"raw integer pseudobulk")
all_counts <- matrix(raw,nrow=25301,ncol=nrow(columns),byrow=TRUE,dimnames=list(symbols,NULL))
feature_rows <- match(bulk$gene_symbol,symbols)
transferred <- !is.na(feature_rows)
bootstrap_grid <- as.matrix(expand.grid(rep(list(1:5),5)))
weights <- t(apply(bootstrap_grid,1,function(x) tabulate(x,nbins=5)/5))
assert(nrow(weights)==3125,"selected descriptive bootstrap draws")
results <- model_summaries <- normalization_list <- background_list <- nuisance_list <- df_list <- list()
model_calls <- 0L
for(label in labels) {
 result <- data.frame(gene_symbol=bulk$gene_symbol,population=label,bulk_direction=bulk$bulk_direction,
  bulk_logFC=bulk$logFC,bulk_q=bulk$adj.P.Val,candidate_transferred=transferred,transferred_feature_row=feature_rows,
  population_metadata_eligible=label %in% eligibility$eligible_population_names,
  population_technical_released=label %in% release$released_population_names,population_filter_pass=NA,
  tested=FALSE,statistically_evaluable=FALSE,nonmalignant_logFC=NA_real_,nonmalignant_logCPM=NA_real_,
  nonmalignant_F=NA_real_,raw_P=NA_real_,joint_BH_q=NA_real_,nonmalignant_direction=NA_character_,
  direction_compatible_with_bulk=NA,N_status="not_evaluable",not_evaluable_reason="POPULATION_METADATA_COUNT_ELIGIBILITY_FAIL",
  mean_patient_log2CPM_delta=NA_real_,descriptive_interval_low=NA_real_,descriptive_interval_high=NA_real_,
  patient_1_delta=NA_real_,patient_5_delta=NA_real_,patient_6_delta=NA_real_,patient_7_delta=NA_real_,patient_8_delta=NA_real_)
 if(!label %in% release$released_population_names) {
  if(label %in% eligibility$eligible_population_names) result$not_evaluable_reason <- "POPULATION_ZERO_LIBRARY_TECHNICAL_FAIL"
  results[[label]] <- result
  model_summaries[[label]] <- list(population=label,released=FALSE,tested_N=0,reason=unique(result$not_evaluable_reason))
  next
 }
 indices <- which(columns$population==label)
 s <- columns[indices,]
 assert(length(indices)==10 && identical(s$patient_id,rep(patients,each=2)) &&
  identical(s$tissue,rep(c("Primary","Node"),5)),paste(label,"sample identities"))
 counts <- all_counts[,indices,drop=FALSE]; colnames(counts) <- s$sample_id
 assert(all(colSums(counts)>0),paste(label,"nonzero libraries"))
 design <- model.matrix(~patient+tissue,data.frame(patient=factor(s$patient_id,levels=patients),tissue=factor(s$tissue,levels=c("Primary","Node"))))
 rownames(design) <- s$sample_id
 rank <- qr(design)$rank
 assert(rank==6 && rank==ncol(design),paste(label,"full paired rank"))
 safe_name <- sprintf("population_%02d",match(label,labels))
 write.csv(data.frame(sample_id=s$sample_id,patient=s$patient_id,tissue=s$tissue,design,check.names=FALSE),file.path(ev,paste0(safe_name,"_DESIGN_MATRIX.csv")),row.names=FALSE)
 y <- calcNormFactors(DGEList(counts=counts),method="TMM")
 raw_lib <- y$samples$lib.size; factors <- y$samples$norm.factors; effective <- raw_lib*factors
 normalization_list[[label]] <- data.frame(population=label,sample_id=s$sample_id,patient=s$patient_id,tissue=s$tissue,
  raw_library_size=raw_lib,TMM_factor=factors,effective_library_size=effective,log_offset=log(effective))
 keep <- filterByExpr(y,design=design,lib.size=raw_lib,min.count=10,min.total.count=15,large.n=10,min.prop=0.7)
 assert(length(keep)==25301 && !anyNA(keep) && sum(keep)>0,paste(label,"full-background filter"))
 background_list[[label]] <- data.frame(population=label,feature_row_one_based=seq_along(symbols),gene_symbol=symbols,filter_pass=keep)
 passes <- rep(FALSE,890); passes[transferred] <- keep[feature_rows[transferred]]
 test_indices <- which(passes)
 result$population_filter_pass <- passes
 result$not_evaluable_reason <- ifelse(transferred,"POPULATION_SELECTED_ABUNDANCE_FILTER_FAIL","NOT_TRANSFERRED_EXACT_SYMBOL")
 planned <- result[,c("gene_symbol","population","candidate_transferred","population_filter_pass")]
 write.csv(planned,file.path(ev,paste0(safe_name,"_PRE_FIT_CANDIDATE_EVALUABILITY.csv")),row.names=FALSE)
 if(length(test_indices)==0) {
  results[[label]] <- result
  model_summaries[[label]] <- list(population=label,released=TRUE,background_N=25301,filter_pass_N=sum(keep),
   design_rank=rank,transferred_N=sum(transferred),candidate_filter_pass_N=0,tested_N=0,fit_performed=FALSE)
  next
 }
 yf <- y[keep,,keep.lib.sizes=TRUE]
 assert(identical(yf$samples$lib.size,raw_lib) && identical(yf$samples$norm.factors,factors),"fixed per-population TMM")
 fit_indices <- match(feature_rows[test_indices],which(keep)); assert(!anyNA(fit_indices),"candidate rows")
 fit <- glmQLFit(yf,design=design,legacy=FALSE,abundance.trend=TRUE,robust=FALSE,
  winsor.tail.p=c(0.05,0.10),dispersion=NULL,top.proportion=NULL)
 model_calls <- model_calls+1L
 subfit <- fit[fit_indices, ]
 preservation <- list()
 for(n in c("dispersion","prior.n","prior.df","s2.post","s2.prior","df.prior","df.residual",
  "df.residual.zeros","df.residual.adj","deviance","deviance.adj","iter","failed","AveLogCPM")) {
  original <- fit[[n]]
  expected <- if(!is.null(original) && length(original)==nrow(fit$counts)) original[fit_indices] else original
  preservation[[n]] <- data.frame(population=label,quantity=n,exactly_preserved=identical(unname(subfit[[n]]),unname(expected)))
 }
 for(n in c("counts","offset","coefficients","fitted.values","unshrunk.coefficients","leverage","unit.deviance.adj","unit.df.adj")) {
  original <- fit[[n]]
  ok <- if(is.null(original)) is.null(subfit[[n]]) else identical(unname(as.matrix(subfit[[n]])),unname(as.matrix(original)[fit_indices,,drop=FALSE]))
  preservation[[n]] <- data.frame(population=label,quantity=n,exactly_preserved=ok)
 }
 for(n in c("design","samples","prior.count")) preservation[[n]] <- data.frame(population=label,quantity=n,exactly_preserved=identical(subfit[[n]],fit[[n]]))
 check <- do.call(rbind,preservation); nuisance_list[[label]] <- check
 write.csv(check,file.path(ev,paste0(safe_name,"_NUISANCE_SUBSET_CHECK.csv")),row.names=FALSE)
 assert(all(check$exactly_preserved),paste(label,"fitted nuisance preservation"))
 assert(!is.null(fit$df.residual.adj),"nonlegacy QL fitted object")
 full_df <- pmin(fit$df.prior+fit$df.residual.adj,sum(fit$df.residual))[fit_indices]
 subset_df <- pmin(subfit$df.prior+subfit$df.residual.adj,sum(subfit$df.residual))
 assert(identical(unname(full_df),unname(subset_df)),paste(label,"candidate subset changes effective QLF df"))
 df_list[[label]] <- data.frame(population=label,gene_symbol=bulk$gene_symbol[test_indices],
  full_background_effective_df=full_df,candidate_subset_effective_df=subset_df,equal=full_df==subset_df)
 qlf <- glmQLFTest(subfit,coef="tissueNode",poisson.bound=FALSE)
 t <- qlf$table
 assert(nrow(t)==length(test_indices) && identical(rownames(t),bulk$gene_symbol[test_indices]),"candidate-only output identity")
 assert(isTRUE(all.equal(unname(qlf$df.total),unname(subset_df),tolerance=0)),"QLF actual effective df")
 failed <- if(is.null(subfit$failed)) rep(FALSE,length(test_indices)) else as.logical(subfit$failed)
 assert(length(failed)==length(test_indices) && !anyNA(failed),"fit convergence flags")
 result$tested[test_indices] <- TRUE
 result$statistically_evaluable[test_indices] <- is.finite(t$logFC) & is.finite(t$PValue) & !failed
 result$nonmalignant_logFC[test_indices] <- t$logFC
 result$nonmalignant_logCPM[test_indices] <- t$logCPM
 result$nonmalignant_F[test_indices] <- t$F
 result$raw_P[test_indices] <- t$PValue
 result$not_evaluable_reason[test_indices] <- ifelse(failed,"NONCONVERGENT_FIT",ifelse(result$statistically_evaluable[test_indices],NA_character_,"NONFINITE_COEFFICIENT_OR_P"))
 logcpm <- cpm(counts[feature_rows[test_indices],,drop=FALSE],lib.size=effective,log=TRUE,prior.count=2)
 deltas <- logcpm[,match(paste0(patients,"_Node"),s$sample_id),drop=FALSE]-logcpm[,match(paste0(patients,"_Primary"),s$sample_id),drop=FALSE]
 assert(all(is.finite(deltas)),"paired descriptive effects")
 intervals <- t(apply(deltas %*% t(weights),1,quantile,probs=c(.025,.975),type=8,names=FALSE))
 result$mean_patient_log2CPM_delta[test_indices] <- rowMeans(deltas)
 result$descriptive_interval_low[test_indices] <- intervals[,1]
 result$descriptive_interval_high[test_indices] <- intervals[,2]
 result[test_indices,paste0(patients,"_delta")] <- deltas
 results[[label]] <- result
 model_summaries[[label]] <- list(population=label,released=TRUE,background_N=25301,filter_pass_N=sum(keep),
  filter_fail_N=sum(!keep),design_rank=rank,design_columns=colnames(design),target="tissueNode",n_patients=5,n_samples=10,
  transferred_N=sum(transferred),nontransferred_N=sum(!transferred),candidate_filter_pass_N=length(test_indices),
  candidate_filter_fail_N=sum(transferred & !passes),tested_N=length(test_indices),fit_performed=TRUE,
  subset_nuisance_preserved=TRUE,effective_QLF_df_unchanged=TRUE,full_background_df_cap=sum(fit$df.residual),
  candidate_subset_df_cap=sum(subfit$df.residual),candidate_effective_df_range=range(subset_df))
 cat("POPULATION candidate-only QLF PASS:",label,"; background",sum(keep),"; tests",length(test_indices),"\n")
}
combined <- do.call(rbind,results)
assert(nrow(combined)==8900 && !anyDuplicated(paste(combined$population,combined$gene_symbol)),"complete ten-population denominator")
tested_rows <- which(combined$tested)
# The only BH call in this stage spans all actual candidate x population tests together.
combined$joint_BH_q[tested_rows] <- p.adjust(combined$raw_P[tested_rows],method="BH",n=length(tested_rows))
good <- combined$statistically_evaluable & is.finite(combined$joint_BH_q)
combined$statistically_evaluable <- good
combined$nonmalignant_direction[good] <- ifelse(combined$nonmalignant_logFC[good]>0,"Node_up",ifelse(combined$nonmalignant_logFC[good]<0,"Node_down","Zero"))
combined$direction_compatible_with_bulk[good] <- sign(combined$nonmalignant_logFC[good])!=0 & sign(combined$nonmalignant_logFC[good])==sign(combined$bulk_logFC[good])
combined$N_status[good] <- ifelse(combined$direction_compatible_with_bulk[good] & combined$joint_BH_q[good]<.05,"N+","N0")
combined$not_evaluable_reason[good] <- NA_character_
# Malignant results are attached only after the full independent test family and adjustment are complete.
malignant <- read_csv(file.path(root,"bulk_malignant_ALL_CANDIDATES.csv"))
assert(identical(malignant$gene_symbol,bulk$gene_symbol),"immutable malignant identities")
combined$malignant_S_class <- malignant$support_class[match(combined$gene_symbol,malignant$gene_symbol)]
write.csv(combined,file.path(root,"bulk_NONMALIGNANT_ALL_CANDIDATES.csv"),row.names=FALSE,na="NA")
positive <- combined[combined$N_status=="N+",]
counts_by_pop <- do.call(rbind,lapply(labels,function(label) {
 r <- combined[combined$population==label,]
 data.frame(population=label,metadata_eligible=label %in% eligibility$eligible_population_names,technical_released=label %in% release$released_population_names,
  candidate_rows=nrow(r),transferred_N=sum(r$candidate_transferred),filter_pass_N=sum(r$population_filter_pass %in% TRUE),
  tested_N=sum(r$tested),evaluable_N=sum(r$statistically_evaluable),Nplus_N=sum(r$N_status=="N+"),N0_N=sum(r$N_status=="N0"),
  not_evaluable_N=sum(r$N_status=="not_evaluable"),joint_BH_family_N=length(tested_rows))
}))
write.csv(counts_by_pop,file.path(root,"nonmalignant_NPLUS_COUNTS_BY_POPULATION.csv"),row.names=FALSE)
s3 <- malignant[malignant$support_class=="S3",]
assert(nrow(s3)==229,"selected S3 denominator")
s3_summary <- do.call(rbind,lapply(s3$gene_symbol,function(gene) {
 names <- positive$population[positive$gene_symbol==gene]
 data.frame(gene=gene,S3=TRUE,number_Nplus_populations=length(names),Nplus_population_names=paste(names,collapse=";"),
  class_B_eligible=if(length(names)>0) "YES" else "NO")
}))
write.csv(s3_summary,file.path(root,"nonmalignant_S3_NONMALIGNANT_ATTRIBUTION_SUMMARY.csv"),row.names=FALSE)
write.csv(do.call(rbind,normalization_list),file.path(ev,"PER_POPULATION_TMM_AND_OFFSETS.csv"),row.names=FALSE)
write.csv(do.call(rbind,background_list),file.path(ev,"PER_POPULATION_BACKGROUND_FILTER_DETAILS.csv"),row.names=FALSE)
write.csv(do.call(rbind,nuisance_list),file.path(ev,"ALL_NUISANCE_SUBSET_CHECKS.csv"),row.names=FALSE)
write.csv(do.call(rbind,df_list),file.path(ev,"ALL_QLF_DF_BOUNDARY_CHECKS.csv"),row.names=FALSE)
