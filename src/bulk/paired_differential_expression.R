# Paired primary-versus-node analysis of author-processed expression.
options(stringsAsFactors=FALSE,contrasts=c("contr.treatment","contr.poly"))
Sys.setlocale("LC_COLLATE","C")
root <- normalizePath(Sys.getenv("BCPN_WORKSPACE"),winslash="/",mustWork=TRUE)
out <- file.path(root,"bulk");dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(nzchar(Sys.getenv("BCPN_R_LIBRARY"))) .libPaths(c(Sys.getenv("BCPN_R_LIBRARY"),.libPaths()))
read_csv <- function(path)read.csv(path,check.names=FALSE,fileEncoding="UTF-8-BOM")
write_new_csv <- function(object,path) {
 if(file.exists(path))stop("Output already exists: ",path)
 write.csv(object,path,row.names=FALSE,fileEncoding="UTF-8",na="NA")
}
samples <- read_csv(file.path(root,"inputs/bulk/samples.csv"))
mapping <- read_csv(file.path(root,"inputs/bulk/probe_mapping.csv"))
patients <- unique(samples$patient)
u1 <- readLines(file.path(root,"inputs/bulk/genes.txt"))
mapping <- mapping[match(u1,mapping$gene),,drop=FALSE]
stopifnot(nrow(samples)==72,length(patients)==36,!anyDuplicated(samples$GSM),
 length(u1)==12502,!anyDuplicated(u1),!anyDuplicated(mapping$selected_probe))
con <- gzfile(file.path(root,"inputs/bulk/GSE57968_expression.txt.gz"),"rt")
matrix_table <- read.delim(con,header=TRUE,check.names=FALSE,comment.char="!",quote="\"")
close(con)
probe_index <- match(mapping$selected_probe,matrix_table[[1]])
sample_index <- match(samples$GSM,names(matrix_table))
stopifnot(!anyNA(probe_index),!anyNA(sample_index))
x <- as.matrix(matrix_table[probe_index,sample_index,drop=FALSE]);storage.mode(x)<-"double"
rownames(x)<-u1;colnames(x)<-samples$GSM
stopifnot(identical(dim(x),c(12502L,72L)),all(is.finite(x)))
samples$patient <- factor(samples$patient,levels=patients)
samples$tissue <- factor(samples$tissue,levels=c("Primary","Node"))
stopifnot(!anyNA(samples$patient),!anyNA(samples$tissue),all(table(samples$patient,samples$tissue)==1L))
design <- model.matrix(~patient+tissue,samples);rownames(design)<-samples$GSM
stopifnot(qr(design)$rank==37L,nrow(design)==72L,"tissueNode"%in%colnames(design))
rank <- qr(design)$rank
fit <- limma::lmFit(x,design)
fit <- limma::eBayes(fit,robust=FALSE,trend=FALSE)
coef_index <- match("tissueNode",colnames(fit$coefficients))
tab <- limma::topTable(fit,coef=coef_index,number=Inf,sort.by="none",adjust.method="BH",confint=TRUE)
stopifnot(nrow(tab)==12502,identical(rownames(tab),u1),all(is.finite(tab$P.Value)),
 identical(unname(tab$adj.P.Val),unname(p.adjust(tab$P.Value,method="BH",n=12502))))
full <- data.frame(gene_symbol=u1,probe_id=mapping$selected_probe,logFC=tab$logFC,AveExpr=tab$AveExpr,t=tab$t,
 P.Value=tab$P.Value,adj.P.Val=tab$adj.P.Val,B=tab$B,
 ordinary_SE=fit$stdev.unscaled[,coef_index]*fit$sigma,
 moderated_SE=fit$stdev.unscaled[,coef_index]*sqrt(fit$s2.post),CI.L=tab$CI.L,CI.R=tab$CI.R)
stopifnot(all(is.finite(as.matrix(full[,-c(1,2)]))))
write_new_csv(full,file.path(root,"bulk/all_genes.csv"))
candidates <- full[full$adj.P.Val<0.05,,drop=FALSE]
stopifnot(all(candidates$logFC!=0))
candidates$bulk_direction <- ifelse(candidates$logFC>0,"Node_up","Node_down")
write_new_csv(candidates,file.path(root,"bulk/candidates.csv"))
txt_path <- file.path(root,"bulk/candidates.txt")
stopifnot(!file.exists(txt_path))
writeLines(candidates$gene_symbol,txt_path,useBytes=TRUE)
warnings <- character()
hallmark_performed <- FALSE
hallmark_sets <- 0L
hallmark_nonestimable <- 0L
{
 pathways <- fgsea::gmtPathways(file.path(root,'resources/h.all.v2025.1.Hs.symbols.gmt'))
 stopifnot(length(pathways)==50L,!anyDuplicated(names(pathways)))
 pathways <- lapply(pathways,function(genes) intersect(unique(genes),u1))
 sizes <- lengths(pathways)
 eligible <- sizes>=15 & sizes<=500
 write_new_csv(data.frame(pathway=names(pathways),size_in_U1=sizes,evaluable=eligible),file.path(out,"HALLMARK_SET_ELIGIBILITY.csv"))
 stats <- setNames(full$t,full$gene_symbol)
 stats <- stats[order(-stats,names(stats),method="radix")]
 stopifnot(length(stats)==12502,!anyDuplicated(names(stats)),all(is.finite(stats)))
 write_new_csv(data.frame(gene_symbol=names(stats),moderated_t=unname(stats)),file.path(out,"COMPLETE_HALLMARK_RANKING.csv"))
 set.seed(20261005)
 gsea <- withCallingHandlers(fgsea::fgseaMultilevel(pathways=pathways[eligible],stats=stats,
   minSize=15,maxSize=500,eps=0,scoreType="std",gseaParam=1,sampleSize=101,nPermSimple=1000,nproc=1),
   warning=function(w) { warnings <<- c(warnings,conditionMessage(w)); cat("FGSEA WARNING:",conditionMessage(w),"\n"); invokeRestart("muffleWarning") })
 gsea <- as.data.frame(gsea)
 stopifnot(nrow(gsea)==sum(eligible),setequal(gsea$pathway,names(pathways)[eligible]))
 gsea <- gsea[match(names(pathways)[eligible],gsea$pathway),,drop=FALSE]
 gsea$padj <- p.adjust(gsea$pval,method="BH",n=nrow(gsea))
 gsea$leadingEdge <- vapply(gsea$leadingEdge,paste,collapse=";",FUN.VALUE=character(1))
 gsea$estimability <- ifelse(is.na(gsea$pval),"P_NOT_ESTIMABLE_UNDER_SELECTED_METHOD","ESTIMABLE")
 write_new_csv(gsea,file.path(root,"bulk_HALLMARK_ALL_TESTED.csv"))
 hallmark_performed <- TRUE
 hallmark_sets <- nrow(gsea)
 hallmark_nonestimable <- sum(is.na(gsea$pval))
 cat("SECONDARY HALLMARK COMPLETE; ATTEMPTED SETS=",hallmark_sets,"\n",sep="")
}
nc <- nrow(candidates)
lopo <- candidates[,c("gene_symbol","probe_id","logFC","bulk_direction"),drop=FALSE]
names(lopo)[names(lopo)=="logFC"] <- "full_data_logFC"
effects <- matrix(NA_real_,nrow=nc,ncol=36,dimnames=list(candidates$gene_symbol,paste0("LOPO_omit_",patients)))
lopo_design_checks <- data.frame(omitted_patient=character(),remaining_pairs=integer(),samples=integer(),design_rank=integer(),design_columns=integer(),target_coefficient=character())
lopo_refits <- 0L
if(nc>0) {
 for(i in seq_along(patients)) {
   keep <- samples$patient!=patients[i]
   subset <- samples[keep,,drop=FALSE]
   subset$patient <- factor(as.character(subset$patient),levels=patients[-i])
   subset$tissue <- factor(as.character(subset$tissue),levels=c("Primary","Node"))
   design_i <- model.matrix(~ patient + tissue,data=subset)
   rank_i <- qr(design_i)$rank
   stopifnot(nrow(design_i)==70,ncol(design_i)==36,rank_i==36,"tissueNode" %in% colnames(design_i),is.null(limma::nonEstimable(design_i)))
   fit_i <- limma::lmFit(x[,keep,drop=FALSE],design_i)
   effects[,i] <- fit_i$coefficients[candidates$gene_symbol,"tissueNode"]
   lopo_refits <- lopo_refits+1L
   lopo_design_checks <- rbind(lopo_design_checks,data.frame(omitted_patient=patients[i],remaining_pairs=35L,samples=70L,design_rank=rank_i,design_columns=ncol(design_i),target_coefficient="tissueNode"))
   cat("LOPO completed ",i,"/36; omit ",patients[i],"\n",sep="")
 }
 stopifnot(all(is.finite(effects)))
 lopo$number_retaining_full_data_sign <- rowSums(sign(effects)==sign(candidates$logFC))
 lopo$proportion_retaining_full_data_sign <- lopo$number_retaining_full_data_sign/36
 lopo$LOPO_min <- apply(effects,1,min)
 lopo$LOPO_median <- apply(effects,1,median)
 lopo$LOPO_max <- apply(effects,1,max)
} else {
 lopo$number_retaining_full_data_sign <- integer()
 lopo$proportion_retaining_full_data_sign <- numeric()
 lopo$LOPO_min <- numeric(); lopo$LOPO_median <- numeric(); lopo$LOPO_max <- numeric()
}
lopo <- cbind(lopo,as.data.frame(effects))
write_new_csv(lopo,file.path(root,"bulk_BULK_LOPO_DESCRIPTIVE.csv"))
write_new_csv(lopo_design_checks,file.path(out,"lopo_design_checks.csv"))
