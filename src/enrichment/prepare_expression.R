source(file.path(Sys.getenv('BCPN_REPOSITORY'), 'src/enrichment/runtime.R'))
full <- read_csv(file.path(root,'bulk/all_genes.csv'))
cand <- read_csv(file.path(root,'bulk/candidates.csv'))
samples <- read_csv(file.path(root,'inputs/bulk/samples.csv'))
probes <- read_csv(file.path(root,'inputs/bulk/probe_mapping.csv'))
con <- gzfile(file.path(root,'inputs/bulk/GSE57968_expression.txt.gz'),'rt')
raw <- read.delim(con,check.names=FALSE,comment.char='!',quote='"');close(con)
tab <- table(samples$patient,samples$tissue)
stopifnot(nrow(samples)==72,!anyDuplicated(samples$GSM),nrow(tab)==36,all(tab==1),
          nrow(full)==12502,!anyDuplicated(full$gene_symbol),nrow(cand)==890,
          sum(cand$bulk_direction=='Node_up')==328,sum(cand$bulk_direction=='Node_down')==562,
          !anyDuplicated(cand$gene_symbol),setequal(cand$gene_symbol,full$gene_symbol[full$adj.P.Val<.05]))
probes <- probes[match(full$gene_symbol,probes$gene),]
ix <- match(probes$selected_probe,raw[[1]]);sx <- match(samples$GSM,names(raw))
stopifnot(!anyNA(ix),!anyNA(sx),all(probes$selected_probe==full$probe_id))
x <- as.matrix(raw[ix,sx]);storage.mode(x)<-'double'
rownames(x)<-full$gene_symbol;colnames(x)<-samples$GSM
stopifnot(all(is.finite(x)),max(abs(rowMeans(x)-full$AveExpr))<1e-8)
write_csv(samples,'work/SAMPLE_AND_PAIR_MAP.csv')

mapping <- read_csv(file.path(out,'work/GENE_ID_MAPPING_DETAILS.csv'))
stopifnot(identical(mapping$gene_symbol,full$gene_symbol))
write_csv(mapping,'work/GENE_ID_MAPPING_DETAILS.csv')
saveRDS(list(x=x,samples=samples,full=full,candidates=cand,mapping=mapping),file.path(out,'work/INPUTS.rds'))
