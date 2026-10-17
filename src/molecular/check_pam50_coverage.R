source(file.path(Sys.getenv('BCPN_REPOSITORY'), 'src/molecular/runtime.R'))
input <- readRDS(file.path(root,'enrichment/work/INPUTS.rds'))
data('pam50.robust',package='genefu')
cat('MODEL DETAILS\n');str(pam50.robust,max.level=1)
model <- pam50.robust;map <- model$centroids.map
cat('CENTROID MAP\n');print(head(map))
ids <- as.character(map[,'EntrezGene.ID'])
ann <- input$mapping;usable <- ann$annotation_usable & !is.na(ann$entrez_id)
duplicate <- duplicated(ann$entrez_id[usable])|duplicated(ann$entrez_id[usable],fromLast=TRUE)
keep <- which(usable)[!duplicate];m <- match(ids,ann$entrez_id[keep]);i <- keep[m]
details <- data.frame(centroid_probe=as.character(map[,'probe']),EntrezGene.ID=ids,
 selected_gene=ann$gene_symbol[i],approved_symbol=ann$approved_symbol[i],selected_probe=input$full$probe_id[i],
 available=!is.na(i),mapping_basis='Inherited unambiguous selected HGNC Entrez mapping; no new probe selection')
stopifnot(nrow(details)==50,nrow(model$centroids)==50,identical(as.character(map[,'probe']),rownames(model$centroids)))
if(any(!is.na(i)))stopifnot(all(apply(input$x[i[!is.na(i)],,drop=FALSE],1,sd)>0))
write_csv(details,'work/PAM50_GENE_COVERAGE.csv')
save_json(list(total_genes=50,available=sum(details$available),missing_entrez=ids[!details$available],
 missing_centroid_probes=details$centroid_probe[!details$available],minimum=43,
 status=if(sum(details$available)==43)'PASS_RESEARCH_IMPLEMENTATION' else 'SKIP',
 model_std=model$std,model_correlation=model$method.cor,model_rescale_q=model$rescale.q,
 implementation='Released genefu 2.40.0 released package functions/model',
 optional_ROR='SKIP: no loaded clinical-assay implementation/calibration'), 'work/PAM50_QC.json')
saveRDS(list(model=model,index=i,details=details),file.path(out,'work/PAM50_INPUT_MAPPING.rds'))
cat('PAM50 coverage',sum(details$available),'/50; missing',paste(details$centroid_probe[!details$available],collapse=', '),'\n')
