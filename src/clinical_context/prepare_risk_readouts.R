source(file.path(Sys.getenv('BCPN_REPOSITORY'), 'src/clinical_context/runtime.R'))
stopifnot(!file.exists(file.path(out,'work/PREPARED.rds')))
input <- readRDS(file.path(root,'enrichment/work/INPUTS.rds'))
previous <- readRDS(file.path(root,'molecular/work/ANALYSIS_RESULTS.rds'))
pam <- read_csv(file.path(root,'molecular/molecular_PAM50_SAMPLE_RESULTS.csv'))
pm <- previous$pair_map
if(is.null(pm))pm <- readRDS(file.path(root,'enrichment/work/TME_RESULTS.rds'))$pair_map
stopifnot(nrow(input$x)==12502,ncol(input$x)==72,nrow(pm)==36,
 identical(colnames(input$x),input$samples$GSM),sum(read_csv(file.path(root,'molecular/molecular_PAM50_DISCORDANCE.csv'))$PAM50_discordant)==8)
ann <- input$mapping
cat('Mapping columns:',paste(names(ann),collapse=', '),'\n')
usable <- ann$annotation_usable & !is.na(ann$entrez_id)
dups <- duplicated(ann$entrez_id[usable]) | duplicated(ann$entrez_id[usable],fromLast=TRUE)
keep <- which(usable)[!dups]
primary <- match(pm$Primary,colnames(input$x))
node <- match(pm$Node,colnames(input$x))
stopifnot(!anyNA(primary),!anyNA(node),length(unique(c(primary,node)))==72)
tables <- list(GGI=read_csv(system.file('extdata/sotiriou2006_sig_ggi.csv',package='genefu')),
 GENE70=read_csv(system.file('extdata/veer2002_sig_gene70.csv',package='genefu')),
 RS=read_csv(system.file('extdata/paik2004_sig_oncotypedx.csv',package='genefu')))
details <- list(); definitions <- list()
for(name in names(tables)){
 tab <- tables[[name]]
 if(name=='RS')tab <- tab[tab$group!='reference',,drop=FALSE]
 id <- as.character(tab$EntrezGene.ID)
 if(name=='GGI')coefficient <- ifelse(tab$grade==1,-1,1)
 if(name=='GENE70')coefficient <- -sign(tab$average.good.prognosis.profile)
 if(name=='RS')coefficient <- tab$weight
 idx <- keep[match(id,as.character(ann$entrez_id[keep]))]
 sd0 <- vapply(idx,function(i)if(is.na(i))NA_real_ else sd(input$x[i,primary]),0.0)
 available <- !is.na(idx)&is.finite(sd0)&sd0>0&is.finite(coefficient)
 conflicting <- vapply(id,function(v)if(is.na(v))FALSE else length(unique(coefficient[!is.na(id)&id==v]))>1,TRUE)
 available <- available&!conflicting
 symbol <- if('HUGO.gene.symbol'%in%names(tab))tab$HUGO.gene.symbol else tab$symbol
 probe <- if('probe'%in%names(tab))tab$probe else tab$probe.affy
 aa <- data.frame(signature=name,source_probe=probe,source_symbol=symbol,EntrezGene.ID=id,
 direction_coefficient=coefficient,selected_gene=ann$gene_symbol[idx],approved_symbol=ann$approved_symbol[idx],
 selected_probe=input$full$probe_id[idx],available=available,conflicting_direction=conflicting,
 Primary_gene_mean=vapply(idx,function(i)if(is.na(i))NA_real_ else mean(input$x[i,primary]),0.0),Primary_gene_SD=sd0)
 details[[name]] <- aa
 use <- which(available&!duplicated(id))
 definitions[[name]] <- list(index=idx[use],ids=id[use],symbols=ann$approved_symbol[idx[use]],
 coefficients=coefficient[use],groups=if(name=='RS')tab$group[use] else rep(name,length(use)),
 full_source_rows=nrow(tab),available_source_rows=sum(available),unique_source_genes=length(unique(id[!is.na(id)])),
 used_unique_genes=length(use),mean=aa$Primary_gene_mean[use],sd=sd0[use])
 cat(name,':',sum(available),'/',nrow(tab),'source features;',length(use),'unique used genes\n')
}
write_csv(do.call(rbind,details),'work/RISK_GENE_MAPPING_QC.csv')
saveRDS(list(input=input,primary=primary,node=node,pm=pm,pam=pam,definitions=definitions,
 nodal_sets=previous$sets[c('ECM','BCR','E2F','G2M','Adipogenesis')]),file.path(out,'work/PREPARED.rds'))
coverage <- do.call(rbind,lapply(names(definitions),function(n){d<-definitions[[n]];data.frame(signature=n,
 source_features=d$full_source_rows,mapped_features=d$available_source_rows,
 unique_source_genes=d$unique_source_genes,used_unique_genes=d$used_unique_genes,
 source_feature_coverage_pct=100*d$available_source_rows/d$full_source_rows,
 status=if(n=='GENE70'&&d$available_source_rows/d$full_source_rows<.8)'SKIP_BELOW_80_PERCENT_OF_ORIGINAL_70_FEATURES' else 'READY_RESEARCH_APPROXIMATION')}))
write_csv(coverage,'work/RISK_SIGNATURE_COVERAGE.csv')
write_csv(data.frame(score=names(previous$sets[c('ECM','BCR','E2F','G2M','Adipogenesis')]),
 full_members=vapply(previous$sets[c('ECM','BCR','E2F','G2M','Adipogenesis')],paste,'',collapse=';')),'work/NODAL_SIGNATURE_DEFINITIONS.csv')
capture.output(sessionInfo(),file=file.path(out,'work/PRE_ANALYSIS_R_SESSION.txt'))
print(coverage)
