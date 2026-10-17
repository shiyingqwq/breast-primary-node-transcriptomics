source(file.path(Sys.getenv('BCPN_REPOSITORY'), 'src/enrichment/runtime.R'))
library(fgsea)
input <- readRDS(file.path(out,'work/INPUTS.rds'))
full <- input$full;mapping <- input$mapping;cand <- input$candidates
parse_gmt <- function(path)fgsea::gmtPathways(path)
go <- parse_gmt(file.path(src,'GO_BP.v2025.1.Hs.gmt'))
reactome <- parse_gmt(file.path(src,'REACTOME.v2025.1.Hs.gmt'))
reliable <- mapping[mapping$annotation_usable,]
to_original <- function(paths)lapply(paths,function(genes)unique(reliable$gene_symbol[reliable$approved_symbol%in%genes]))
go <- to_original(go);reactome <- to_original(reactome)
kegg_link <- read.delim(file.path(src,'KEGG_GENE_PATHWAY.tsv'),header=FALSE,col.names=c('gene','pathway'))
kegg_names <- read.delim(file.path(src,'KEGG_PATHWAYS.tsv'),header=FALSE,col.names=c('pathway','description'))
kegg_names$pathway <- sub('^path:','',kegg_names$pathway)
kegg_link$pathway <- sub('^path:','',kegg_link$pathway)
kegg_entrez <- split(sub('^hsa:','',kegg_link$gene),kegg_link$pathway)
kegg <- lapply(kegg_entrez,function(ids)unique(reliable$gene_symbol[!is.na(reliable$entrez_id)&reliable$entrez_id%in%ids]))
collections <- list(GO_BP=go,KEGG=kegg,REACTOME=reactome)
description <- function(db,id) {
  if(db=='KEGG')return(kegg_names$description[match(id,kegg_names$pathway)])
  gsub('_',' ',sub('^(GOBP_|REACTOME_)','',id))
}
ora <- function(genes,db,set_name,minSize=10,maxSize=500) {
  paths <- collections[[db]]
  background <- intersect(full$gene_symbol,unique(unlist(paths)))
  genes <- intersect(unique(genes),background)
  sizes <- lengths(paths);paths <- paths[sizes>=minSize&sizes<=maxSize]
  rows <- lapply(names(paths),function(id) {
    members <- intersect(paths[[id]],background);hit <- intersect(genes,members)
    p <- if(length(hit))phyper(length(hit)-1,length(members),length(background)-length(members),
                              length(genes),lower.tail=FALSE) else 1
    data.frame(database=db,set=set_name,pathway=id,description=description(db,id),
      raw_P=p,gene_count=length(hit),gene_ratio=paste0(length(hit),'/',length(genes)),
      gene_ratio_numeric=length(hit)/length(genes),background_ratio=paste0(length(members),'/',length(background)),
      contributing_genes=paste(hit,collapse=';'),input_tested_background_N=12502,
      effective_annotated_background_N=length(background),effective_query_N=length(genes),
      pathway_tested_size=length(members))
  })
  result <- do.call(rbind,rows)
  result$BH_q <- p.adjust(result$raw_P,'BH')
  result$significant <- result$BH_q<.05
  result$BH_family <- paste(set_name,db,sep=' x ')
  result[order(result$BH_q,result$raw_P,result$pathway),]
}
sets <- list(ALL=cand$gene_symbol,NODE_UP=cand$gene_symbol[cand$bulk_direction=='Node_up'],
             NODE_DOWN=cand$gene_symbol[cand$bulk_direction=='Node_down'])
for(n in names(sets)) {
  results <- do.call(rbind,lapply(names(collections),function(db)ora(sets[[n]],db,n)))
  write_csv(results,paste0('enrichment_ORA_',n,'.csv'))
  cat('ORA',n,'significant',sum(results$significant),'\n');flush.console()
}
coverage <- do.call(rbind,lapply(names(collections),function(db){
  annotated <- unique(unlist(collections[[db]]))
  data.frame(database=db,pathways=length(collections[[db]]),input_background_N=12502,
    annotated_background_N=length(intersect(annotated,full$gene_symbol)),
    mapped_candidates=length(intersect(annotated,cand$gene_symbol)),unannotated_candidates=
      paste(setdiff(cand$gene_symbol,annotated),collapse=';'))
}))
write_csv(coverage,'work/DATABASE_GENE_COVERAGE.csv')
saveRDS(collections,file.path(out,'work/COLLECTIONS.rds'))
saveRDS(ora,file.path(out,'work/ORA_FUNCTION.rds'))
file.copy(file.path(root,'bulk_HALLMARK_ALL_TESTED.csv'),file.path(out,'enrichment_GSEA_HALLMARK.csv'),overwrite=TRUE)
rank <- setNames(full$t,full$gene_symbol);rank <- rank[order(-rank,names(rank),method='radix')]
stopifnot(length(rank)==12502,all(is.finite(rank)),!anyDuplicated(names(rank)))
write_csv(data.frame(gene=names(rank),moderated_t=unname(rank)),'work/GSEA_RANKING.csv')
for(db in names(collections)) {
  cat('GSEA START',db,'\n');flush.console()
  set.seed(20261006)
  g <- as.data.frame(fgsea::fgseaMultilevel(pathways=collections[[db]],stats=rank,
    minSize=15,maxSize=500,eps=1e-10,scoreType='std',gseaParam=1,sampleSize=101,nPermSimple=1000,nproc=1))
  g$padj <- p.adjust(g$pval,'BH',n=nrow(g))
  g$leadingEdge <- vapply(g$leadingEdge,paste,collapse=';',FUN.VALUE=character(1))
  g$description <- vapply(g$pathway,function(id)description(db,id),'')
  g$database <- db;g$ranking <- 'Original 12502 limma moderated t; positive Node-minus-Primary'
  g$method <- 'fgseaMultilevel; min15 max500 eps1e-10; seed20261006'
  g <- g[order(g$padj,g$pval,g$pathway),]
  write_csv(g,paste0('enrichment_GSEA_',db,'.csv'))
  cat('GSEA DONE',db,'tested',nrow(g),'significant',sum(g$padj<.05,na.rm=TRUE),'\n');flush.console()
}
save_json(list(ORA_method='R upper-tail hypergeometric; database-specific annotated tested background',
  ORA_minSize=10,ORA_maxSize=500,ORA_BH='separate set x database; includes zero-hit eligible terms',
  GSEA_rank_N=12502,GSEA_statistic='moderated t',GSEA_eps=1e-10,
  selected_Hallmark_reused=TRUE,Hallmark_rerun=FALSE), 'work/ENRICHMENT_COMPLETE.json')
capture.output(sessionInfo(),file=file.path(out,'work/ENRICHMENT_SESSION.txt'))
