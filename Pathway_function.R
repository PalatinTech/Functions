library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(EnhancedVolcano)

# This file will have functions to perform pathway enrichment analysis
#input is a dataframe with log2FC values and gene names
KEGGpath<-function(input)
{
  hs <- org.Hs.eg.db
  allgenes<-toupper(input$Gene.name)
  
  allgenes_out<-AnnotationDbi::select(hs, keys=allgenes, columns = c("ENTREZID","SYMBOL"),keytype = "SYMBOL")
  head(allgenes_out)
  
  tmp<-data.frame(log2FC=input$log2FoldChange,SYMBOL=toupper(input$Gene.name))
  a<-merge(tmp,allgenes_out, by=c('SYMBOL'))
  a<-a[!is.na(a$ENTREZID),]
  genelist<-a$log2FC
  names(genelist)<-a$ENTREZID
  genelist
  kk2 <- gseKEGG(geneList     = sort(genelist,decreasing = TRUE),
                 organism     = 'hsa',
                 minGSSize    = 10,
                 pvalueCutoff = 0.2,
                 verbose      = FALSE)
  return(kk2)
}

#GSEA Molecular SIgnature database
GSEApath<-function(input)
{
  #running MsigDB
  #getting all hallmark gene sets
  msigdbr_t2g = msigdbr(species = "Homo sapiens",category = "H")%>%
    dplyr::select(gs_name,gene_symbol)
  head(msigdbr_t2g)
  
  #Using GSEA for gene set analysis
  #prepare input
  glist<-input$log2FoldChange
  names(glist)<-toupper((input$Gene.name))
  glist = sort(glist, decreasing = TRUE)
  em2<-GSEA(glist,TERM2GENE = msigdbr_t2g,pAdjustMethod = "none",pvalueCutoff = "0.2")
  clusterProfiler::dotplot(em2)
  return(em2)
}


#Gene Set Enrichment Analysis of Gene Ontology
GSEgopath<-fnciton(input)
{
  glist<-input$log2FoldChange
  names(glist)<-toupper((input$Gene.name))
  glist = sort(glist, decreasing = TRUE)
  gse<-gseGO(geneList = glist,
             ont="BP",# biological pathway
             keyType = "SYMBOL",
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             by="fgsea",
             pAdjustMethod = "BH")
  clusterProfiler::dotplot(gse)
  return(gse)
}


