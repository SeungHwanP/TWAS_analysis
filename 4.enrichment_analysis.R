rm(list=ls())
library(enrichR)
library(data.table)
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}
dbs <- listEnrichrDbs()
sort(dbs$libraryName)

ress <- NULL
dbs <- c("KEGG_2021_Human","GWAS_Catalog_2023","DisGeNET","Reactome_Pathways_2024",
         "Human_Phenotype_Ontology","WikiPathways_2024_Human",
         "GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021")

for(phe1 in c("COVID_Sev(HGI)","COVID_Inf(HGI)")){
  for(phe2 in c("T2D","BMI","WC","WHR")){
    genes <- fread(glue("~/3.TWAS/3-2.result_geneinfo_{phe1}_{phe2}.csv"))
    genes <- genes[genes$`Gene type` != "",]
    table(genes$`Gene type`)
    genes <- rbind(genes[grep("IncRNA",genes$`Gene type`),],
                   genes[grep("protein_coding",genes$`Gene type`),])
    genes <- unique(genes$symbol)
    enriched <- enrichr(genes, dbs)
    
    make_enriched <- lapply(1:length(dbs),function(x) enriched[[x]][,c("Term","Overlap","P.value","Adjusted.P.value","Genes")])
    
    res <- lapply(1:length(dbs),function(i){
      res <- make_enriched[[i]]
      tmp <- strsplit(res$Overlap,"/",fixed=T)
      res$N_overlap <- as.numeric(lapply(1:length(tmp),function(x) tmp[[x]][1]))
      res$N_genes <- as.numeric(lapply(1:length(tmp),function(x) tmp[[x]][2]))
      
      res <- res[,c("Term","N_genes","N_overlap","P.value","Adjusted.P.value","Genes")]
      res <- res[res$Adjusted.P.value < 0.05,]
      colnames(res) <- c("GeneSet","N_genes","N_overlap","p","Benjamini-Hochberg False Discovery Rate","genes")
      if(nrow(res)!=0){
        data.table(cbind(Category=dbs[i],res))
      }else{
        NULL
      }
    })
    
    lapply(1:9,function(x) dim(res[[x]]))
    null_data <- data.frame(matrix(NA,nrow=1,ncol=ncol(res[[1]])))
    colnames(null_data) <- colnames(res[[1]])
    ress <- rbind(res[[1]],
                  res[[2]],
                  res[[3]],
                  res[[4]],
                  res[[5]],
                  res[[6]],
                  res[[7]],
                  res[[8]],
                  res[[9]])
    
    ress$Category <- gsub("GO_Biological_Process_2021","GO Biological Process",ress$Category)
    ress$Category <- gsub("GO_Cellular_Component_2021","GO Cellular Component",ress$Category)
    ress$Category <- gsub("GO_Molecular_Function_2021","GO Molecular Function",ress$Category)
    ress$Category <- gsub("GWAS_Catalog_2023","GWAS Catalog",ress$Category)
    ress$Category <- gsub("KEGG_2021_Human","KEGG",ress$Category)
    ress$Category <- gsub("Reactome_Pathways_2024","Reactome Pathways",ress$Category)
    
    fwrite(ress,glue("~/4.Enrichment/enrichr_{phe1}_{phe2}_cut.txt"),row.names=F,quote=F,sep="\t",col.names=T)
    fwrite(ress,glue("~/4.Enrichment/enrichr_{phe1}_{phe2}_cut.csv"),row.names=F,quote=F,sep=",",col.names=T)
    
  }
}
