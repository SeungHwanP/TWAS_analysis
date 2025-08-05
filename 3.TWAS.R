#############
# make sumstats
#############



#############
# Run TWAS
#############
setwd("~/result_pleio")
rm(list=ls())
library(glue)
library(data.table)
folders <- list.files("~/GTEx/",".tar.gz")
for(dir in dirs){
  print(dir)
  folder <- folders[1]
  i <- 1
  for(folder in folders){
    print(folder)
    system(glue("mkdir -p {dir}/{folder}"))
    chr <- 1
    numbers <- NULL
    for(chr in 1:22){
      if(length(system(glue("cat {dir}/{folder}/Fusion_{chr}.log |grep complete"),intern=T)) <1){
        Rscript <- glue("/usr/local/bin/Rscript ~/fusion_twas/FUSION.assoc_test.R --sumstats {dir}/1.sumstat.sumstats.gz --weights ~/GTEx/{folder}.pos --weights_dir ~/GTEx/ --ref_ld_chr ~/LDREF/chr --chr {chr} --out {dir}/{folder}/Fusion_{chr}.dat")
        if(is.null(numbers)){
          sbatch <- glue("/usr/bin/sbatch -n 1 --mem=30G -o {dir}/{folder}/Fusion_{chr}.log -J {j}_{i}_{chr} --wrap='{Rscript}'")
        }else{
          sbatch <- glue("/usr/bin/sbatch -n 1 -d afterok:{numbers} --mem=30G -o {dir}/{folder}/Fusion_{chr}.log -J {j}_{i}_{chr} --wrap='{Rscript}'")
        }
        numbers <- system(sbatch,intern=T)
        numbers <- as.numeric(gsub("Submitted batch job ","",numbers))
      }
    }
    i <- i+1
  }
  j <- j+1
}

#############
# merge data
#############
rm(list=ls())
library(data.table)
library(glue)
dirs <- list.dirs(".",recursive = F)
gene_list <- fread("~/fusion_twas/glist-hg19")
ress <- NULL
gene_name <- fread("~/fusion_twas/hgnc_complete_set.txt")[,c("ensembl_gene_id","symbol")]
folders <- list.files("~/GTEx/",".tar.gz")
folders <- gsub(".tar.gz","",folders)

dir <- dirs[1]
for(dir in dirs){
  tmp_dir <- strsplit(dir,"/",fixed=T)[[1]][2]
  print(glue("======={tmp_dir}======="))
  folders <- list.dirs(glue("{tmp_dir}"),recursive=F)[-2]
  folder <- folders[1]
  for(folder in folders){
    #print(folder)
    folders2 <- list.dirs(folder,recursive=F)
    ress <- NULL
    folder2 <- folders2[1]
    i <- 1
    for(folder2 in folders2){
      print(folder2)
      res <- Reduce(rbind,lapply(list.files(glue("{folder2}"),"dat",full.names = T),function(x) fread(x)))[,-c("PANEL","FILE")]
      res <- na.omit(res)
      N_tissue <- length(unique(res$ID))
      tissue <- strsplit(folder2,"GTExv8.EUR.")[[1]][2]
      res$tissue <- tissue
      res$N_tissue <- N_tissue
      tmp <- strsplit(res$ID,".",fixed=T)
      res$ID <- unlist(lapply(1:length(tmp),function(x) tmp[[x]][1]))
      res <- merge(res,gene_name,by.x="ID",by.y="ensembl_gene_id",all.x=T,sort=F)
      res <- na.omit(res)
      fwrite(res,glue("{folder}/3.{i}_{tissue}.csv"),quote=F,sep=",",col.names=T,row.names=F)
      ress <- rbind(ress,res)
      i <- i+1
    }
    ress <- ress[order(ress$symbol),]
    fwrite(ress,glue("{folder}/3.47_result.csv"),quote=F,sep=",",col.names=T,row.names=F)
  }
}


