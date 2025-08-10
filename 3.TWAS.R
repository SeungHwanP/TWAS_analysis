#############
# make sumstats
#############
setwd("~/result_pleio")
rm(list=ls())
library(data.table)
library(dplyr)
library(glue)

folders <- list.dirs("~/2.PLEIO",recursive = F)
folders <- grep("/2.",folders,value=T)

#SNPid-info
snp_list <- fread("~/2.PLEIO/SNP_list.txt")
snp_list1 <- snp_list2 <- snp_list;rm(snp_list)
snp_list1 <- snp_list1 %>% mutate(index = paste0(CHR,":", BP, "_", ALLELE1, "_", ALLELE0))
snp_list2 <- snp_list2 %>% mutate(index = paste0(CHR,":", BP, "_", ALLELE0, "_", ALLELE1))
#bim
bim <- Reduce(rbind,lapply(1:22,function(x) fread(glue("~/LDREF/chr{x}.bim"))))
bim <- bim %>% mutate(index = paste0(V1,":", V4, "_", V5, "_", V6))
#merge between snp_list and bim
snp_list1 <- snp_list1 %>% filter(index %in% bim$index)
snp_list2 <- snp_list2 %>% filter(index %in% bim$index)

for(folder in folders){
  print(folder)
  setwd(folder)
  txts <- list.files(".","txt")
  txts <- txts[!grepl("SNP_list",txts)]
  for(txt in txts){
    print(txt)
    pleio <- fread(glue("zcat {txt}"))
    
    #data merge
    snp_list1_tmp <- merge(snp_list1, pleio, by='SNP',sort=F)
    snp_list2_tmp <- merge(snp_list2, pleio, by='SNP',sort=F)
    snp_list_tmp <- rbind(snp_list1_tmp,snp_list2_tmp)
    snp_list_tmp <- na.omit(snp_list_tmp)
    res <- merge(snp_list_tmp,bim[,c("index","V2","V1","V4")],by.x="SNP",by.y="V2",sort=F)[,-1]
    res <- res[order(res$CHR,res$BP),]
    name <- gsub(".txt.gz","",txt)
    system(glue("mkdir -p ~/3.TWAS/{folder}"))
    fwrite(res,glue("~/3.TWAS/{folder}/0.full.txt"),row.names=F,quote=F,sep="\t",col.names=T)
    # check each pheno
    name1 <- paste0("~/2.PLEIO/",gsub("1.","0.",name),".txt")
    phe1 <- fread(name1)
    phe1
    sum(duplicated(phe1$SNP))
    name2 <- paste0("~/2.PLEIO/",gsub("~/3.TWAS/","",folder),".txt")
    phe2 <- fread(name2)
    sum(duplicated(phe2$SNP))
    
    
    # make '.sumstats.gz'
    snp_list_tmp <- snp_list_tmp[order(snp_list_tmp$CHR,snp_list_tmp$BP),]
    sumstat <- data.table(SNP=snp_list_tmp$SNP,
                          A1=snp_list_tmp$ALLELE1, A2=snp_list_tmp$ALLELE0,
                          Z=snp_list_tmp$pleio_stat, N=459119)
    fwrite(sumstat,glue("~/3.TWAS/{folder}/1.sumstat.sumstats.gz"),row.names=F,quote=F,sep="\t",col.names=T)
    rm(sumstat,res)
  }
}


#############
# Run TWAS
#############
setwd("~/result_pleio")
rm(list=ls())
library(glue)
library(data.table)
folders <- list.files("~/GTEx/",".tar.gz")
folders <- gsub("~/GTEx/","",folders)
for(folder in folders){
  print(folder)
  system(glue("mkdir -p ~/3.TWAS/{folder}"))
  chr <- 1
  i <- 1
  numbers <- NULL
  for(chr in 1:22){
    if(length(system(glue("cat ~/3.TWAS/{folder}/Fusion_{chr}.log |grep complete"),intern=T)) <1){
      Rscript <- glue("/usr/local/bin/Rscript ~/fusion_twas/FUSION.assoc_test.R --sumstats ~/3.TWAS/{folder}/1.sumstat.sumstats.gz --weights ~/GTEx/{folder}.pos --weights_dir ~/GTEx/ --ref_ld_chr ~/LDREF/chr --chr {chr} --out ~/3.TWAS/{folder}/Fusion_{chr}.dat")
      if(is.null(numbers)){
        sbatch <- glue("/usr/bin/sbatch -n 1 --mem=30G -o ~/3.TWAS/{folder}/Fusion_{chr}.log -J {j}_{i}_{chr} --wrap='{Rscript}'")
      }else{
        sbatch <- glue("/usr/bin/sbatch -n 1 -d afterok:{numbers} --mem=30G -o ~/3.TWAS/{folder}/Fusion_{chr}.log -J {chr} --wrap='{Rscript}'")
      }
      numbers <- system(sbatch,intern=T)
      numbers <- as.numeric(gsub("Submitted batch job ","",numbers))
    }
  }
  i <- i+1
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


