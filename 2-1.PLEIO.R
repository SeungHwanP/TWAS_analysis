#############
# make PLEIO input data (GWAS)
#############
rm(list=ls())
library(data.table)
library(glue)

for(phe in c("UKB_COVID_Inf","UKB_COVID_Sev","UKB_T2D_adjBMI","UKB_BMI","UKB_WHR","UKB_WC"){
  bolt <- fread(glue("zcat ~/1.GWAS/BOLT_{phe}.bgen.gz"))
  pleio = data.frame(SNP=bolt$SNP,
                   A1=bolt$ALLELE1, 
                   A2=bolt$ALLELE0, 
                   Z=bolt$BETA/bolt$SE,
                   P=bolt$P_BOLT_LMM,
                   N=bolt$N)
  fwrite(pleio,glue("~/2.PLEIO/0.{phe}.txt"),row.names=F,quote=F,sep="\t",col.names=T)
}


#############
# make PLEIO input data (Nature data)
#############
rm(list=ls())
library(data.table)
library(glue)

for(phe in c("UKB_COVID_Inf","UKB_COVID_Sev","UKB_T2D_adjBMI","UKB_BMI","UKB_WHR","UKB_WC"){
  bolt <- fread(glue("zcat ~/1.GWAS/BOLT_{phe}.bgen.gz"))
  pleio = data.frame(SNP=bolt$SNP,
                   A1=bolt$ALLELE1, 
                   A2=bolt$ALLELE0, 
                   Z=bolt$BETA/bolt$SE,
                   P=bolt$P_BOLT_LMM,
                   N=bolt$N)
  fwrite(pleio,glue("~/2.PLEIO/0.{phe}.txt"),row.names=F,quote=F,sep="\t",col.names=T)
}




#############
# make PLEIO input data (data list)
#############
UKB_T2D_adjBMI <- c(glue("~/2.PLEIO/0.UKB_T2D_adjBMI.txt"), "binary", 0.0805, 0.0805, "T2D_adjBMI")
UKB_BMI <- c(glue("~/2.PLEIO/0.UKB_BMI.txt"), "quantitative", "nan", "nan", "BMI")
UKB_WHR <- c(glue("~/2.PLEIO/0.UKB_WHR.txt"), "quantitative", "nan", "nan", "WHR")
UKB_WC <- c(glue("~/2.PLEIO/0.UKB_WC.txt"), "quantitative", "nan", "nan", "WC")
UKB_COVID_Inf <- c(glue("~/2.PLEIO/0.UKB_COVID_Inf.txt"), "binary", 0.220577, 0.220577, "UKB_COVID_Inf")
UKB_COVID_Sev <- c(glue("~/2.PLEIO/0.UKB_COVID_Sev.txt"), "binary", 0.0162877, 0.0162877, "UKB_COVID_Sev")
HGI_COVID_Inf <- c(glue("~/2.PLEIO/0.HGI_COVID_Inf.txt"), "binary", 0.047198921, 0.047198921, "HGI_COVID_Inf")
HGI_COVID_Sev <- c(glue("~/2.PLEIO/0.HGI_COVID_Sev.txt"), "binary", 0.015519796, 0.015519796, "HGI_COVID_Sev")


T2D_adjBMI_UKB_COVID_Sev <- rbind(UKB_T2D_adjBMI, UKB_COVID_Sev)
colnames(T2D_adjBMI_UKB_COVID_Sev) <- c("FILE", "TYPE", "SPREV", "PPREV", "NAME")
write.table(T2D_adjBMI_UKB_COVID_Sev,
            "~/2.PLEIO/UKB_T2D_adjBMI/UKB_T2D_adjBMI_UKB_COVID_Sev_input_list.txt", row.names=F,quote=F,sep="\t")

T2D_adjBMI_HGI_COVID_Sev <- rbind(UKB_T2D_adjBMI, HGI_COVID_Sev)
colnames(T2D_adjBMI_HGI_COVID_Sev) <- c("FILE", "TYPE", "SPREV", "PPREV", "NAME")
T2D_adjBMI_HGI_COVID_Sev
write.table(T2D_adjBMI_HGI_COVID_Sev,
            "~/2.PLEIO/UKB_T2D_adjBMI/UKB_T2D_adjBMI_HGI_COVID_Sev_input_list.txt", row.names=F,quote=F,sep="\t")

T2D_adjBMI_UKB_COVID_Inf <- rbind(UKB_T2D_adjBMI, UKB_COVID_Inf)
colnames(T2D_adjBMI_UKB_COVID_Inf) <- c("FILE", "TYPE", "SPREV", "PPREV", "NAME")
write.table(T2D_adjBMI_UKB_COVID_Inf,
            "~/2.PLEIO/UKB_T2D_adjBMI/UKB_T2D_adjBMI_UKB_COVID_Inf_input_list.txt", row.names=F,quote=F,sep="\t")

T2D_adjBMI_HGI_COVID_Inf <- rbind(UKB_T2D_adjBMI, HGI_COVID_Inf)
colnames(T2D_adjBMI_HGI_COVID_Inf) <- c("FILE", "TYPE", "SPREV", "PPREV", "NAME")
T2D_adjBMI_HGI_COVID_Inf
write.table(T2D_adjBMI_HGI_COVID_Inf,
            "~/2.PLEIO/UKB_T2D_adjBMI/UKB_T2D_adjBMI_HGI_COVID_Inf_input_list.txt", row.names=F,quote=F,sep="\t")

