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
library(dplyr)
library(glue)

## COVID_Severe(HGI)
HGI <- fread(glue("~/1.GWAS/COVID19_HGI_B2_ALL_eur_leave23andme_20220403_GRCh37.tsv.gz")) 
HGI$N <- HGI$all_inv_var_meta_cases + HGI$all_inv_var_meta_controls
HGI$Z <- HGI$all_inv_var_meta_beta/HGI$all_inv_var_meta_sebeta

pleio = data.frame(SNP=HGI$rsid,
                   A1=HGI$ALT, 
                   A2=HGI$REF, 
                   Z=HGI$Z,
                   P=HGI$all_inv_var_meta_p,
                   N=HGI$N)


ukb_df <- fread(glue("~/2.PLEIO/0.UKB_T2D.txt"))

df_merge <- merge(ukb_df, pleio, by.x = c("SNP", "A1", "A2"), by.y=c("SNP", "A1", "A2"))
df_merge2 <- merge(ukb_df, pleio, by.x = c("SNP", "A1", "A2"), by.y=c("SNP", "A2", "A1"))
df_merge2$Z.y <- -df_merge2$Z.y

df_last <- bind_rows(df_merge, df_merge2)
df_last <- df_last %>% select(SNP, A1, A2, Z.y, P.y, N.y) %>% rename(Z=Z.y, P=P.y, N=N.y)
dim(df_last)
head(df_last)

# pleio input 파일 형식 저장 
write.table(df_last, glue("~/2.PLEIO/0.HGI_COVID_Sev.txt"), row.names=F,quote=F,sep="\t")

## COVID_Infection(HGI)
HGI <- fread(glue("~/1.GWAS/COVID19_HGI_C2_ALL_eur_leave23andme_20220403_GRCh37.tsv.gz")) 
HGI$N <- HGI$all_inv_var_meta_cases + HGI$all_inv_var_meta_controls
HGI$Z <- HGI$all_inv_var_meta_beta/HGI$all_inv_var_meta_sebeta

pleio = data.frame(SNP=HGI$rsid,
                   A1=HGI$ALT, 
                   A2=HGI$REF, 
                   Z=HGI$Z,
                   P=HGI$all_inv_var_meta_p,
                   N=HGI$N)


ukb_df <- fread(glue("~/2.PLEIO/0.UKB_T2D.txt"))

df_merge <- merge(ukb_df, pleio, by.x = c("SNP", "A1", "A2"), by.y=c("SNP", "A1", "A2"))
df_merge2 <- merge(ukb_df, pleio, by.x = c("SNP", "A1", "A2"), by.y=c("SNP", "A2", "A1"))
df_merge2$Z.y <- -df_merge2$Z.y

df_last <- bind_rows(df_merge, df_merge2)
df_last <- df_last %>% select(SNP, A1, A2, Z.y, P.y, N.y) %>% rename(Z=Z.y, P=P.y, N=N.y)
dim(df_last)
head(df_last)

# pleio input 파일 형식 저장 
write.table(df_last, glue("~/2.PLEIO/0.HGI_COVID_Inf.txt"), row.names=F,quote=F,sep="\t")


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

## T2D_adjBMI
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


## T2D_BMI
T2D_BMI_UKB_COVID_Sev <- rbind(UKB_T2D_BMI, UKB_COVID_Sev)
colnames(T2D_BMI_UKB_COVID_Sev) <- c("FILE", "TYPE", "SPREV", "PPREV", "NAME")
write.table(T2D_BMI_UKB_COVID_Sev,
            "~/2.PLEIO/UKB_T2D_BMI/UKB_T2D_BMI_UKB_COVID_Sev_input_list.txt", row.names=F,quote=F,sep="\t")

T2D_BMI_HGI_COVID_Sev <- rbind(UKB_T2D_BMI, HGI_COVID_Sev)
colnames(T2D_BMI_HGI_COVID_Sev) <- c("FILE", "TYPE", "SPREV", "PPREV", "NAME")
T2D_BMI_HGI_COVID_Sev
write.table(T2D_BMI_HGI_COVID_Sev,
            "~/2.PLEIO/UKB_T2D_BMI/UKB_T2D_BMI_HGI_COVID_Sev_input_list.txt", row.names=F,quote=F,sep="\t")

T2D_BMI_UKB_COVID_Inf <- rbind(UKB_T2D_BMI, UKB_COVID_Inf)
colnames(T2D_BMI_UKB_COVID_Inf) <- c("FILE", "TYPE", "SPREV", "PPREV", "NAME")
write.table(T2D_BMI_UKB_COVID_Inf,
            "~/2.PLEIO/UKB_T2D_BMI/UKB_T2D_BMI_UKB_COVID_Inf_input_list.txt", row.names=F,quote=F,sep="\t")

T2D_BMI_HGI_COVID_Inf <- rbind(UKB_T2D_BMI, HGI_COVID_Inf)
colnames(T2D_BMI_HGI_COVID_Inf) <- c("FILE", "TYPE", "SPREV", "PPREV", "NAME")
T2D_BMI_HGI_COVID_Inf
write.table(T2D_BMI_HGI_COVID_Inf,
            "~/2.PLEIO/UKB_T2D_BMI/UKB_T2D_BMI_HGI_COVID_Inf_input_list.txt", row.names=F,quote=F,sep="\t")


## T2D_WHR
T2D_WHR_UKB_COVID_Sev <- rbind(UKB_T2D_WHR, UKB_COVID_Sev)
colnames(T2D_WHR_UKB_COVID_Sev) <- c("FILE", "TYPE", "SPREV", "PPREV", "NAME")
write.table(T2D_WHR_UKB_COVID_Sev,
            "~/2.PLEIO/UKB_T2D_WHR/UKB_T2D_WHR_UKB_COVID_Sev_input_list.txt", row.names=F,quote=F,sep="\t")

T2D_WHR_HGI_COVID_Sev <- rbind(UKB_T2D_WHR, HGI_COVID_Sev)
colnames(T2D_WHR_HGI_COVID_Sev) <- c("FILE", "TYPE", "SPREV", "PPREV", "NAME")
T2D_WHR_HGI_COVID_Sev
write.table(T2D_WHR_HGI_COVID_Sev,
            "~/2.PLEIO/UKB_T2D_WHR/UKB_T2D_WHR_HGI_COVID_Sev_input_list.txt", row.names=F,quote=F,sep="\t")

T2D_WHR_UKB_COVID_Inf <- rbind(UKB_T2D_WHR, UKB_COVID_Inf)
colnames(T2D_WHR_UKB_COVID_Inf) <- c("FILE", "TYPE", "SPREV", "PPREV", "NAME")
write.table(T2D_WHR_UKB_COVID_Inf,
            "~/2.PLEIO/UKB_T2D_WHR/UKB_T2D_WHR_UKB_COVID_Inf_input_list.txt", row.names=F,quote=F,sep="\t")

T2D_WHR_HGI_COVID_Inf <- rbind(UKB_T2D_WHR, HGI_COVID_Inf)
colnames(T2D_WHR_HGI_COVID_Inf) <- c("FILE", "TYPE", "SPREV", "PPREV", "NAME")
T2D_WHR_HGI_COVID_Inf
write.table(T2D_WHR_HGI_COVID_Inf,
            "~/2.PLEIO/UKB_T2D_WHR/UKB_T2D_WHR_HGI_COVID_Inf_input_list.txt", row.names=F,quote=F,sep="\t")


## T2D_WC
T2D_WC_UKB_COVID_Sev <- rbind(UKB_T2D_WC, UKB_COVID_Sev)
colnames(T2D_WC_UKB_COVID_Sev) <- c("FILE", "TYPE", "SPREV", "PPREV", "NAME")
write.table(T2D_WC_UKB_COVID_Sev,
            "~/2.PLEIO/UKB_T2D_WC/UKB_T2D_WC_UKB_COVID_Sev_input_list.txt", row.names=F,quote=F,sep="\t")

T2D_WC_HGI_COVID_Sev <- rbind(UKB_T2D_WC, HGI_COVID_Sev)
colnames(T2D_WC_HGI_COVID_Sev) <- c("FILE", "TYPE", "SPREV", "PPREV", "NAME")
T2D_WC_HGI_COVID_Sev
write.table(T2D_WC_HGI_COVID_Sev,
            "~/2.PLEIO/UKB_T2D_WC/UKB_T2D_WC_HGI_COVID_Sev_input_list.txt", row.names=F,quote=F,sep="\t")

T2D_WC_UKB_COVID_Inf <- rbind(UKB_T2D_WC, UKB_COVID_Inf)
colnames(T2D_WC_UKB_COVID_Inf) <- c("FILE", "TYPE", "SPREV", "PPREV", "NAME")
write.table(T2D_WC_UKB_COVID_Inf,
            "~/2.PLEIO/UKB_T2D_WC/UKB_T2D_WC_UKB_COVID_Inf_input_list.txt", row.names=F,quote=F,sep="\t")

T2D_WC_HGI_COVID_Inf <- rbind(UKB_T2D_WC, HGI_COVID_Inf)
colnames(T2D_WC_HGI_COVID_Inf) <- c("FILE", "TYPE", "SPREV", "PPREV", "NAME")
T2D_WC_HGI_COVID_Inf
write.table(T2D_WC_HGI_COVID_Inf,
            "~/2.PLEIO/UKB_T2D_WC/UKB_T2D_WC_HGI_COVID_Inf_input_list.txt", row.names=F,quote=F,sep="\t")


