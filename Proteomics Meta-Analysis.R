Prot_dir=("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Proteome")

####Create function to format results for METAL software####
create_summary <- function(toptable = NULL,
                           dataset_label= NULL,
                           directory = getwd())
{
  PTN <- rownames(toptable)
  ALLELE1 <- rep(1,nrow(toptable))
  ALLELE2 <- rep(2,nrow(toptable))
  TESTSTAT <- toptable$t
  PVALUE <- toptable$P.Value
  EFFECTSIZE <- toptable$logFC
  SE <- toptable$SE
  results = data.frame(PTN,
                       ALLELE1,
                       ALLELE2,
                       TESTSTAT,
                       PVALUE,
                       EFFECTSIZE,
                       SE)
  write.table(results,
              file=paste0(directory,"/",dataset_label,".tbl"),
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE,
              sep="\t")
}

#### Gene SMART ####
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Proteome/GeneSMART")
library(tidyverse)

Proteomics=read.csv("MetaData_norm_and_batch.csv")
rownames(Proteomics)=Proteomics$Gene.names
Proteomics=Proteomics[,33:187]
colnames(Proteomics)
Proteomics=Proteomics[,-c(12,16,30,45,60,75,89,104,111,120,133,134,148)]#Remove technical replicates
Proteomics=Proteomics%>%
  dplyr::rename("SG156_PRE"="SG166_PRE", "SG166_PRE"="SG166_PRE.1")

pheno=read.csv("Pheno_proteomics.csv")
rownames(pheno)=pheno$X

colnames(Proteomics)==rownames(pheno)


#run limma
library(limma)

#change time for continuous
design=model.matrix(~VO2max_rel*sex+Timepoint+Age_PRE+intervention,data=pheno) #time as continuous

corfit_PRE_Prot <- duplicateCorrelation(Proteomics,design,block=pheno$ID)
cor=corfit_PRE_Prot$consensus.correlation
#Then, we need to input this inter-correlation into the linear model to take into account the paired design
fit = lmFit(Proteomics, design,block=pheno$ID,correlation=cor) #I have chosen the option "robust" which will downplay the influence of outliers. This takes significantly more time for your computer to run but is more robust and is always my preference.

#get info on eBayes in the vignette of the limma package
fit2 <- eBayes(fit)

coef = "VO2max_rel"
results <- topTable(fit2,
                    coef=coef,
                    number=Inf,
                    p.value=1)
SE <- fit2$sigma * fit2$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('GeneSMART_pvalhist_PTNs.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for VO2max Proteins",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "GeneSMART_VO2max_PTNs",
               directory = getwd())

coef = "Age_PRE"
results <- topTable(fit2,
                    coef=coef,
                    number=Inf,
                    p.value=1)

SE <- fit2$sigma * fit2$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('GeneSMART_pvalhist_PTN_Age.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Age Proteins",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "GeneSMART_Age_Proteins",
               directory = getwd())

coef = "sexM"
results <- topTable(fit2,
                    coef=coef,
                    number=Inf,
                    p.value=1)

SE <- fit2$sigma * fit2$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('GeneSMART_pvalhist_PTN_sex.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Sex PTNs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "GeneSMART_Sex",
               directory = getwd())

coef = "VO2max_rel:sexM"
results <- topTable(fit2,
                    coef=coef,
                    number=Inf,
                    p.value=1)

SE <- fit2$sigma * fit2$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('GeneSMART_pvalhist_Proteins_sex_vo2max.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Sex|VO2max DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "GeneSMART_Sex_VO2max_proteins",
               directory = getwd())


#### Deshmukh ####
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Proteome/Proteomics_Deshmukh")

library(readxl)
library(tidyverse)

pheno=read_excel("ACE meta for atul - 09.05.20.xlsx", skip=2)
pheno=pheno[1:16, 1:9]
colnames(pheno)=c("ID","Timepiont","Age","Height", "Weight", "HRmax", "VO2max","VO2max_rel", "aPPO (s)")

#Load proteomics data
Data_prot=read_xlsx("elife-69802-supp1-v1.xlsx")



####do imputations but skip pre-processing (already done by the author)####
library(DEP)
library(GenomeInfoDb)
library(SummarizedExperiment)  


#Make a data frame with only label, condition and replicate
pheno2=data.frame(label=colnames(Data_prot2[3:18]), condition=pheno$Timepiont, replicate=pheno$ID)

#Make_unique for meake_se
Data_prot[Data_prot=="Filtered"] <- NA
Data_prot2=as.data.frame(sapply(Data_prot[4:19], as.numeric))
Data_prot2=cbind("Genes"=Data_prot$Genes, "Accession"=Data_prot$ProteinAccessions, Data_prot2)
data_unique=make_unique(Data_prot2, "Genes", "Accession", delim = ";" )

#make data with only intensity cols
intensity_cols=c(3:18)

pheno2$replicate=as.numeric(pheno2$replicate)

pheno2$label=c("Pre_1","Pre_2","Pre_3", "Pre_4","Pre_5","Pre_6","Pre_7","Pre_8",
               "Post_1", "Post_2","Post_3","Post_4","Post_5","Post_6","Post_7","Post_8")

#Summarize experiment using make_se #this log transforms the data
data_se=DEP::make_se(data_unique, intensity_cols, pheno2)

data_filtered=filter_missval(data_se, thr=0.5)

data_filtered2=assay(data_filtered)

#QC Plots

plot_numbers(data_filtered)
plot_frequency(data_filtered)
plot_missval(data_filtered)

#Imputations and VSN normalisation
data_imp<-DEP:::impute(data_filtered,fun="man",shift=1.8,scale=0.3)

data_norm<-DEP::normalize_vsn(data_imp)

#Box plots before and after norm
plot_normalization(data_imp)
plot_normalization(data_norm)

#PCA before and After normalisation
plot_pca(data_imp, n=2480, indicate = "condition") +
  labs(title = "Before VSN Normalisation") # n is total number of proteins in the test data.

plot_pca(data_norm, n=2480, indicate = "condition", label = FALSE)+
  labs(title = "After VSN Normalisation")

Prot1=assay(data_norm)
Prot2=assay(data_imp)
colnames(Prot1)=c("1_PRE","2_PRE","3_PRE", "4_PRE", "5_PRE", "6_PRE", "7_PRE", "8_PRE",
                      "1_POST","2_POST","3_POST", "4_POST", "5_POST", "6_POST", "7_POST", "8_POST")

####run limma####
library(limma)

#change time for continuous
design=model.matrix(~VO2max_rel+Timepiont+Age,data=pheno) #time as continuous

corfit_PRE_Prot <- duplicateCorrelation(Prot1,design,block=pheno$ID)
#Then, we need to input this inter-correlation into the linear model to take into account the paired design
fit_PRE_Prot = lmFit(Prot1, design,maxit=1000,block=pheno$ID,correlation=corfit_PRE_Prot$consensus) #I have chosen the option "robust" which will downplay the influence of outliers. This takes significantly more time for your computer to run but is more robust and is always my preference.

#get info on eBayes in the vignette of the limma package
fit2_PRE_Prot <- eBayes(fit_PRE_Prot)

signif_VO2max_PRE_prot=limma::topTable(fit2_PRE_Prot, coef="VO2max_rel", adjust="BH",number=Inf,  p.value = 0.05)
signif_Timepoint_prot=limma::topTable(fit2_PRE_Prot, coef="Timepiont", adjust="BH",number=Inf,  p.value = 0.05)

Prot_VO2max_PRE=limma::topTable(fit2_PRE_Prot, coef="VO2max_rel", adjust="BH",number=Inf,  p.value = 1)
Timepoint_prot=limma::topTable(fit2_PRE_Prot, coef="Timepiont", adjust="BH",number=Inf,  p.value = 1)

SE <- fit2_PRE_Prot$sigma * fit2_PRE_Prot$stdev.unscaled
SE <- SE[rownames(Prot2),"Timepiont"]
Timepoint_prot <- cbind(Timepoint_prot,SE)

#Save results
Gene <- rownames(Timepoint_prot) #name of marker (CpG)
ALLELE1 <- rep(1,nrow(Timepoint_prot)) #fake allele info (not important)
ALLELE2 <- rep(2,nrow(Timepoint_prot)) #fake allele info (not important)
TESTSTAT <- Timepoint_prot$t #t-statistic
PVALUE <- Timepoint_prot$P.Value #p-value
EFFECTSIZE <- Timepoint_prot$logFC #effect size
SE <- Timepoint_prot$SE #standard error

Timepoint_prot = data.frame(Gene,
                            ALLELE1,
                            ALLELE2,
                            TESTSTAT,
                            PVALUE,
                            EFFECTSIZE,
                            SE)


setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Proteome/Metal Proteins")
write.table(Timepoint_prot,
            file="Proteomics_Timepoint_prot_atul.tbl",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")  



#Run Bacon

library(bacon)

#List all files ending in "tbl" in the folder
directory=c("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Proteome/Metal Proteins")

setwd(directory)
files <- list.files()[grep(".tbl",list.files())]

#Run each file in bacon using a loop
for (f in files)
{
  file <- read_tsv(f) #Read the file
  f <- sub("\\.tbl", "", f) #Obtain the name of the dataset without the ".tbl" at the end
  print(f) #show which dataset we are currently analysing
  bc <- bacon(teststatistics = NULL, #run bacon on effect sizes and standard errors
              effectsizes = file$EFFECTSIZE,
              standarderrors = file$SE)
  
  tiff(paste0('QQ-plot_',f,'.tiff'), #save the graph of raw and adjusted effect sizes and standard errors
       width =4,
       height = 2.5,
       units = 'in',
       res=600)
  print(plot(bc, type="qq")) #q-q plot shows the deviation of p-value distribution from null hypothesis. The observed P values for each CpG are sorted from largest to smallest and plotted against expected values from a theoretical ??2-distribution.
  dev.off()
  
  print(c(inflation(bc), #show the inflation factor for this dataset
          bias(bc))) #show the bias factor for this dataset
  
  file$EFFECTSIZE_CORR <- es(bc)[,1] #add a column to the file corresponding to the corrected effect size
  file$SE_CORR <- se(bc)[,1] #add a column to the file corresponding to the corrected standard error
  file$PVALUE_CORR <- pval(bc)[,1] #add a column to the file corresponding to the corrected p-value
  write.table(file, #save the file
              file=paste0(directory,"/",f,"_corrected.tbl"),
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE,
              sep="\t")
}

#Once this is done run METAL 
#Then load results for downstream analyses
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Proteome")
Meta_Proteomics=read.delim("METAANALYSIS1.tbl")

#Select only genes present on 1 or more studies
Meta_Proteomics_subset=Meta_Proteomics%>%
  filter(HetDf >=0)
Meta_Proteomics_subset=Meta_Proteomics_subset%>%
  mutate(FDR=p.adjust(Meta_Proteomics_subset$P.value, method = "fdr"))%>%
  mutate(t=Effect/StdErr)

signif_meta_Proteomics=Meta_Proteomics_subset[Meta_Proteomics_subset$FDR<0.05,]
signif_meta_Proteomics2=signif_meta_Proteomics[,c(1,4,5,6,12)]

write.csv(signif_meta_Proteomics2, "Metanalyses Proteomics significant results.csv")

sum(signif_meta_Proteomics$Effect >0)*100/nrow(signif_meta_Proteomics) #58.1
sum(signif_meta_Proteomics$Effect <0)*100/nrow(signif_meta_Proteomics) #41.9

#Build volcano plot with all the results
####PLOT 1: Volcano plot
library(ggplot2)
library(ggpubr)

#Add in different color the points that are significant
library(dplyr)
Proteins=Meta_Proteomics_subset$MarkerName

#Volcano Plot for p-values 0.005 threshold
#add a column of NAs
all_results = mutate(Meta_Proteomics_subset, diff="Not Sig")

# if Estimate > 0 and pvalue < 0.05, set as "UP" 
all_results$diff[all_results$Effect > 0.0 & all_results$FDR <0.05] <- "Up"
# if log2Foldchange < -0.0 and pvalue < 0.1, set as "DOWN"
all_results$diff[all_results$Effect < -0.0 & all_results$FDR <0.05] <- "Down"

all_results=all_results%>%
  dplyr::rename(logFC=Effect)

library(ggpubr)
library(ggrepel)

pp=ggplot(all_results, aes(logFC, -log10(`P.value`), col=diff)) +
  geom_point() +
  scale_color_manual(values=c("#1AEDD8","grey","#1BD2F7"))+
  labs(x=expression("Protein logFC per VO2max"), y="-log10(p-value)")+
  #geom_label_repel(data=filter(all_results[1:10,], adj.P.Val<0.05&!is.na(X)), aes(label=X), max.overlaps = 10)+
  theme_classic()+
  theme(legend.position = "none")

pp




tiff('PROT_meta analysis volcano.tiff',
     width =100,
     height = 80,
     units = 'mm',
     res=300)
pp
dev.off()

#Show histogram of effect sizes for DMPs
tiff('Distribution of effect size in protein_FINAL.tiff',
     width =100,
     height = 65,
     units = 'mm',
     res=300)
ggplot(data = all_results%>%dplyr::filter(FDR<0.05),
       aes(x=`logFC`,
           fill=diff)) +
  geom_histogram(colour = "grey")+
  scale_fill_manual(values=c("#1BD2F7","#1AEDD8"))+
  lims(x=c(min(Meta_Proteomics_subset$Effect),
           max(Meta_Proteomics_subset$Effect)))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

####Forest plot for top Proteins####
####Prep####
library(readxl)
dataset_summary <- read_excel('//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Proteome/PWAS datasets.xlsx',
                              na = c("","NA"))

#Create empty list
L.names <- all_results$MarkerName
L <- vector("list", length(L.names))
names(L) <- L.names

#Initiate list
tib <- tibble(Genes = all_results$MarkerName ,
              Dataset = rep("Meta-analysis",nrow(all_results)),
              ES = all_results$logFC,
              SE = all_results$StdErr,
              PVAL = all_results$`P.value`,
              FDR = all_results$FDR,
              Type = rep("Meta-analysis",nrow(all_results)))
tib <- tib %>%
  mutate_if(is.numeric,
            signif,
            digits = 2)

#Initiate the list from meta-analysis
L <- split(tib, seq(nrow(tib))) #splits each row of the tibble into a list component
names(L) <- all_results$MarkerName


colnames(dataset_summary)=c("Dataset", "Muscle type", "n")

#List tbl files
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Proteome/Metal Proteins")
files <- list.files(pattern=".tbl",
                    recursive = TRUE)

#Add each study to each component of the list
for (s in files)
{
  print(s)
  file <- read_tsv(s)
  
  #Add FDR to the table
  file <- file %>%
    mutate(FDR = p.adjust(PVALUE_CORR))
  
  #Add Effect 
  #Create tibble
  summary <- tibble(Genes = file$Gene,
                    Dataset = rep(s,nrow(file)),
                    ES = file$EFFECTSIZE_CORR,
                    SE = file$SE_CORR,
                    PVAL = file$PVALUE_CORR,
                    FDR = file$FDR,
                    Type = rep("Individual study",nrow(file)))
  summary <- summary %>%
    mutate_if(is.numeric,
              signif,
              digits = 2)
  
  summary <- summary %>%
    filter(Genes %in% names(L))
  subL <- split(summary, seq(nrow(summary)))
  names(subL) = summary$Genes
  
  #Merge the pieces of the two lists that are in common
  L2 <- L[names(subL)]
  L3 <- L[setdiff(names(L),names(subL))]
  listinter <- Map(bind_rows,
                   L2,
                   subL)
  L <- c(listinter,
         L3)
}

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Proteome/Metal Proteins/Tables and Figures")

saveRDS(L,
        file = "METAL_muscle_VO2max.rds")

####Forestplot - start here####
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Proteome/Metal Proteins/Tables and Figures")

L <- read_rds("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Proteome/Metal Proteins/Tables and Figures/METAL_muscle_VO2max.rds")

#1
top_down1 <- all_results %>%
  filter(logFC< 0) %>%
  arrange(FDR)%>%
  dplyr::slice(1)%>%
  pull(MarkerName) #PLCL1
top_down1=L[[top_down1]]
top_down1$Dataset=sub("\\..*","",top_down1$Dataset)

tib_down1 <- left_join(top_down1,
                       dataset_summary)


tib_down1=tib_down1 %>%
  dplyr::select("Genes","Dataset", "ES", "SE", "PVAL","FDR","n")%>%
  rename("logFC"="ES")

#2
top_down2 <- all_results %>%
  filter(logFC < 0) %>%
  arrange(FDR)%>%
  dplyr::slice(2)%>%
  pull(MarkerName) #OSBP
top_down2=L[[top_down2]]
top_down2$Dataset=sub("\\..*","",top_down2$Dataset)

tib_down2 <- left_join(top_down2,
                       dataset_summary)

tib_down2=tib_down2 %>%
  dplyr::select("Genes","Dataset", "ES", "SE", "PVAL","FDR","n")%>%
  rename("logFC"="ES")


#3
top_down3 <- all_results %>%
  filter(logFC < 0) %>%
  arrange(FDR)%>%
  dplyr::slice(3)%>%
  pull(MarkerName) #MTFR1L
top_down3=L[[top_down3]]
top_down3$Dataset=sub("\\..*","",top_down3$Dataset)

tib_down3 <- left_join(top_down3,
                       dataset_summary)%>%
  arrange(-`n`)

tib_down3=tib_down3 %>%
  dplyr::select("Genes","Dataset", "ES", "SE", "PVAL","FDR","n")%>%
  rename("logFC"="ES")

#4
top_down4 <-  all_results %>%
  filter(logFC < 0) %>%
  arrange(FDR)%>%
  dplyr::slice(4)%>%
  pull(MarkerName) #FXR1
top_down4=L[[top_down4]]
top_down4$Dataset=sub("\\..*","",top_down4$Dataset)

tib_down4 <- left_join(top_down4,
                       dataset_summary)%>%
  arrange(-`n`)

tib_down4=tib_down4 %>%
  dplyr::select("Genes","Dataset", "ES", "SE", "PVAL","FDR","n")%>%
  rename("logFC"="ES")

#5
top_down5 <- all_results %>%
  filter(logFC < 0) %>%
  arrange(FDR)%>%
  dplyr::slice(5)%>%
  pull(MarkerName) #RPL15
top_down5=L[[top_down5]]
top_down5$Dataset=sub("\\..*","",top_down5$Dataset)

tib_down5 <- left_join(top_down5,
                       dataset_summary)%>%
  arrange(-`n`)

tib_down5=tib_down5 %>%
  dplyr::select("Genes","Dataset", "ES", "SE", "PVAL","FDR","n")%>%
  rename("logFC"="ES")


top5_down=rbind(tib_down1,tib_down2,tib_down3,tib_down4,tib_down5)


#1
top_up1 <- all_results %>%
  filter(logFC > 0) %>%
  arrange(FDR)%>%
  dplyr::slice(1)%>%
  pull(MarkerName) #NTMT1
top_up1=L[[top_up1]]
top_up1$Dataset=sub("\\..*","",top_up1$Dataset)

tib_up1 <- left_join(top_up1,
                     dataset_summary)%>%
  arrange(-`n`)

tib_up1=tib_up1 %>%
  dplyr::select("Genes","Dataset", "ES", "SE", "PVAL","FDR","n")%>%
  rename("logFC"="ES")

#2
top_up2 <- all_results %>%
  filter(logFC > 0) %>%
  arrange(FDR)%>%
  dplyr::slice(2)%>%
  pull(MarkerName) #COQ10A
top_up2=L[[top_up2]]
top_up2$Dataset=sub("\\..*","",top_up2$Dataset)

tib_up2 <- left_join(top_up2,
                     dataset_summary)%>%
  arrange(-`n`)

tib_up2=tib_up2 %>%
  dplyr::select("Genes","Dataset", "ES", "SE", "PVAL","FDR","n")%>%
  rename("logFC"="ES")

#3
top_up3 <-all_results %>%
  filter(logFC > 0) %>%
  arrange(FDR)%>%
  dplyr::slice(3)%>%
  pull(MarkerName) #MRPS7
top_up3=L[[top_up3]]
top_up3$Dataset=sub("\\..*","",top_up3$Dataset)
tib_up3 <- left_join(top_up3,
                     dataset_summary)%>%
  arrange(-`n`)

tib_up3=tib_up3 %>%
  dplyr::select("Genes","Dataset", "ES", "SE", "PVAL","FDR","n")%>%
  rename("logFC"="ES")

#4
top_up4 <- all_results %>%
  filter(logFC > 0) %>%
  arrange(FDR)%>%
  dplyr::slice(4)%>%
  pull(MarkerName) #Immunoglobulin alpha-2 heavy chain
top_up4=L[[top_up4]]
top_up4$Dataset=sub("\\..*","",top_up4$Dataset)

tib_up4 <- left_join(top_up4,
                     dataset_summary)%>%
  arrange(-`n`)

tib_up4=tib_up4 %>%
  dplyr::select("Genes","Dataset", "ES", "SE", "PVAL","FDR","n")%>%
  rename("logFC"="ES")

#5
top_up5 <- all_results %>%
  filter(logFC > 0) %>%
  arrange(FDR)%>%
  dplyr::slice(5)%>%
  pull(MarkerName) #GARS
top_up5=L[[top_up5]]
top_up5$Dataset=sub("\\..*","",top_up5$Dataset)

tib_up5 <- left_join(top_up5,
                     dataset_summary)%>%
  arrange(-`n`)

tib_up5=tib_up5 %>%
  dplyr::select("Genes","Dataset", "ES", "SE", "PVAL","FDR","n")%>%
  rename("logFC"="ES")

top5_up=rbind(tib_up1,tib_up2,tib_up3,tib_up4,tib_up5)

top5_up<-top5_up%>%
  mutate(Dataset=factor(Dataset, levels= c("Meta-analysis", "Gene SMART", "PXD023084")))

fullTopCPGs=rbind(top5_up, top5_down)
fullTopCPGs<-fullTopCPGs%>%
  mutate(Dataset=factor(Dataset, levels= c("Meta-analysis", "Gene SMART", "PXD023084")))
  

tiff("Top5 Proteins.tiff",
     height = 11,
     width = 10,
     unit = 'in',
     res = 200)
ggforestplot::forestplot(
  df = fullTopCPGs,
  name = Genes,
  estimate = logFC,
  se = SE,
  pvalue = FDR,
  psignif = 0.05,
  colour = Dataset,
  shape = Dataset,
  #xlab = "Odds ratio per change in VO2max (95% CI) per 1−SD decrease in log2FC mRNA expression",
  title = "Top 5 down- and up- expressed Proteins associated with VO2max",
  logodds = F
)+
  ggplot2::scale_shape_manual(
    values = c(23L, 21L, 21L, 21L, 21L,21L,21L,21L,21L,21L),
    labels = c("Meta-analysis", "Gene SMART", "PXD023084")
  )+
  ggplot2::scale_colour_manual(values = c("#FF9A97","#1AEDD8","#1BD2F7"))
dev.off()

####pathway analysis####
library("mitch")
library(tidyverse)

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Proteome/Metal Proteins/Tables and Figures")
prot=read_tsv("METAANALYSIS1.TBL")

Prot=prot%>%
  mutate(FDR=p.adjust(prot$`P-value`, method = "fdr"))%>%
  filter(FDR<0.05)%>%
  dplyr::select("MarkerName", "Effect")%>%
  rename(Gene="MarkerName", EFFECTSIZE="Effect")

prot1 <- prot%>%
  mutate(t=Effect/StdErr)

colnames(prot1)[1] <- "gene"
head(prot1)
str(prot1)

download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip", destfile="ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
genesets <- gmt_import("ReactomePathways.gmt")

prot1 <- prot1[,c("gene","t")]
dim(prot1)

l <- list("prot"=prot1)
m <- mitch_import(x = l,geneIDcol = "gene",DEtype = "limma")


res <- mitch_calc(x=m,genesets = genesets ,priority = "significance",minsetsize = 5, resrows = 25)
head(res$enrichment_result,50)
mitch_report(res = res, outfile = "prot_significance_reactome.html")
mitch_plots(res,outfile="prot_significance_reactome.pdf")

#save results as table
res=res$enrichment_result
res_sig=res[res$p.adjustANOVA<0.05 ,]

#Reproduce plot using results table
res1 = mutate(res, sig=ifelse(res$p.adjustANOVA <0.05, "FDR<0.05", "Not Sig"))

sigdir <- res1$sig
sigdir[res1$sig=="FDR<0.05"&
         res1$s.dist <0]="Down"
sigdir[res1$sig=="FDR<0.05"&
         res1$s.dist>0]="Up"
res1=res1 %>%
  mutate(sigdir = sigdir)

write.table(res1, "results enrichment proteome reactome.txt")

library("ggrepel")

tiff('pathway VO2max proteome Reactome.tiff', width =10, height = 7, units = 'in', res=300)
pvsef=ggplot(res1, aes(x=`s.dist`, y=-log(`p.adjustANOVA`)))+
  geom_point(aes(col=sigdir))+
  scale_color_manual(values=c("#1BD2F7","gray","#1AEDD8"), name = " ")+
  theme_minimal()+
  labs(y="-log(p.adjustedANOVA)", x="s score")+
  geom_label_repel(data=filter(res1[1:10,]), aes(label=set), max.overlaps = 7, color=c("#1AEDD8","#1AEDD8", "#1AEDD8","#1AEDD8","#1AEDD8","#1AEDD8","#1BD2F7", "#1AEDD8","#1BD2F7" ,"#1BD2F7" ))
pvsef
dev.off()


