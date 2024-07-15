setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise/Methylation")

####Create function to format results for METAL software####
create_summary <- function(toptable = NULL,
                           dataset_label= NULL,
                           directory = getwd())
{
  CPG <- rownames(toptable)
  ALLELE1 <- rep(1,nrow(toptable))
  ALLELE2 <- rep(2,nrow(toptable))
  TESTSTAT <- toptable$t
  PVALUE <- toptable$P.Value
  EFFECTSIZE <- toptable$logFC
  SE <- toptable$SE
  results = data.frame(CPG,
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


####E-MTAB-11282####
#Load data (8 weeks of MICT)
B <- read.table("E-MTAB-11282_normalized_batch_corrected.txt")
library(minfi)
M <- logit2(B)
library(tidyverse)
pheno <- read_csv("E-MTAB-11282 phenotypes.csv")
nrow(pheno)
mean(pheno$Age)
sd(pheno$Age)
range(pheno$Age)
sum(pheno$Gender=="Male")/nrow(pheno)
pheno <- pheno %>%
  mutate(Timepoint = str_extract(`Sample ID`,
                                 "(?<= ).*"),
         ID = paste0("I",str_extract(`Sample ID`,
                                     ".*(?= )")))
sd(pheno$`VO2 max`)
mean(pheno$`VO2 max`)

#EWAS
design=model.matrix(~Age+
                      Gender*`VO2 max`+
                      Timepoint+
                      `Lean or Obese`,
                    pheno)
library(limma)
corfit <- duplicateCorrelation(M,
                               design,
                               block=pheno$ID)
cor = 0.150
fit1_M <- lmFit(M,
                design,
                block=pheno$ID,
                correlation=cor)
fit2_M <- eBayes(fit1_M)
fit1_B <- lmFit(B,
                design,
                block=pheno$ID,
                correlation=cor)
fit2_B <- eBayes(fit1_B)

coef = "TimepointPRE"
results <- topTable(fit2_M,
                    coef=coef,
                    number=Inf,
                    p.value=1)
results_B <- topTable(fit2_B,
                      coef=coef,
                      number=Inf,
                      p.value=1)
results$logFC <- results_B[rownames(results),"logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('E-MTAB-11282_pvalhist_DMPs.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Training DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "E-MTAB-11282_training",
               directory = getwd())

coef = "Age"
results <- topTable(fit2_M,
                    coef=coef,
                    number=Inf,
                    p.value=1)
results_B <- topTable(fit2_B,
                      coef=coef,
                      number=Inf,
                      p.value=1)
results$logFC <- results_B[rownames(results),"logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('E-MTAB-11282_pvalhist_DMPs_Age.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Age DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "E-MTAB-11282_Age",
               directory = getwd())


#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="E-MTAB-11282_M_res_AGE.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")


####Gene SMART 1####
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise/Methylation")

#Load data
M <- read.table("GeneSMART_M.txt")
library(minfi)
B <- ilogit2(M)
library(tidyverse)
pheno <- read.delim("GeneSMART_phenotypes.txt")
nrow(pheno)
mean(pheno$Age)
sd(pheno$Age)
range(pheno$Age)
sum(pheno$Sex=="M")/nrow(pheno)
mean(pheno$VO2max_rel, na.rm=T)
sd(pheno$VO2max_rel, na.rm=T)

#clean NA values
pheno=pheno[-234,]
M=M[,-234]
B=B[,-234]
pheno=pheno[-198,]
M=M[,-198]
B=B[,-198]
#replace missing VO2max
pheno[107,10]=37.8
pheno[121,10]=49.49
pheno[133,10]=41.95
pheno[184,10]=37.04

View(pheno)

#EWAS
design=model.matrix(~Age+
                      Sex*VO2max_rel+
                      Timepoint+
                      Batch+Intervention,
                    pheno)
library(limma)
corfit <- duplicateCorrelation(M,
                               design,
                               block=pheno$ID)
cor = 0.231
fit1_M <- lmFit(M,
                design,
                block=pheno$ID,
                correlation=cor)
fit2_M <- eBayes(fit1_M)
fit1_B <- lmFit(B,
                design,
                block=pheno$ID,
                correlation=cor)
fit2_B <- eBayes(fit1_B)

coef = "Timepoint4WP"
results <- topTable(fit2_M,
                    coef=coef,
                    number=Inf,
                    p.value=1)
results_B <- topTable(fit2_B,
                      coef=coef,
                      number=Inf,
                      p.value=1)
results$logFC <- results_B[rownames(results),"logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('Gene_SMART_pvalhist_DMPs_training.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for training DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "Gene SMART_training",
               directory = getwd())

coef = "Age"
results <- topTable(fit2_M,
                    coef=coef,
                    number=Inf,
                    p.value=1)
results_B <- topTable(fit2_B,
                      coef=coef,
                      number=Inf,
                      p.value=1)
results$logFC <- results_B[rownames(results),"logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('Gene_SMART_pvalhist_DMPs_Age.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Age DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "Gene_SMART_Age",
               directory = getwd())

#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="Gene_SMART_M_res_AGE.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")

####EXACT####
#Load data
B <- read.table("EXACT_normalized_batch_corrected.txt")
library(minfi)
M <- logit2(B)
library(tidyverse)
pheno <- read_csv("E192093_SampleSheet_EXACT_EPIC_48samples_230719.csv",
                  skip = 7)%>%
  mutate(ID = substr(Sample_Name,1,4))
nrow(pheno)
mean(pheno$Age)
sd(pheno$Age)
range(pheno$Age)
mean(pheno$VO2max)
sd(pheno$VO2max)


#EWAS
design=model.matrix(~Age+
                      Timepoint+
                      VO2max,
                    pheno)
library(limma)
corfit <- duplicateCorrelation(M,
                               design,
                               block=pheno$ID)
#cor = corfit$consensus.correlation #0.155
cor = 0.155
fit1_M <- lmFit(M,
                design,
                block=pheno$ID,
                correlation=cor)
fit2_M <- eBayes(fit1_M)
fit1_B <- lmFit(B,
                design,
                block=pheno$ID,
                correlation=cor)
fit2_B <- eBayes(fit1_B)

coef = "TimepointPost Exercise"
results <- topTable(fit2_M,
                    coef=coef,
                    number=Inf,
                    p.value=1)
results_B <- topTable(fit2_B,
                      coef=coef,
                      number=Inf,
                      p.value=1)
results$logFC <- results_B[rownames(results),"logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('EXACT_pvalhist_DMPs_Timepoint.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Timepoint DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "EXACT",
               directory = getwd())

coef = "Age"
results <- topTable(fit2_M,
                    coef=coef,
                    number=Inf,
                    p.value=1)
results_B <- topTable(fit2_B,
                      coef=coef,
                      number=Inf,
                      p.value=1)
results$logFC <- results_B[rownames(results),"logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('EXACT_pvalhist_DMPs_Age.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Age DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "EXACT_Age",
               directory = getwd())

#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="EXACT_M_res_AGE.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")


####GSE213029####
B <- read.table("GSE213029_normalized_batch_corrected_beta.txt")
library(minfi)
M <- logit2(B)
library(tidyverse)
pheno <- read_csv("GSE213029 phenotypes.csv")%>%
  mutate(Condition = str_replace(Condition,
                                 " .*$",
                                 ""))
nrow(pheno)
mean(pheno$Age)
sd(pheno$Age)
range(pheno$Age)
mean(pheno$VO2max)
sd(pheno$VO2max)

#EWAS
design=model.matrix(~Age+
                      Condition+
                      Time+
                      VO2max,
                    pheno)
library(limma)
corfit <- duplicateCorrelation(M,
                               design,
                               block=pheno$`Original Participant No`)
#cor = corfit$consensus.correlation #0.221
cor = 0.221
fit1_M <- lmFit(M,
                design,
                block=pheno$`Original Participant No`,
                correlation=cor)
fit2_M <- eBayes(fit1_M)
fit1_B <- lmFit(B,
                design,
                block=pheno$`Original Participant No`,
                correlation=cor)
fit2_B <- eBayes(fit1_B)

coef = "TimePre"
results <- limma::topTable(fit2_M,
                           coef=coef,
                           number=Inf,
                           p.value=1)
results_B <- limma::topTable(fit2_B,
                             coef=coef,
                             number=Inf,
                             p.value=1)
results$logFC <- results_B[rownames(results),"logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('GSE213029_pvalhist_timepoint_DMPs.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for timepoint DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "GSE213029_training",
               directory = getwd())

coef = "Age"
results <- topTable(fit2_M,
                    coef=coef,
                    number=Inf,
                    p.value=1)
results_B <- topTable(fit2_B,
                      coef=coef,
                      number=Inf,
                      p.value=1)
results$logFC <- results_B[rownames(results),"logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('GSE213029_pvalhist_DMPs_Age.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Age DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "GSE213029_Age",
               directory = getwd())

#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="GSE213029_M_res_Age.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")

####GSE60655####
B <- read.table("GSE60655_B.txt")
library(minfi)
M <- logit2(B)
library(tidyverse)
pheno <- read_delim("GSE60655_phenotypes.txt")

nrow(pheno)
mean(pheno$Age)
sd(pheno$Age)
range(pheno$Age)
sum(pheno$Sex=="Male")/nrow(pheno)
mean(pheno$VO2max)
sd(pheno$VO2max)

#EWAS
design=model.matrix(~Age+
                      Sex*VO2max+
                      Timepoint,
                    pheno)
library(limma)
corfit <- duplicateCorrelation(M,
                               design,
                               block=pheno$ID)
cor = corfit$consensus.correlation #0.216
cor = 0.216
fit1_M <- lmFit(M,
                design,
                block=pheno$ID,
                correlation=cor)
fit2_M <- eBayes(fit1_M)
fit1_B <- lmFit(B,
                design,
                block=pheno$ID,
                correlation=cor)
fit2_B <- eBayes(fit1_B)

coef = "TimepointT2"
results <- limma::topTable(fit2_M,
                           coef=coef,
                           number=Inf,
                           p.value=1)
results_B <- limma::topTable(fit2_B,
                             coef=coef,
                             number=Inf,
                             p.value=1)
results$logFC <- results_B[rownames(results),"logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('GSE60655_pvalhist_DMPs_training.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for training DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "GSE60655_timepoint",
               directory = getwd())

coef = "Age"
results <- topTable(fit2_M,
                    coef=coef,
                    number=Inf,
                    p.value=1)
results_B <- topTable(fit2_B,
                      coef=coef,
                      number=Inf,
                      p.value=1)
results$logFC <- results_B[rownames(results),"logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('GSE60655_pvalhist_DMPs_Age.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Age DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "GSE60655_Age",
               directory = getwd())

#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="GSE60655_M_res_AGE.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")

####trainSMART####
B <- read.table("beta trainsmart.txt")
library(minfi)
M <- logit2(B)
library(tidyverse)
library(readxl)
pheno <- read_xlsx("phenotype-train-smart.xlsx")%>%
  rename(Timepoint="Pre=T0_Post=T2")

#reorder beta and M table to match pheno ID
new_order=c(pheno$Sample_ID)
new_order=gsub("\\s",".", new_order)
B=B[,new_order]
M=M[,new_order]

nrow(pheno)
mean(pheno$Age)
sd(pheno$Age)
range(pheno$Age)
sum(pheno$Sex=="M")/nrow(pheno)
mean(pheno$VO2)
sd(pheno$VO2)

#EWAS
design=model.matrix(~Age+
                      Sex+
                      Timepoint,
                    pheno)
library(limma)
corfit <- duplicateCorrelation(M,
                               design,
                               block=pheno$CODE)
cor = corfit$consensus.correlation #0.1294
cor = 0.1294
fit1_M <- lmFit(M,
                design,
                block=pheno$CODE,
                correlation=cor)
fit2_M <- eBayes(fit1_M)
fit1_B <- lmFit(B,
                design,
                block=pheno$CODE,
                correlation=cor)
fit2_B <- eBayes(fit1_B)

coef = "TimepointT2"

results <- limma::topTable(fit2_M,
                           coef=coef,
                           number=Inf,
                           p.value=1)
results_B <- limma::topTable(fit2_B,
                             coef=coef,
                             number=Inf,
                             p.value=1)
results$logFC <- results_B[rownames(results),"logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('trainSMART_pvalhist_DMPs_training.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for timepoint DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "trainSMART_training",
               directory = getwd())

coef = "Age"
results <- topTable(fit2_M,
                    coef=coef,
                    number=Inf,
                    p.value=1)
results_B <- topTable(fit2_B,
                      coef=coef,
                      number=Inf,
                      p.value=1)
results$logFC <- results_B[rownames(results),"logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('trainSMART_pvalhist_DMPs_Age.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Age DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "trainSMART_Age",
               directory = getwd())


#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="trainSMART_M_res_Age.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")

####GSE114763####
M <- read.table("GSE114763_M.txt")
library(minfi)
B <- ilogit2(M)
library(tidyverse)
library(readxl)
pheno <- read.table("GSE114763_phenotypes.txt", row.names=NULL)
pheno2=pheno[c(3,5,8,10,13,15,17,20,22,24,26,29,31,34,36,38),]

pheno2=unite(pheno2, Sample_ID, "row.names", "Sample_Name", sep=".")

#reorder beta and M table to match pheno ID
new_order=c(pheno2$Sample_ID)
B=B[,new_order]
M=M[,new_order]

nrow(pheno)
mean(pheno$Age)
sd(pheno$Age)
range(pheno$Age)
sum(pheno$Sex=="M")/nrow(pheno)
mean(pheno$VO2)
sd(pheno$VO2)

#EWAS
design=model.matrix(~Age+
                      Timepoint,
                    pheno2)
library(limma)
corfit <- duplicateCorrelation(M,
                               design,
                               block=pheno2$ID)
cor = corfit$consensus.correlation #0.0829
cor = 0.1294
fit1_M <- lmFit(M,
                design,
                block=pheno2$ID,
                correlation=cor)
fit2_M <- eBayes(fit1_M)
fit1_B <- lmFit(B,
                design,
                block=pheno2$ID,
                correlation=cor)
fit2_B <- eBayes(fit1_B)

coef = "TimepointLOADING"

results <- limma::topTable(fit2_M,
                           coef=coef,
                           number=Inf,
                           p.value=1)
results_B <- limma::topTable(fit2_B,
                             coef=coef,
                             number=Inf,
                             p.value=1)
results$logFC <- results_B[rownames(results),"logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('GSE114763_pvalhist_DMPs_training.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for timepoint DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "GSE114763_training",
               directory = getwd())


####Deakin####
B <- read.table("beta Deakin.txt")
library(minfi)
M <- logit2(B)
library(tidyverse)
library(readxl)
pheno <- read_xlsx("phenotype-deakin.xlsx")
pheno<-pheno%>%
  separate(Sample_ID, c("Sample_ID", "Timepoint"), " ")%>%
  mutate(Timepoint=fct_relevel(Timepoint, "PRE", "POST"))


#EWAS
design=model.matrix(~Age+
                      Timepoint,
                    pheno)
library(limma)
corfit <- duplicateCorrelation(M,
                               design,
                               block=pheno$Sample_ID)
cor = corfit$consensus.correlation #0.07732
cor = 0.07732
fit1_M <- lmFit(M,
                design,
                block=pheno$Sample_ID,
                correlation=cor)
fit2_M <- eBayes(fit1_M)
fit1_B <- lmFit(B,
                design,
                block=pheno$Sample_ID,
                correlation=cor)
fit2_B <- eBayes(fit1_B)

coef = "TimepointPOST"

results <- limma::topTable(fit2_M,
                           coef=coef,
                           number=Inf,
                           p.value=1)
results_B <- limma::topTable(fit2_B,
                             coef=coef,
                             number=Inf,
                             p.value=1)
results$logFC <- results_B[rownames(results),"logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('Deakin_pvalhist_DMPs_training.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for timepoint DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "Deakin_training",
               directory = getwd())

coef = "Age"
results <- topTable(fit2_M,
                    coef=coef,
                    number=Inf,
                    p.value=1)
results_B <- topTable(fit2_B,
                      coef=coef,
                      number=Inf,
                      p.value=1)
results$logFC <- results_B[rownames(results),"logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('Deakin_pvalhist_DMPs_Age.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Age DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "Deakin_Age",
               directory = getwd())


#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="Deakin_M_res_Age.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")

####epiH####
B <- read.table("EpiH Normalized beta values after batch correction.txt")
library(minfi)
M <- logit2(B)
library(tidyverse)
library(readxl)
pheno <- read_xlsx("Phenotype data - EPI-H samples -For Andrew.xlsx")

pheno=pheno%>%
  mutate(Timepoint=fct_relevel(Condition, "pre", "post"))

#reorder beta and M table to match pheno ID
new_order=c(pheno$ID)
B=B[,new_order]
M=M[,new_order]


nrow(pheno)
mean(pheno$`Age at inclusion (years)`)
sd(pheno$`Age at inclusion (years)`)
range(pheno$`Age at inclusion (years)`)
sum(pheno$Sex=="M")/nrow(pheno)
mean(pheno$`Fitness level (ml/min/kg)`)
sd(pheno$`Fitness level (ml/min/kg)`)


#EWAS
design=model.matrix(~`Age at inclusion (years)`+
                      Sex*`Fitness level (ml/min/kg)`+
                      Timepoint,
                    pheno)
library(limma)
corfit <- duplicateCorrelation(M,
                               design,
                               block=pheno$ID_rep)
cor = corfit$consensus.correlation #0.1
cor=0.1

fit1_M <- lmFit(M,
                design,
                block=pheno$ID,
                correlation=cor)
fit2_M <- eBayes(fit1_M)
fit1_B <- lmFit(B,
                design,
                block=pheno$ID,
                correlation=cor)
fit2_B <- eBayes(fit1_B)

coef = "Timepointpost"

results <- limma::topTable(fit2_M,
                           coef=coef,
                           number=Inf,
                           p.value=1)
results_B <- limma::topTable(fit2_B,
                             coef=coef,
                             number=Inf,
                             p.value=1)
results$logFC <- results_B[rownames(results),"logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('epiH_pvalhist_DMPs_training.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for training DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "epiH_training",
               directory = getwd())

coef = "`Age at inclusion (years)`"
results <- topTable(fit2_M,
                    coef=coef,
                    number=Inf,
                    p.value=1)
results_B <- topTable(fit2_B,
                      coef=coef,
                      number=Inf,
                      p.value=1)
results$logFC <- results_B[rownames(results),"logFC"]
SE <- fit2_B$sigma * fit2_B$stdev.unscaled
results$SE <- SE[rownames(results),coef]

#Save p-value histogram
tiff('EpiH_pvalhist_DMPs_Age.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Age DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "EpiH_Age",
               directory = getwd())


#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="epiH_M_res_Age.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")
####Bacon - inflation correction####
#List all files ending in "tbl" in the folder
directory=c("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise/Methylation/Metal exercise/repeatbacon")
setwd(directory)
library(bacon)
library(tidyverse)

#List tbl files
files <- list.files()[grep(".tbl",list.files())]

for (f in files)
{
  file <- read_tsv(f) #Read the file
  f <- sub("\\.tbl", "", f) #Obtain the name of the dataset without the ".tbl" at the end
  print(f)
  
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
          bias(bc))) 
  
  file$EFFECTSIZE_CORR <- es(bc)[,1]
  file$SE_CORR <- se(bc)[,1]
  file$PVALUE_CORR <- pval(bc)[,1]
  
  #Multiply ES and SE by 100 because METAL rounds very small numbers
  File <- file %>%
    mutate(EFFECTSIZE_CORR100=(EFFECTSIZE_CORR*100))%>%
    mutate(SE_CORR100=(SE_CORR*100))
  
  write.table(File, #save the file
              file=paste0(directory,"/",f,"_corrected.tbl"),
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE,
              sep="\t")
}


#Write code for METAL
#List tbl files
files <- list.files(pattern=".tbl",
                    recursive = TRUE)
file <- "METAL_commands.txt"
write.table(paste("SCHEME STDERR",
                  "MARKER CPG",
                  "EFFECT EFFECTSIZE_CORR100",
                  "SEPARATOR TAB",
                  "STDERR SE_CORR100",
                  "PVALUE PVALUE_CORR",
                  sep = "\r"),
            file = file,
            quote=F,
            row.names=F,
            col.names=F)
#Add datasets to analyse
add_datasets_METAl <- function(d = NULL,
                               f = NULL)
{
  write(paste("PROCESS",
              d,
              sep = "\t"),
        file=f,
        append=TRUE)
}
sapply(files,
       add_datasets_METAl,
       f = file)
#Add last line of code
write("ANALYZE HETEROGENEITY",
      file=file,
      append=TRUE)

####Load tables to change annotation so it matched between arrays####
#Load annotation
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Annotation")
annotation <- read.delim("Annotation.txt")

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise/Methylation/Metal exercise")
E_MTAB <- read_tsv("E-MTAB-11282_training_corrected.TBL")
E_MTAB2 <- left_join(E_MTAB, annotation, join_by(CPG==probeID))
E_MTAB2=E_MTAB2[!duplicated(E_MTAB2$Unique_ID),]
E_MTAB2=E_MTAB2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(E_MTAB2, #save the file
            file="E-MTAB-11282_training_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

Deakin <- read_tsv("Deakin_training_corrected.TBL")
Deakin2 <- left_join(Deakin, annotation, join_by(CPG==probeID))
Deakin2=Deakin2[!duplicated(Deakin2$Unique_ID),]
Deakin2=Deakin2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(Deakin2, #save the file
            file="Deakin_training_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

EXACT <- read_tsv("EXACT_training_corrected.TBL")
EXACT2 <- left_join(EXACT, annotation, join_by(CPG==probeID))
EXACT2=EXACT2[!duplicated(EXACT2$Unique_ID),]
EXACT2=EXACT2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(EXACT2, #save the file
            file="EXACT_training_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

GeneSMART <- read_tsv("Gene SMART_training_corrected.TBL")
GeneSMART2 <- left_join(GeneSMART, annotation, join_by(CPG==probeID))
GeneSMART2=GeneSMART2[!duplicated(GeneSMART2$Unique_ID),]
GeneSMART2=GeneSMART2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(GeneSMART2, #save the file
            file="Gene SMART_training_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

GSE60655 <- read_tsv("GSE60655_timepoint_corrected.TBL")
GSE60655_2 <- left_join(GSE60655, annotation, join_by(CPG==probeID))
GSE60655_2=GSE60655_2[!duplicated(GSE60655_2$Unique_ID),]
GSE60655_2=GSE60655_2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(GSE60655_2, #save the file
            file="GSE60655_timepoint_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

GSE114763 <- read_tsv("GSE114763_training_corrected.TBL")
GSE114763_2 <- left_join(GSE114763, annotation, join_by(CPG==probeID))
GSE114763_2=GSE114763_2[!duplicated(GSE114763_2$Unique_ID),]
GSE114763_2=GSE114763_2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(GSE114763_2, #save the file
            file="GSE114763_training_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

GSE213029 <- read_tsv("GSE213029_training_corrected.TBL")
GSE213029_2 <- left_join(GSE213029, annotation, join_by(CPG==probeID))
GSE213029_2=GSE213029_2[!duplicated(GSE213029_2$Unique_ID),]
GSE213029_2=GSE213029_2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(GSE213029_2, #save the file
            file="GSE213029_training_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

trainSMART <- read_tsv("trainSMART_training_corrected.TBL")
trainSMART_2 <- left_join(trainSMART, annotation, join_by(CPG==probeID))
trainSMART_2=trainSMART_2[!duplicated(trainSMART_2$Unique_ID),]
trainSMART_2=trainSMART_2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(trainSMART_2, #save the file
            file="trainSMART_training_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

EPIK <- read_tsv("EPIK_training_corrected.TBL")
EPIK_2 <- left_join(EPIK, annotation, join_by(CPG==probeID))
EPIK_2=EPIK_2[!duplicated(EPIK_2$Unique_ID),]
EPIK_2=EPIK_2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(EPIK_2, #save the file
            file="EPIK_training_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

EpiH <- read_tsv("epiH_training_corrected.TBL")
EpiH2 <- left_join(EpiH, annotation2, join_by(CPG==probeID))
EpiH2=EpiH2[!duplicated(EpiH2$Unique_ID),]
EpiH2=EpiH2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(EpiH2, #save the file
            file="epiH_training_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

annotation2=unite(annotation, "probeID", "probeID":"EPICv2")


trainSMART <- read_tsv("trainSMART_training_corrected.TBL")
trainSMART2 <- left_join(trainSMART, annotation2, join_by(CPG==probeID))
trainSMART2=trainSMART2[!duplicated(trainSMART2$Unique_ID),]
trainSMART2=trainSMART2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(trainSMART2, #save the file
            file="trainSMART_training_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")


Deakin <- read_tsv("Deakin_training_corrected.TBL")
Deakin2 <- left_join(Deakin, annotation2, join_by(CPG==probeID))
Deakin2=Deakin2[!duplicated(Deakin2$Unique_ID),]
Deakin2=Deakin2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(Deakin2, #save the file
            file="Deakin_training_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

####Load METAL results and annotate####
library(tidyverse)
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets exercise/Methylation/METAL exercise")
meta_res_aerobic <- read_tsv("METAANALYSIS_Aerobic.TBL")
meta_res_resistance <- read_tsv("METAANALYSIS_Resistance.TBL")

#Change name to probeID
meta_res_aerobic <- dplyr::rename(meta_res_aerobic,
                          Unique_ID = MarkerName)

meta_res_resistance <- dplyr::rename(meta_res_resistance,
                                  Unique_ID = MarkerName)

#Add CHR and position to meta analysis results
meta_res_Aerobic_tot <- left_join(x = meta_res_aerobic,
                           y = annotation,
                           by = "Unique_ID")
meta_res_Aerobic_tot<-meta_res_Aerobic_tot[!duplicated(meta_res_Aerobic_tot$Unique_ID), ]

max(meta_res_Aerobic_tot$HetDf)


#Add CHR and position to meta analysis results
meta_res_resistance_tot <- left_join(x = meta_res_resistance,
                                  y = annotation,
                                  by = "Unique_ID")
meta_res_resistance_tot<-meta_res_resistance_tot[!duplicated(meta_res_resistance_tot$Unique_ID), ]

max(meta_res_resistance_tot$HetDf)

#Calculate nb of CpGs included in analysis for each nb of studies
nb_CpGs <- c()
for (i in 1:max(meta_res_Aerobic_tot$HetDf))
{
  nb_CpGs <- c(nb_CpGs,nrow(meta_res_Aerobic_tot %>% filter(HetDf>=(i-1))))
}

cumsum <- tibble(`Number of included studies` = 1:max(meta_res_Aerobic_tot$HetDf),
                 `Number of CpGs present in all studies` = nb_CpGs)

ggplot(cumsum,
       mapping = aes(x = `Number of included studies`,
                     y = `Number of CpGs present in all studies`))+
  geom_point()+
  geom_line()+
  theme_classic()

#Calculate nb of CpGs included in analysis for each nb of studies
nb_CpGs <- c()
for (i in 1:max(meta_res_resistance_tot$HetDf))
{
  nb_CpGs <- c(nb_CpGs,nrow(meta_res_resistance_tot %>% filter(HetDf>=(i-1))))
}

cumsum <- tibble(`Number of included studies` = 1:max(meta_res_resistance_tot$HetDf),
                 `Number of CpGs present in all studies` = nb_CpGs)

ggplot(cumsum,
       mapping = aes(x = `Number of included studies`,
                     y = `Number of CpGs present in all studies`))+
  geom_point()+
  geom_line()+
  theme_classic()

#View results including in at least 1 study
meta_res_Aerobic_robust <- meta_res_Aerobic_tot %>% filter(HetDf>=3)
padj <- p.adjust(meta_res_Aerobic_robust$`P-value`,method="BH")
meta_res_Aerobic_robust <- meta_res_Aerobic_robust %>%
  mutate(`Adjusted p-value` = padj)
meta_res_Aerobic_robust <- meta_res_Aerobic_robust %>%
  mutate(sig = ifelse(`Adjusted p-value`<0.05,
                      "Sig",
                      "Not Sig"))
nrow(meta_res_Aerobic_robust) #635923
nDMPs <- nrow(meta_res_Aerobic_robust %>% filter(sig=="Sig")) #DMPs = 66012

#View results including in at least 1 study
meta_res_resistance_robust <- meta_res_resistance_tot %>% filter(HetDf>=0)
padj <- p.adjust(meta_res_resistance_robust$`P-value`,method="BH")
meta_res_resistance_robust <- meta_res_resistance_robust %>%
  mutate(`Adjusted p-value` = padj)
meta_res_resistance_robust <- meta_res_resistance_robust %>%
  mutate(sig = ifelse(`Adjusted p-value`<0.05,
                      "Sig",
                      "Not Sig"))
nrow(meta_res_resistance_robust) #950907
nDMPs <- nrow(meta_res_resistance_robust %>% filter(sig=="Sig")) #DMPs = 140


#####Arrange
sigdir <- meta_res_Aerobic_robust$sig
sigdir[meta_res_Aerobic_robust$sig=="Sig"&
         meta_res_Aerobic_robust$Effect<0]="Hypo"
sigdir[meta_res_Aerobic_robust$sig=="Sig"&
         meta_res_Aerobic_robust$Effect>0]="Hyper"
meta_res_Aerobic_robust <- meta_res_Aerobic_robust %>%
  mutate(sigdir = sigdir)
meta_res_Aerobic_robust <- meta_res_Aerobic_robust %>%
  arrange(desc(sigdir))

nrow(meta_res_Aerobic_robust %>% filter(sigdir=="Hypo"))/nDMPs #% hypo DMPs = 47% hypo

#Select relevant columns only and rename them
meta_res_Aerobic_robust <- meta_res_Aerobic_robust %>%
  dplyr::select(probeID,
                CpG_chrm,
                CpG_beg,
                CGIposition,
                GeneHancer_interaction,
                genesUniq,
                Effect,
                StdErr,
                `P-value`,
                `Adjusted p-value`,
                HetDf,
                HetISq,
                HetPVal,
                sig,
                sigdir,
                Unique_ID)
meta_res_Aerobic_robust <- meta_res_Aerobic_robust %>%
  mutate(genesUniq = replace_na(genesUniq,""),
         GeneHancer_interaction = replace_na(GeneHancer_interaction,""))

#Replace CpG island names with better words
meta_res_Aerobic_robust$CGIposition[grep("Shore",meta_res_Aerobic_robust$CGIposition)]="Shore"
meta_res_Aerobic_robust$CGIposition[grep("Shelf",meta_res_Aerobic_robust$CGIposition)]="Shelf"
meta_res_Aerobic_robust$CGIposition[is.na(meta_res_Aerobic_robust$CGIposition)]="Open sea"

#Change HetDf by number of included studies
meta_res_Aerobic_robust$HetDf <- meta_res_Aerobic_robust$HetDf +1
#Rename columns
colnames(meta_res_Aerobic_robust) = c("CpG",
                              "Chromosome",
                              "Position (hg38)",
                              "CpG island position",
                              "GeneHancer interaction",
                              "Annotated gene(s)",
                              "Effect size",
                              "SE",
                              "P-value",
                              "FDR",
                              "Number of studies",
                              "Heterogeneity index (I2)",
                              "Heterogeneity p-value",
                              "Significance",
                              "Direction",
                              "Unique ID")

nDMPs <- nrow(meta_res_Aerobic_robust %>% filter(Significance=="Sig")) #3100

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise/Methylation/Metal exercise")
write.table(meta_res_Aerobic_robust,
            file="METAL_muscle_Aerobic.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")
meth_aerobic=read.delim("METAL_muscle_Aerobic.txt")

#####Arrange
sigdir <- meta_res_resistance_robust$sig
sigdir[meta_res_resistance_robust$sig=="Sig"&
         meta_res_resistance_robust$Effect<0]="Hypo"
sigdir[meta_res_resistance_robust$sig=="Sig"&
         meta_res_resistance_robust$Effect>0]="Hyper"
meta_res_resistance_robust <- meta_res_resistance_robust %>%
  mutate(sigdir = sigdir)
meta_res_resistance_robust <- meta_res_resistance_robust %>%
  arrange(desc(sigdir))

nrow(meta_res_resistance_robust %>% filter(sigdir=="Hypo"))/nDMPs #% hypo DMPs = 61% hypo

#Select relevant columns only and rename them
meta_res_resistance_robust <- meta_res_resistance_robust %>%
  dplyr::select(probeID,
                CpG_chrm,
                CpG_beg,
                CGIposition,
                GeneHancer_interaction,
                genesUniq,
                Effect,
                StdErr,
                `P-value`,
                `Adjusted p-value`,
                HetDf,
                HetISq,
                HetPVal,
                sig,
                sigdir,
                Unique_ID)
meta_res_resistance_robust <- meta_res_resistance_robust %>%
  mutate(genesUniq = replace_na(genesUniq,""),
         GeneHancer_interaction = replace_na(GeneHancer_interaction,""))

#Replace CpG island names with better words
meta_res_resistance_robust$CGIposition[grep("Shore",meta_res_resistance_robust$CGIposition)]="Shore"
meta_res_resistance_robust$CGIposition[grep("Shelf",meta_res_resistance_robust$CGIposition)]="Shelf"
meta_res_resistance_robust$CGIposition[is.na(meta_res_resistance_robust$CGIposition)]="Open sea"

#Change HetDf by number of included studies
meta_res_resistance_robust$HetDf <- meta_res_resistance_robust$HetDf +1
#Rename columns
colnames(meta_res_resistance_robust) = c("CpG",
                                      "Chromosome",
                                      "Position (hg38)",
                                      "CpG island position",
                                      "GeneHancer interaction",
                                      "Annotated gene(s)",
                                      "Effect size",
                                      "SE",
                                      "P-value",
                                      "FDR",
                                      "Number of studies",
                                      "Heterogeneity index (I2)",
                                      "Heterogeneity p-value",
                                      "Significance",
                                      "Direction",
                                      "Unique ID")

nDMPs <- nrow(meta_res_resistance_robust %>% filter(Significance=="Sig")) #140

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise/Methylation/Metal exercise")
write.table(meta_res_resistance_robust,
            file="METAL_muscle_resistance.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")
meth_resistance=read.delim("METAL_muscle_resistance.txt")

#### Tables & Figures ####
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise/Methylation/Tables & Figures")
library(tidyverse)
#Load METAL results that have been formatted to annotate each CpG
to_write_A <- meta_res_Aerobic_robust
to_order_A <- tibble(CpG = to_write_A$CpG)

to_write_R <- meta_res_resistance_robust
to_order_R <- tibble(CpG = to_write_R$CpG)

to_write_A <- to_write_A %>%
  dplyr::rename(`% DNAm change after aerobic training` = `Effect size`,
                `P-value Aerobic` = `P-value`,
                `FDR Aerobic` = FDR)%>%
  arrange(`P-value Aerobic`)

DMPs_A <- subset(to_write_A,
               `FDR Aerobic`<0.05)

to_write_R <- to_write_R %>%
  dplyr::rename(`% DNAm change after resistance training` = `Effect size`,
                `P-value Resistance` = `P-value`,
                `FDR Resistance` = FDR)%>%
  arrange(`P-value Resistance`)

DMPs_R <- subset(to_write_R,
                 `FDR Resistance`<0.05)

#Volcano plot for EWAS RE meta-analysis of Aerobic training

tiff('Volcano plot_EWAS RE meta-analysis of Aerobic.tiff',
     width =100,
     height = 80,
     units = 'mm',
     res=300)
ggplot() +
  geom_point(data = to_write_A,
             aes(`% DNAm change after aerobic training`,
                 -log10(`P-value Aerobic`),
                 col=Direction),
             size=0.5)+
  scale_color_manual(values=c("#1BD2F7", "#1AEDD8","grey"))+
  labs(y="-log10(p-value)")+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

#Volcano plot for EWAS RE meta-analysis of Resistance training

tiff('Volcano plot_EWAS RE meta-analysis of Resistance.tiff',
     width =100,
     height = 80,
     units = 'mm',
     res=300)
ggplot() +
  geom_point(data = to_write_R,
             aes(`% DNAm change after resistance training`,
                 -log10(`P-value Resistance`),
                 col=Direction),
             size=0.5)+
  scale_color_manual(values=c("#1BD2F7", "#1AEDD8","grey"))+
  labs(y="-log10(p-value)")+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

#### Integration ####
library(tidyverse)

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise/Methylation/METAL exercise")
meth_Aerobic=read.delim("METAL_muscle_Aerobic.txt")
meth_Resistance=read.delim("METAL_muscle_Resistance.txt")

####Enrichment with cluster profiler####
# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(pheatmap)
library(clusterProfiler) # for PEA analysis
library('org.Hs.eg.db')
library(DOSE)
library(enrichplot) # for visualisations
library(ggupset) # for visualisations

## Function: Adjacency matrix to list ####
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

#Methylation Resistance
meth_R1 <- meth_Resistance%>%
  mutate(t=Effect.size/SE)%>%
  dplyr::select("CpG","Effect.size","SE","P.value","FDR", "t","Direction","Annotated.gene.s.")
dim(meth_R1)
methg_R <- meth_R1$`Annotated.gene.s.`
methg_spl_R <- strsplit(methg_R,";")
meth_R1 <- meth_R1[which(lapply(methg_spl_R,length)>0),]
head(meth_R1)
dim(meth_R1)
x <- apply(meth_R1,1,function(x) {
  GENES <- as.character(x[length(x)])
  GENES <- unlist(strsplit(GENES,";"))
  z <- lapply(GENES,function(G) {
    y <- x[1:(length(x)-1)] 
    y$Gene <- G
    y
  })
  do.call(rbind,z)
})
meth_fix <- do.call(rbind,x)
meth_fix2 <- as.data.frame(apply(meth_fix[,c(2:6)],2,as.numeric))
meth_fix2$gene <- unlist(meth_fix[,"Gene"])
meth_fix2$Direction <- unlist(meth_fix[,"Direction"])
head(meth_fix2)
dim(meth_fix2)
meth_R1 <- meth_fix2[,c("gene","t")]
meth_R1 <- aggregate(. ~ gene, meth_R1, sum)
dim(meth_R1)

# Get the genes that are present in your dataframe
genes_in_meth_R <- meth_fix2$gene
# Read in the .gmt file
pwl2 <- read.gmt("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/ReactomePathways.gmt") 
# Subset to the genes that are present in our dataset
pwl2_R <- pwl2[pwl2$gene %in% genes_in_meth_R,] 
# Save the filtered background gene set
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise")
filename <- 'meth_path_R.RDS'
saveRDS(pwl2_R, filename)

# Remove non-significant genes
meth_path_R <- meth_fix2[meth_fix2$Direction != 'Not Sig', ]
# Substitute names so they are annotated nicely in the heatmap later
# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
meth_results_list_R <- split(meth_path_R, meth_path_R$Direction)

## Run ClusterProfiler -----------------------------------------------
# Settings
name_of_comparison <- 'Resistance' # for our filename
background_genes <- 'reactome' # for our filename
bg_genes_R <- readRDS("meth_path_R.RDS") # read in the background genes
padj_cutoff <- 0.05 # p-adjusted threshold, used to filter out pathways
genecount_cutoff <- 5 # minimum number of genes in the pathway, used to filter out pathways

res_R <- lapply(names(meth_results_list_R),
              function(x) enricher(gene = meth_results_list_R[[x]]$gene,
                                   TERM2GENE = bg_genes_R))
names(res_R) <- names(meth_results_list_R)
#Convert the list of enrichResults for each sample_pattern to a dataframe with the pathways
res_df_R <- lapply(names(res_R), function(x) rbind(res_R[[x]]@result))
names(res_df_R) <- names(res_R)
res_df_R <- do.call(rbind, res_df_R)
head(res_df_R)

#Convert the enrichResults to a dataframe with the pathways
res_df_R <- lapply(names(res_R), function(x) rbind(res_R[[x]]@result))
names(res_df_R) <- names(res_R)
res_df_R <- do.call(rbind, res_df_R)
head(res_df_R)

res_df_R <- res_df_R %>% mutate(minuslog10padj = -log10(p.adjust),
                            diffexpressed = gsub('\\.REACTOME.*$', '', rownames(res_df_R)))

# Subset to those pathways that have p adj < cutoff and gene count > cutoff (you can also do this in the enricher function)
target_pws <- unique(res_df_R$ID[res_df_R$p.adjust < padj_cutoff ]) # select only target pathways have p adjusted < 0.05 and at least 6 genes
res_df_R <- res_df_R[res_df_R$ID %in% target_pws, ]

write.csv(res_df_R, 'meth_resclusterp_Resistance.csv', row.names = FALSE)

res_df_down_R <- res_df_R %>% filter(grepl("Hypo", diffexpressed)) %>% 
  dplyr::select(!c('minuslog10padj', 'diffexpressed')) 
rownames(res_df_down_R) <- res_df_down_R$ID

# For visualisation purposes, let's shorten the pathway names
res_df_down_R$Description <- gsub('(H|h)iv', 'HIV', 
                              gsub('pd 1', 'PD-1',
                                   gsub('ecm', 'ECM', 
                                        gsub('(I|i)nterleukin', 'IL', 
                                             gsub('(R|r)na', 'RNA', 
                                                  gsub('(D|d)na', 'DNA',
                                                       gsub(' i ', ' I ', 
                                                            gsub('(A|a)tp ', 'ATP ', 
                                                                 gsub('(N|n)adh ', 'NADH ', 
                                                                      gsub('(N|n)ad ', 'NAD ',
                                                                           gsub('t cell', 'T cell',
                                                                                gsub('b cell', 'B cell',
                                                                                     gsub('built from .*', ' (...)',
                                                                                          gsub('mhc', 'MHC',
                                                                                               gsub('mhc class i', 'MHC I', 
                                                                                                    gsub('mhc class ii', 'MHC II', 
                                                                                                         stringr::str_to_sentence(
                                                                                                           gsub('_', ' ',  
                                                                                                                gsub('REACTOME_', '', res_df_down_R$Description)))))))))))))))))))

enrichres_down_R <- new("enrichResult",
                    readable = FALSE,
                    result = res_df_down_R,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.2,
                    organism = "human",
                    ontology = "UNKNOWN",
                    gene = meth_fix2$gene,
                    keytype = "UNKNOWN",
                    universe = unique(bg_genes_R$gene),
                    gene2Symbol = character(0),
                    geneSets = bg_genes_R)
class(enrichres_down)


# Dotplot
tiff('dot plot hypo pathways Resistance.tiff',
     width =200,
     height = 300,
     units = 'mm',
     res=300)
dotplot(enrichres_down_R, showCategory = 20) + ggtitle("Top pathways from hypo methylated DMRs")
dev.off()

tiff('dot plot hyper pathways.tiff',
     width =200,
     height = 300,
     units = 'mm',
     res=300)
dotplot(enrichres_up_A, showCategory = 20) + ggtitle("Top pathways from hyper methylated DMRs")
dev.off()


#Methylation Aerobic
meth_A1 <- meth_Aerobic%>%
  mutate(t=Effect.size/SE)%>%
  dplyr::select("CpG","Effect.size","SE","P.value","FDR", "t","Direction","Annotated.gene.s.")
dim(meth_A1)
methg_A <- meth_A1$`Annotated.gene.s.`
methg_spl_A <- strsplit(methg_A,";")
meth_A1 <- meth_A1[which(lapply(methg_spl_A,length)>0),]
head(meth_A1)
dim(meth_A1)
x <- apply(meth_A1,1,function(x) {
  GENES <- as.character(x[length(x)])
  GENES <- unlist(strsplit(GENES,";"))
  z <- lapply(GENES,function(G) {
    y <- x[1:(length(x)-1)] 
    y$Gene <- G
    y
  })
  do.call(rbind,z)
})
meth_fix <- do.call(rbind,x)
meth_fix2 <- as.data.frame(apply(meth_fix[,c(2:6)],2,as.numeric))
meth_fix2$gene <- unlist(meth_fix[,"Gene"])
meth_fix2$Direction <- unlist(meth_fix[,"Direction"])
head(meth_fix2)
dim(meth_fix2)
meth_A1 <- meth_fix2[,c("gene","t")]
meth_A1 <- aggregate(. ~ gene, meth_A1, sum)
dim(meth_A1)

# Get the genes that are present in your dataframe
genes_in_meth_A <- meth_fix2$gene
# Read in the .gmt file
pwl2 <- read.gmt("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/ReactomePathways.gmt") 
# Subset to the genes that are present in our dataset
pwl2_A <- pwl2[pwl2$gene %in% genes_in_meth_A,] 
# Save the filtered background gene set
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise")
filename <- 'meth_path_A.RDS'
saveRDS(pwl2_A, filename)

# Remove non-significant genes
meth_path_A <- meth_fix2[meth_fix2$Direction != 'Not Sig', ]
# Substitute names so they are annotated nicely in the heatmap later
# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
meth_results_list_A <- split(meth_path_A, meth_path_A$Direction)

## Run ClusterProfiler -----------------------------------------------
# Settings
name_of_comparison <- 'Aerobic' # for our filename
background_genes <- 'reactome' # for our filename
bg_genes_A <- readRDS("meth_path_A.RDS") # read in the background genes
padj_cutoff <- 0.05 # p-adjusted threshold, used to filter out pathways
genecount_cutoff <- 5 # minimum number of genes in the pathway, used to filter out pathways

res_A <- lapply(names(meth_results_list_A),
                function(x) enricher(gene = meth_results_list_A[[x]]$gene,
                                     TERM2GENE = bg_genes_A))
names(res_A) <- names(meth_results_list_A)
#Convert the list of enrichResults for each sample_pattern to a dataframe with the pathways
res_df_A <- lapply(names(res_A), function(x) rbind(res_A[[x]]@result))
names(res_df_A) <- names(res_A)
res_df_A <- do.call(rbind, res_df_A)
head(res_df_A)

#Convert the enrichResults to a dataframe with the pathways
res_df_A <- lapply(names(res_A), function(x) rbind(res_A[[x]]@result))
names(res_df_A) <- names(res_A)
res_df_A <- do.call(rbind, res_df_A)
head(res_df_A)

res_df_A <- res_df_A %>% mutate(minuslog10padj = -log10(p.adjust),
                                diffexpressed = gsub('\\.REACTOME.*$', '', rownames(res_df_A)))

# Subset to those pathways that have p adj < cutoff and gene count > cutoff (you can also do this in the enricher function)
target_pws <- unique(res_df_A$ID[res_df_A$p.adjust < padj_cutoff & res_df_A$Count > genecount_cutoff]) # select only target pathways have p adjusted < 0.05 and at least 6 genes
res_df_A <- res_df_A[res_df_A$ID %in% target_pws, ]

write.csv(res_df_A, 'meth_resclusterp_Aerobic.csv', row.names = FALSE)

res_df_down_A <- res_df_A %>% filter(grepl("Hypo", diffexpressed)) %>% 
  dplyr::select(!c('minuslog10padj', 'diffexpressed')) 
rownames(res_df_down_A) <- res_df_down_A$ID

# For visualisation purposes, let's shorten the pathway names
res_df_down_A$Description <- gsub('(H|h)iv', 'HIV', 
                                  gsub('pd 1', 'PD-1',
                                       gsub('ecm', 'ECM', 
                                            gsub('(I|i)nterleukin', 'IL', 
                                                 gsub('(R|r)na', 'RNA', 
                                                      gsub('(D|d)na', 'DNA',
                                                           gsub(' i ', ' I ', 
                                                                gsub('(A|a)tp ', 'ATP ', 
                                                                     gsub('(N|n)adh ', 'NADH ', 
                                                                          gsub('(N|n)ad ', 'NAD ',
                                                                               gsub('t cell', 'T cell',
                                                                                    gsub('b cell', 'B cell',
                                                                                         gsub('built from .*', ' (...)',
                                                                                              gsub('mhc', 'MHC',
                                                                                                   gsub('mhc class i', 'MHC I', 
                                                                                                        gsub('mhc class ii', 'MHC II', 
                                                                                                             stringr::str_to_sentence(
                                                                                                               gsub('_', ' ',  
                                                                                                                    gsub('REACTOME_', '', res_df_down_A$Description)))))))))))))))))))

enrichres_down_A <- new("enrichResult",
                        readable = FALSE,
                        result = res_df_down_A,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.2,
                        organism = "human",
                        ontology = "UNKNOWN",
                        gene = meth_fix2$gene,
                        keytype = "UNKNOWN",
                        universe = unique(bg_genes_A$gene),
                        gene2Symbol = character(0),
                        geneSets = bg_genes_A)
class(enrichres_down)


# Dotplot
tiff('dot plot hypo pathways.tiff',
     width =200,
     height = 300,
     units = 'mm',
     res=300)
dotplot(enrichres_down_A, showCategory = 20) + ggtitle("Top pathways from hypo methylated DMRs")
dev.off()

tiff('dot plot hyper pathways.tiff',
     width =200,
     height = 300,
     units = 'mm',
     res=300)
dotplot(enrichres_up_A, showCategory = 20) + ggtitle("Top pathways from hyper methylated DMRs")
dev.off()


####ven diagram####

DNAm_A <- meth_aerobic %>%
  filter(FDR<0.05)%>%
  dplyr::rename(Gene = `Annotated.gene.s.`)%>%
  separate_rows(Gene,
                sep=";")%>%
  dplyr::select(Gene,
                `Effect.size`)

DNAm_R <- meth_resistance %>%
  filter(FDR<0.05)%>%
  dplyr::rename(Gene = `Annotated.gene.s.`)%>%
  separate_rows(Gene,
                sep=";")%>%
  dplyr::select(Gene,
                `Effect.size`)


merged <- left_join(DNAm_R,
                   DNAm_A)%>%
  drop_na()

merged=drop_na(merged)

#Venn diagram
library("ggVennDiagram")
library("ggplot2")
library("sf")
library("ggvenn")

x <- list("DMGs Aerobic" = genesZscore,
          "DMGs Resistance" = DNAm_R$Gene)
value=c(DNAm_A$Gene,DNAm_R$Gene)
value=unique(value)

data_venn <- data.frame(value = value,    # Create example data frame
                        "DMGs Aerobic" = FALSE,
                        "DMGs Resistance" = FALSE)
data_venn$`DMGs.Resistance` <- data_venn$value %in% x$`DMGs Resistance`
data_venn$`DMGs.Aerobic` <- data_venn$value %in% x$`DMGs Aerobic`

ggplot(data_venn) +
  geom_venn(aes(A=`DMGs.Aerobic`, B=`DMGs.Resistance`), fill_color = c("#1AEDD8", "#1BD2F7"), 
            stroke_size = 1, set_name_size = 5, fill_alpha = 0.3, text_size = 5, 
            show_percentage = FALSE, text_color = "#333333")+
  theme_classic()+
  theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank())

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise/Methylation/Tables & Figures")
tiff('Venn diagram.tiff',
     width =100,
     height = 100,
     units = 'mm',
     res=300)
ggplot(data_venn, aes(A=DMGs, B=DEGs, C=DEPs)) +
  geom_venn(fill_color = c("#1AEDD8", "#1BD2F7", "#CC99FF"), 
            stroke_size = 1, set_name_size = 5, fill_alpha = 0.3, text_size = 5, 
            show_percentage = FALSE, text_color = "#333333")+
  theme_classic()+
  theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank())

dev.off()

####Use Mitch to contrast results####
#change table for mitch
library("mitch")
library(tidyverse)

#add genes to methylation results
meth_A1 <- meth_aerobic%>%
  mutate(t=Effect.size/SE)%>%
  select("CpG","Effect.size","SE","P.value","FDR", "t","Annotated.gene.s.")
dim(meth_A1)
methg_A <- meth_A1$`Annotated.gene.s.`
methg_spl_A <- strsplit(methg_A,";")
meth_A1 <- meth_A1[which(lapply(methg_spl_A,length)>0),]
head(meth_A1)
dim(meth_A1)
x <- apply(meth_A1,1,function(x) {
  GENES <- as.character(x[length(x)])
  GENES <- unlist(strsplit(GENES,";"))
  z <- lapply(GENES,function(G) {
    y <- x[1:(length(x)-1)] 
    y$Gene <- G
    y
  })
  do.call(rbind,z)
})
meth_fix_A <- do.call(rbind,x)
meth_fix_A2 <- as.data.frame(apply(meth_fix_A[,c(2:6)],2,as.numeric))
meth_fix_A2$gene <- unlist(meth_fix_A[,"Gene"])
head(meth_fix_A2)
dim(meth_fix_A2)

#add genes to methylation results - resistance
meth_R1 <- meth_resistance%>%
  mutate(t=Effect.size/SE)%>%
  select("CpG","Effect.size","SE","P.value","FDR", "t","Annotated.gene.s.")
dim(meth_R1)
methg_R <- meth_R1$`Annotated.gene.s.`
methg_spl_R <- strsplit(methg_R,";")
meth_R1 <- meth_R1[which(lapply(methg_spl_R,length)>0),]
head(meth_R1)
dim(meth_R1)
x <- apply(meth_R1,1,function(x) {
  GENES <- as.character(x[length(x)])
  GENES <- unlist(strsplit(GENES,";"))
  z <- lapply(GENES,function(G) {
    y <- x[1:(length(x)-1)] 
    y$Gene <- G
    y
  })
  do.call(rbind,z)
})
meth_fix_R <- do.call(rbind,x)
meth_fix_R2 <- as.data.frame(apply(meth_fix_R[,c(2:6)],2,as.numeric))
meth_fix_R2$gene <- unlist(meth_fix_R[,"Gene"])
head(meth_fix_R2)
dim(meth_fix_R2)


#Now import with mitch
meth_A1 <- meth_fix_A2[,c("gene","t")]
meth_A1 <- aggregate(. ~ gene, meth_A1, sum)
dim(meth_A1)

meth_R1 <- meth_fix_R2[,c("gene","t")]
meth_R1 <- aggregate(. ~ gene, meth_R1, sum)
dim(meth_R1)

l <- list("meth aerobic"=meth_A1, "meth resistance"=meth_R1)
m <- mitch_import(x = l,geneIDcol = "gene",DEtype = "limma")

#Fetch gene sets 
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip", destfile="ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
genesets <- gmt_import("ReactomePathways.gmt")

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/MSigDB files/MSigDB/msigdb_v2022.1.Hs_GMTs")
genesets2<- gmt_import("c5.all.v2023.2.Hs.symbols.gmt")

#Calculate enrichment priority as significance
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise")

res <- mitch_calc(x=m,genesets = genesets2 ,priority = "significance",minsetsize = 5, resrows = 25)
head(res$enrichment_result,50)
mitch_report(res = res, outfile = "myreport_significance_hpo_Aerobic vs res.html")
mitch_plots(res,outfile="mycharts_significance_hpo_Aerobic vs res.pdf")

#save results as table
res=res$enrichment_result
res_sig=res[res$p.adjustMANOVA<0.05 ,]

#Reproduce plot using results table
res1 = mutate(res, sig=ifelse(res$p.adjustMANOVA<0.05, "FDR<0.05", "Not Sig"))

library("ggrepel")

tiff('Effect size versus significance pathways Aerobic vs Resistance hpo.tiff', width =8, height = 6, units = 'in', res=600)
pvsef=ggplot(res1, aes(x=`s.dist`, y=-log(`p.adjustMANOVA`)))+
  geom_point(aes(col=sig))+
  scale_color_manual(values=c("#1BD2F7", "gray"), name = " ")+
  theme_minimal()+
  labs(y="-log(p.adjustedMANOVA)", x="s.dist")
pvsef
dev.off()

#Calculate enrichment priority effect size
res <- mitch_calc(x=m,genesets = genesets ,priority = "effect",minsetsize = 5, resrows = 25)
head(res$enrichment_result,50)
mitch_report(res = res, outfile = "myreport_effect_hpo.html")
mitch_plots(res,outfile="mycharts_effect_hpo.pdf")

#Reproduce HeatMap
d=ncol(res$input_profile)

if (d==1) {
  formatted<-t(as.data.frame(res$analysis_metrics[ c(1,2,3,4,5 ) ]))
} else if (d==2) {
  unformatted<-t(as.data.frame(res$analysis_metrics[ c(2,3,4,5,11,12 )  ]))
  formatted<-unformatted
  formatted[1:4]<-as.character(round(as.numeric(unformatted[1:4]) , digits=0))
  formatted[5:6]<-as.character(round(as.numeric(unformatted[5:6]) , digits=5))
} else if (d>2) {
  formatted<-t(as.data.frame(res$analysis_metrics[ c(2,3,4,5 )  ]))
}

colnames(formatted)="Profile metrics"

formatted %>%
  kbl(caption="Profiling data metrics") %>%
  kable_styling("hover", full_width = FALSE)



hmapx<-head(res$enrichment_result[,4:(4+d-1)],14)
rownames(hmapx)<-head(res$enrichment_result$set,14)
colnames(hmapx)<-gsub("^s.","",colnames(hmapx))
my_palette <- colorRampPalette(c("#CC99FF", "#1BD2F7","white","#1AEDD8"))(n = 10)
hmapx=as.matrix(hmapx)
colnames(hmapx)=c("DNAm Aerobic", "mRNA", "Prot")


tiff('Heatplot integration.tiff',
     width =200,
     height = 250,
     units = 'mm',
     res=600)
heatmap.2(hmapx,scale="none",Rowv = T, Colv = FALSE ,margin=c(10, 25),
          cexRow=0.8,trace="none",cexCol=1,col=my_palette,
          dendrogram="none")
dev.off()


####DMRs for Aerobic####
library(DMRcate)
library(IRanges)

meth_aerobic_sig=meth_aerobic%>%
  filter(FDR<0.05)

CpGs <-meth_aerobic_sig$CpG # the ones with a p-value below 0.05
Unique_IDs <- meth_aerobic_sig$`Unique.ID`


library(tidyverse)
annotation_overlap_only=filter(annotation, Unique_ID %in% Unique_IDs)
annotation_overlap_only=annotation_overlap_only[!duplicated(annotation_overlap_only$Unique_ID),]


annotation_overlap_only=left_join(annotation_overlap_only, meth_aerobic_sig, by=join_by(Unique_ID==`Unique.ID`))
annotation_overlap_only=mutate(annotation_overlap_only, t=annotation_overlap_only$`Effect.size`/annotation_overlap_only$SE)

rownames(annotation_overlap_only)=annotation_overlap_only$Unique_ID

annotated <- GenomicRanges::GRanges(annotation_overlap_only[Unique_IDs,"CpG_chrm" ], #chromosome
                                    IRanges(annotation_overlap_only[Unique_IDs,"CpG_beg"],annotation_overlap_only[Unique_IDs,"CpG_end"]), #add location #put the same location twice so it recognises as a region
                                    stat = annotation_overlap_only[Unique_IDs,"t"], #t-statistic
                                    diff = annotation_overlap_only[Unique_IDs,"Effect.size"], #effect size from beta 
                                    ind.fdr = annotation_overlap_only[Unique_IDs,"FDR"], #adjusted p-value
                                    is.sig = annotation_overlap_only[Unique_IDs,"FDR"] < 0.05) #p-value threshold #logical operator TRUE of FALSE so it recognises later on

names(annotated) <- annotation_overlap_only$Unique_ID #name each line with the CpG
annotated <- sort(annotated) #sorting the CpGs by position on the DNA
annotated <- new("CpGannotated", ranges = annotated) #create a "CpGannotated" object

DMR <- dmrcate(annotated, lambda=1000, C=2, min.cpgs = 2) #C is an statistical factor, optimal at 2 #minimal number of CpGs that wou would consider to have a DMR default=2
DMR

#obtain results
resultsRanges <- extractRanges(DMR, genome = "hg38") #different gene annotations this is the most updated one
resultsRanges_dataframe=as.data.frame (resultsRanges)

#Obtain better annotation of the DMRs
chromHMM <- read.delim("E:/Macsue - Hard drive/PhD Victoria University/R Analyses/R tutorial - Sarah/Methylation Analyses/Chromatin states annotation.txt")
rownames(chromHMM) <- paste(chromHMM$`STATE.NO.`,chromHMM$MNEMONIC,sep="_")
#Replace CpG island names with better words
annotation$CGIposition[grep("Shore",annotation$CGIposition)]="Shore"
annotation$CGIposition[grep("Shelf",annotation$CGIposition)]="Shelf"
annotation$CGIposition[is.na(annotation$CGIposition)]="Open sea"
annotation_GR <- makeGRangesFromDataFrame(annotation,
                                          keep.extra.columns=TRUE,
                                          ignore.strand=FALSE,
                                          seqinfo=NULL,
                                          seqnames.field=c("CpG_chrm"),
                                          start.field="CpG_beg",
                                          end.field=c("CpG_end"),
                                          strand.field="probe_strand",
                                          starts.in.df.are.0based=FALSE)
genesidx <- as.data.frame(findOverlaps(resultsRanges, annotation_GR))
genesover <- tapply(genesidx$subjectHits, genesidx$queryHits,
                    function(x) annotation_GR$genesUniq_with_enh[x])
cpgislands <- tapply(genesidx$subjectHits, genesidx$queryHits,
                     function(x) annotation_GR$CGIposition[x])
chromstatemale <- tapply(genesidx$subjectHits, genesidx$queryHits,
                         function(x) annotation_GR$E107[x])
chromstatefemale <- tapply(genesidx$subjectHits, genesidx$queryHits,
                           function(x) annotation_GR$E108[x])

op.A <- sapply(genesover, function(l) paste(unique(unlist(strsplit(l,split=";"))), collapse = ";"))
name.A <- names(genesover)
m.A <- as.numeric(name.A)
M <- length(resultsRanges)
better_annotation <- rep("", M)
better_annotation[m.A] <- op.A
op.A <- sapply(cpgislands, function(l) paste(unique(l), collapse = ";"))
name.A <- names(cpgislands)
m.A <- as.numeric(name.A)
M <- length(resultsRanges)
cpg_island <- rep("", M)
cpg_island[m.A] <- op.A
op.A.male <- sapply(chromstatemale, function(l) paste(unique(unlist(strsplit(chromHMM[l,"DESCRIPTION"],split=";"))), collapse = ";"))
name.A <- names(chromstatemale)
m.A <- as.numeric(name.A)
M <- length(resultsRanges)
chromstate_male <- rep("", M)
chromstate_male[m.A] <- op.A.male

op.A.female <- sapply(chromstatefemale, function(l) paste(unique(unlist(strsplit(chromHMM[l,"DESCRIPTION"],split=";"))), collapse = ";"))
name.A <- names(chromstatefemale)
m.A <- as.numeric(name.A)
M <- length(resultsRanges)
chromstate_female <- rep("", M)
chromstate_female[m.A] <- op.A.female
#Change to tibble and format
results.ranges <- as_tibble(resultsRanges)
results.ranges$better_annotation <- better_annotation
results.ranges$cpg_island <- cpg_island
results.ranges$chromstate_male <- chromstate_male
results.ranges$chromstate_female <- chromstate_female
DMRs <- results.ranges%>%
  dplyr::filter(Stouffer<0.05&
                  HMFDR <0.05&
                  Fisher<0.05) %>%
  dplyr::select(seqnames,
                start,
                end,
                width,
                no.cpgs,
                cpg_island,
                chromstate_male,
                chromstate_female,
                better_annotation,
                maxdiff,
                meandiff,
                Stouffer,
                HMFDR,
                Fisher)
colnames(DMRs) <- c("Chromosome",
                    "Position start (hg38)",
                    "Position end (hg38)",
                    "Length of DMR (base pairs)",
                    "Number of CpGs in DMR",
                    "CpG island position",
                    "Chromatin state in male skeletal muscle",
                    "Chromatin state in female skeletal muscle",
                    "Annotated gene(s)",
                    "Maximum effect size in DMR",
                    "Mean effect size in DMR",
                    "Stouffer score",
                    "Harmonic mean of the individual component FDRs",
                    "Fisher multiple comparison statistic")
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise/Methylation/Tables & Figures")
write.csv(DMRs,
          file="DMRs with Aerobic training.csv")
DMRs_A=read.csv("DMRs with Aerobic training.csv")

library(tidyverse)

top_DMRs_A=DMRs_A%>%
  dplyr::filter(Stouffer.score<0.05&
                  Harmonic.mean.of.the.individual.component.FDRs <0.05&
                  Fisher.multiple.comparison.statistic<0.05)

#Find how many unique genes in DMRs
genes=top_DMRs_A$Annotated.gene.s.
genes=c(str_split(genes, ";"))
genes=unlist(genes)
genes=c(str_split(genes, ","))
genes2=unlist(genes)
genesZscore=unique(genes2)
view(genesZscore) #8752 unique genes in the DMRs

#Calculate % of hypo and hyper DMRs
sum(DMRs_A$Mean.effect.size.in.DMR >0)*100/nrow(DMRs_A) #58.76
sum(DMRs_A$Mean.effect.size.in.DMR <0)*100/nrow(DMRs_A) #41.23


####DMRs for Resistance####
library(DMRcate)
library(IRanges)

meth_resistance_sig=meth_resistance%>%
  filter(FDR<0.05)

CpGs <-meth_resistance_sig$CpG # the ones with a p-value below 0.05
Unique_IDs <- meth_resistance_sig$`Unique.ID`


library(tidyverse)
annotation_overlap_only=filter(annotation, Unique_ID %in% Unique_IDs)
annotation_overlap_only=annotation_overlap_only[!duplicated(annotation_overlap_only$Unique_ID),]


annotation_overlap_only=left_join(annotation_overlap_only, meth_resistance_sig, by=join_by(Unique_ID==`Unique.ID`))
annotation_overlap_only=mutate(annotation_overlap_only, t=annotation_overlap_only$`Effect.size`/annotation_overlap_only$SE)

rownames(annotation_overlap_only)=annotation_overlap_only$Unique_ID

annotated <- GenomicRanges::GRanges(annotation_overlap_only[Unique_IDs,"CpG_chrm" ], #chromosome
                                    IRanges(annotation_overlap_only[Unique_IDs,"CpG_beg"],annotation_overlap_only[Unique_IDs,"CpG_end"]), #add location #put the same location twice so it recognises as a region
                                    stat = annotation_overlap_only[Unique_IDs,"t"], #t-statistic
                                    diff = annotation_overlap_only[Unique_IDs,"Effect.size"], #effect size from beta 
                                    ind.fdr = annotation_overlap_only[Unique_IDs,"FDR"], #adjusted p-value
                                    is.sig = annotation_overlap_only[Unique_IDs,"FDR"] < 0.05) #p-value threshold #logical operator TRUE of FALSE so it recognises later on

names(annotated) <- annotation_overlap_only$Unique_ID #name each line with the CpG
annotated <- sort(annotated) #sorting the CpGs by position on the DNA
annotated <- new("CpGannotated", ranges = annotated) #create a "CpGannotated" object

DMR <- dmrcate(annotated, lambda=1000, C=2, min.cpgs = 2) #C is an statistical factor, optimal at 2 #minimal number of CpGs that wou would consider to have a DMR default=2
DMR

#obtain results
resultsRanges <- extractRanges(DMR, genome = "hg38") #different gene annotations this is the most updated one
resultsRanges_dataframe=as.data.frame (resultsRanges)

#Obtain better annotation of the DMRs
chromHMM <- read.delim("E:/Macsue - Hard drive/PhD Victoria University/R Analyses/R tutorial - Sarah/Methylation Analyses/Chromatin states annotation.txt")
rownames(chromHMM) <- paste(chromHMM$`STATE.NO.`,chromHMM$MNEMONIC,sep="_")
#Replace CpG island names with better words
annotation$CGIposition[grep("Shore",annotation$CGIposition)]="Shore"
annotation$CGIposition[grep("Shelf",annotation$CGIposition)]="Shelf"
annotation$CGIposition[is.na(annotation$CGIposition)]="Open sea"
annotation_GR <- makeGRangesFromDataFrame(annotation,
                                          keep.extra.columns=TRUE,
                                          ignore.strand=FALSE,
                                          seqinfo=NULL,
                                          seqnames.field=c("CpG_chrm"),
                                          start.field="CpG_beg",
                                          end.field=c("CpG_end"),
                                          strand.field="probe_strand",
                                          starts.in.df.are.0based=FALSE)
genesidx <- as.data.frame(findOverlaps(resultsRanges, annotation_GR))
genesover <- tapply(genesidx$subjectHits, genesidx$queryHits,
                    function(x) annotation_GR$genesUniq_with_enh[x])
cpgislands <- tapply(genesidx$subjectHits, genesidx$queryHits,
                     function(x) annotation_GR$CGIposition[x])
chromstatemale <- tapply(genesidx$subjectHits, genesidx$queryHits,
                         function(x) annotation_GR$E107[x])
chromstatefemale <- tapply(genesidx$subjectHits, genesidx$queryHits,
                           function(x) annotation_GR$E108[x])

op.A <- sapply(genesover, function(l) paste(unique(unlist(strsplit(l,split=";"))), collapse = ";"))
name.A <- names(genesover)
m.A <- as.numeric(name.A)
M <- length(resultsRanges)
better_annotation <- rep("", M)
better_annotation[m.A] <- op.A
op.A <- sapply(cpgislands, function(l) paste(unique(l), collapse = ";"))
name.A <- names(cpgislands)
m.A <- as.numeric(name.A)
M <- length(resultsRanges)
cpg_island <- rep("", M)
cpg_island[m.A] <- op.A
op.A.male <- sapply(chromstatemale, function(l) paste(unique(unlist(strsplit(chromHMM[l,"DESCRIPTION"],split=";"))), collapse = ";"))
name.A <- names(chromstatemale)
m.A <- as.numeric(name.A)
M <- length(resultsRanges)
chromstate_male <- rep("", M)
chromstate_male[m.A] <- op.A.male

op.A.female <- sapply(chromstatefemale, function(l) paste(unique(unlist(strsplit(chromHMM[l,"DESCRIPTION"],split=";"))), collapse = ";"))
name.A <- names(chromstatefemale)
m.A <- as.numeric(name.A)
M <- length(resultsRanges)
chromstate_female <- rep("", M)
chromstate_female[m.A] <- op.A.female
#Change to tibble and format
results.ranges <- as_tibble(resultsRanges)
results.ranges$better_annotation <- better_annotation
results.ranges$cpg_island <- cpg_island
results.ranges$chromstate_male <- chromstate_male
results.ranges$chromstate_female <- chromstate_female
DMRs <- results.ranges%>%
  dplyr::filter(Stouffer<0.05&
                  HMFDR <0.05&
                  Fisher<0.05) %>%
  dplyr::select(seqnames,
                start,
                end,
                width,
                no.cpgs,
                cpg_island,
                chromstate_male,
                chromstate_female,
                better_annotation,
                maxdiff,
                meandiff,
                Stouffer,
                HMFDR,
                Fisher)
colnames(DMRs) <- c("Chromosome",
                    "Position start (hg38)",
                    "Position end (hg38)",
                    "Length of DMR (base pairs)",
                    "Number of CpGs in DMR",
                    "CpG island position",
                    "Chromatin state in male skeletal muscle",
                    "Chromatin state in female skeletal muscle",
                    "Annotated gene(s)",
                    "Maximum effect size in DMR",
                    "Mean effect size in DMR",
                    "Stouffer score",
                    "Harmonic mean of the individual component FDRs",
                    "Fisher multiple comparison statistic")
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise/Methylation/Tables & Figures")
write.csv(DMRs,
          file="DMRs with Aerobic training.csv")
DMRs_A=read.csv("DMRs with Aerobic training.csv")

library(tidyverse)

top_DMRs_A=DMRs_A%>%
  dplyr::filter(Stouffer.score<0.05&
                  Harmonic.mean.of.the.individual.component.FDRs <0.05&
                  Fisher.multiple.comparison.statistic<0.05)

#Find how many unique genes in DMRs
genes=meth_resistance_sig$Annotated.gene.s.
genes=c(str_split(genes, ";"))
genes=unlist(genes)
genes=c(str_split(genes, ","))
genes2=unlist(genes)
genesZscore_R=unique(genes2)
view(genesZscore) #8752 unique genes in the DMRs Aerobic
view(genesZscore_R) #157 unique genes in the DMRs Aerobic

#Calculate % of hypo and hyper DMRs
sum(DMRs_A$Mean.effect.size.in.DMR >0)*100/nrow(DMRs_A) #58.76
sum(DMRs_A$Mean.effect.size.in.DMR <0)*100/nrow(DMRs_A) #41.23




####Transcriptomics meta-analyses Extrameta####
library(readxl)
setwd("E:/Macsue - Hard drive/PhD Victoria University/Papers to work on/Methylation Meta-Analyses/Transcriptomics")
Extrameta=read.delim("Extrameta all genes results.txt")

Extrameta$FDR=p.adjust(Extrameta$p)

Significant_Extrameta=Extrameta[Extrameta$FDR<0.05,]
Significant_mRNA=unique(Significant_Extrameta$gene)

Extrameta=Extrameta%>%
  mutate(t=ES/SE)

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise")
write.csv(Significant_Extrameta, "meta-analyses extrameta significant.csv")

#Calculate % of up down mRNAs
sum(Significant_Extrameta$ES >0)*100/nrow(Significant_Extrameta)
sum(Significant_Extrameta$ES <0)*100/nrow(Significant_Extrameta)

####PLOT 1: Volcano plot
library(ggplot2)
library(ggpubr)

#Add in different color the points that are significant
library(dplyr)
mRNAs=Extrameta$gene

#Volcano Plot for p-values 0.005 threshold
#add a column of NAs
all_results = mutate(Extrameta, diff=">0.05")

# if Estimate > 0 and pvalue < 0.05, set as "UP" 
all_results$diff[all_results$ES > 0.0 & all_results$FDR <0.05] <- "increase"
# if log2Foldchange < -0.0 and pvalue < 0.1, set as "DOWN"
all_results$diff[all_results$ES < -0.0 & all_results$FDR <0.05] <- "decrease"

all_results=data.frame(Extrameta,all_results)

all_results=dplyr::arrange(all_results, p)
all_results=all_results[47:14309,]
gene=all_results$gene

tiff('Volcano_signif mRNA after exercise.tiff', width =7, height = 5, units = 'in', res=600)
ggplot() +
  geom_point(data = all_results,
             aes(ES,
                 -log10(`p`),
                 col=diff), size=0.85)+
  xlim(-1,1)+
  scale_color_manual(values=c("grey","#1BD2F7", "#1AEDD8"))+
  labs(y="-log10(p-value)")+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

####Enrichment with cluster profiler####
# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(pheatmap)
library(clusterProfiler) # for PEA analysis
library('org.Hs.eg.db')
library(DOSE)
library(enrichplot) # for visualisations
library(ggupset) # for visualisations

## Function: Adjacency matrix to list ####
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

# Get the genes that are present in your dataframe
genes_in_transcriptome <- all_results$gene
# Read in the .gmt file
pwl2 <- read.gmt("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/ReactomePathways.gmt") 
# Subset to the genes that are present in our dataset
pwl2_R <- pwl2[pwl2$gene %in% genes_in_transcriptome,] 
# Save the filtered background gene set
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise")
filename <- 'transcriptome_path_R.RDS'
saveRDS(pwl2_R, filename)

# Remove non-significant genes
trans_path <- all_results[all_results$diff != '>0.05', ]
# Substitute names so they are annotated nicely in the heatmap later
# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
trans_results_list <- split(trans_path, trans_path$diff)

## Run ClusterProfiler -----------------------------------------------
# Settings
name_of_comparison <- 'Transcriptome' # for our filename
background_genes <- 'reactome' # for our filename
bg_genes_tras <- readRDS("transcriptome_path_R.RDS") # read in the background genes
padj_cutoff <- 0.05 # p-adjusted threshold, used to filter out pathways
genecount_cutoff <- 5 # minimum number of genes in the pathway, used to filter out pathways

res_trans <- lapply(names(trans_results_list),
                function(x) enricher(gene = trans_results_list[[x]]$gene,
                                     TERM2GENE = bg_genes_tras))
names(res_trans) <- names(trans_results_list)
#Convert the list of enrichResults for each sample_pattern to a dataframe with the pathways
res_df_trans <- lapply(names(res_trans), function(x) rbind(res_trans[[x]]@result))
names(res_df_trans) <- names(res_trans)
res_df_trans <- do.call(rbind, res_df_trans)
head(res_df_trans)

res_df_trans <- res_df_trans %>% mutate(minuslog10padj = -log10(p.adjust),
                                diffexpressed = gsub('\\.REACTOME.*$', '', rownames(res_df_trans)))

# Subset to those pathways that have p adj < cutoff and gene count > cutoff (you can also do this in the enricher function)
target_pws <- unique(res_df_trans$ID[res_df_trans$p.adjust < padj_cutoff ]) # select only target pathways have p adjusted < 0.05 and at least 6 genes
res_df_trans <- res_df_trans[res_df_trans$ID %in% target_pws, ]

write.csv(res_df_trans, 'transcriptome_resclusterp.csv', row.names = FALSE)

res_df_up_trans <- res_df_trans %>% filter(grepl("increase", diffexpressed)) %>% 
  dplyr::select(!c('minuslog10padj', 'diffexpressed')) 
rownames(res_df_up_trans) <- res_df_up_trans$ID

# For visualisation purposes, let's shorten the pathway names
res_df_up_trans$Description <- gsub('(H|h)iv', 'HIV', 
                                  gsub('pd 1', 'PD-1',
                                       gsub('ecm', 'ECM', 
                                            gsub('(I|i)nterleukin', 'IL', 
                                                 gsub('(R|r)na', 'RNA', 
                                                      gsub('(D|d)na', 'DNA',
                                                           gsub(' i ', ' I ', 
                                                                gsub('(A|a)tp ', 'ATP ', 
                                                                     gsub('(N|n)adh ', 'NADH ', 
                                                                          gsub('(N|n)ad ', 'NAD ',
                                                                               gsub('t cell', 'T cell',
                                                                                    gsub('b cell', 'B cell',
                                                                                         gsub('built from .*', ' (...)',
                                                                                              gsub('mhc', 'MHC',
                                                                                                   gsub('mhc class i', 'MHC I', 
                                                                                                        gsub('mhc class ii', 'MHC II', 
                                                                                                             stringr::str_to_sentence(
                                                                                                               gsub('_', ' ',  
                                                                                                                    gsub('REACTOME_', '', res_df_up_trans$Description)))))))))))))))))))

enrichres_up_trans <- new("enrichResult",
                        readable = FALSE,
                        result = res_df_up_trans,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.2,
                        organism = "human",
                        ontology = "UNKNOWN",
                        gene = all_results$gene,
                        keytype = "UNKNOWN",
                        universe = unique(bg_genes_tras$gene),
                        gene2Symbol = character(0),
                        geneSets = bg_genes_tras)
class(enrichres_up_trans)


# Dotplot
tiff('dot plot down pathways transcriptome.tiff',
     width =200,
     height = 300,
     units = 'mm',
     res=300)
dotplot(enrichres_down_trans, showCategory = 20) + ggtitle("Top pathways from down expressed DEGs")
dev.off()

tiff('dot plot up pathways transcriptome.tiff',
     width =200,
     height = 300,
     units = 'mm',
     res=300)
dotplot(enrichres_up_trans, showCategory = 20) + ggtitle("Top pathways from up expressed DEGs")
dev.off()



####Meta-Analysis Proteomics for Exercise####
setwd("E:/Macsue - Hard drive/PhD Victoria University/PhD related/Proteomics project/Pre-processing and analyses for all samples")
Proteomics=read.csv("MetaData_norm_and_batch.csv")
setwd("E:/Macsue - Hard drive/PhD Victoria University/Papers to work on/Methylation Meta-Analyses/New analyses - VO2max and aerobic only")
combined_pheno=read.csv("combined_pheno.csv")
rownames(combined_pheno)=combined_pheno$X

#Change proteomics data so gene names are rows and columns IDs to fit limma
rownames(Proteomics)=Proteomics$Gene.names
Proteomics=Proteomics[,33:187]
Proteomics=Proteomics[,-c(12,16,29,45,60,75,89,104,111,120,133,134,148)]#Remove technical replicates
Proteomics=Proteomics%>%
  dplyr::rename("SG156_PRE"="SG166_PRE", "SG162_PRE"="SG162_PRE.1","SG166_PRE"="SG166_PRE.1")

Pheno_Proteomics=combined_pheno[colnames(Proteomics),]
colnames(Proteomics)==rownames(Pheno_Proteomics)

Pheno_Proteomics=Pheno_Proteomics%>%
  mutate(Timepoint=fct_relevel(Timepoint, "PRE", "4WP", "8WP","12WP"))

#run limma
library(limma)

#change time for continuous
design=model.matrix(~VO2max_rel+Timepoint+Age+Type_I+intervention+sex,data=Pheno_Proteomics) #time as continuous

corfit_PRE_Prot <- duplicateCorrelation(Proteomics,design,block=Pheno_Proteomics$ID)
#Then, we need to input this inter-correlation into the linear model to take into account the paired design
fit_PRE_Prot = lmFit(Proteomics, design,maxit=1000,block=Pheno_Proteomics$ID,correlation=corfit_PRE_Prot$consensus) #I have chosen the option "robust" which will downplay the influence of outliers. This takes significantly more time for your computer to run but is more robust and is always my preference.

#get info on eBayes in the vignette of the limma package
fit2_PRE_Prot <- eBayes(fit_PRE_Prot)

signif_VO2max_PRE_prot=topTable(fit2_PRE_Prot, coef="VO2max_rel", adjust="BH",number=Inf,  p.value = 0.005)
signif_timepoint_PRE_prot=limma::topTable(fit2_PRE_Prot, coef="Timepoint4WP", adjust="BH",number=Inf,  p.value = 0.05)
signif_intervention_PRE_prot=topTable(fit2_PRE_Prot, coef="intervention", adjust="BH",number=Inf,  p.value = 0.005)

Prot_VO2max_PRE=topTable(fit2_PRE_Prot, coef="VO2max_rel", adjust="BH",number=Inf,  p.value = 1)
Prot_timepoint_PRE=limma::topTable(fit2_PRE_Prot, coef="Timepoint4WP", adjust="BH",number=Inf,  p.value = 1)
Prot_intervention_PRE=topTable(fit2_PRE_Prot, coef="intervention", adjust="BH",number=Inf,  p.value = 1)

SE <- fit2_PRE_Prot$sigma * fit2_PRE_Prot$stdev.unscaled
SE <- SE[rownames(Proteomics),"Timepoint4WP"]
Proteomics_timepoint_PRE <- cbind(Prot_timepoint_PRE,SE)

#Save results
Gene <- rownames(Proteomics_timepoint_PRE) #name of marker (CpG)
ALLELE1 <- rep(1,nrow(Proteomics_timepoint_PRE)) #fake allele info (not important)
ALLELE2 <- rep(2,nrow(Proteomics_timepoint_PRE)) #fake allele info (not important)
TESTSTAT <- Proteomics_timepoint_PRE$t #t-statistic
PVALUE <- Proteomics_timepoint_PRE$P.Value #p-value
EFFECTSIZE <- Proteomics_timepoint_PRE$logFC #effect size
SE <- Proteomics_timepoint_PRE$SE #standard error

Proteomics_timepoint_PRE = data.frame(Gene,
                                      ALLELE1,
                                      ALLELE2,
                                      TESTSTAT,
                                      PVALUE,
                                      EFFECTSIZE,
                                      SE)

setwd("E:/Macsue - Hard drive/PhD Victoria University/Papers to work on/Methylation Meta-Analyses/New analyses - VO2max and aerobic only/Final Codes")
write.table(Proteomics_timepoint_PRE,
            file="Proteomics_Timepoint.tbl",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

#Load second data and pre-process/analyse
setwd("E:/Macsue - Hard drive/PhD Victoria University/Papers to work on/Methylation Meta-Analyses/New analyses - VO2max and aerobic only/Proteomics_Deshmukh")

library(readxl)

prot_pheno=read_excel("ACE meta for atul - 09.05.20.xlsx", skip=2)
prot_pheno=prot_pheno[1:16, 1:9]
colnames(prot_pheno)=c("ID","Timepiont","Age","Height", "Weight", "HRmax", "VO2max","VO2max_rel", "aPPO (s)")

#Load proteomics data
Data_prot=read_delim("Spectronaut_v14_humanmuscle_ProteinQuant_Report.xls")
Data_prot=Data_prot[,c(3,4,134:181)]
Data_prot=Data_prot[,1:18]
colnames(Data_prot)=c("Accession","Genes", "1_PRE","2_PRE","3_PRE", "4_PRE", "5_PRE", "6_PRE", "7_PRE", "8_PRE",
                      "1_POST","2_POST","3_POST", "4_POST", "5_POST", "6_POST", "7_POST", "8_POST")

library(DEP)
library(SummarizedExperiment)  


#Make a data frame with only label, condition and replicate
prot_pheno2=data.frame(label=colnames(Data_prot[3:18]), condition=prot_pheno$Timepiont, replicate=prot_pheno$ID)

#Make_unique for meake_se
Data_prot[Data_prot=="Filtered"] <- NA
Data_prot2=as.data.frame(sapply(Data_prot[,3:18], as.numeric))
Data_prot2=cbind("Genes"=Data_prot$Genes, "Accession"=Data_prot$Accession, Data_prot2)
data_unique=make_unique(Data_prot2, "Genes", "Accession", delim = ";" )

#make data with only intensity cols
intensity_cols=c(3:18)

prot_pheno2$replicate=as.numeric(prot_pheno2$replicate)

#Summarize experiment using make_se #this log transforms the data
data_se=DEP::make_se(data_unique, intensity_cols, prot_pheno2)

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
plot_pca(data_imp, n=2317, indicate = "condition") +
  labs(title = "Before VSN Normalisation") # n is total number of proteins in the test data.

plot_pca(data_norm, n=2317, indicate = "condition", label = FALSE)+
  labs(title = "After VSN Normalisation")

Prot2=assay(data_norm)
colnames(data_norm)=c("1_PRE","2_PRE","3_PRE", "4_PRE", "5_PRE", "6_PRE", "7_PRE", "8_PRE",
                      "1_POST","2_POST","3_POST", "4_POST", "5_POST", "6_POST", "7_POST", "8_POST")

#run limma
library(limma)

#change time for continuous
design=model.matrix(~VO2max_rel+Timepiont+Age,data=prot_pheno) #time as continuous

corfit_PRE_Prot <- duplicateCorrelation(Prot2,design,block=prot_pheno$ID)
#Then, we need to input this inter-correlation into the linear model to take into account the paired design
fit_PRE_Prot = lmFit(Prot2, design,maxit=1000,block=prot_pheno$ID,correlation=corfit_PRE_Prot$consensus) #I have chosen the option "robust" which will downplay the influence of outliers. This takes significantly more time for your computer to run but is more robust and is always my preference.

#get info on eBayes in the vignette of the limma package
fit2_PRE_Prot <- eBayes(fit_PRE_Prot)

signif_VO2max_PRE_prot=topTable(fit2_PRE_Prot, coef="VO2max_rel", adjust="BH",number=Inf,  p.value = 0.05)
signif_Timepoint_prot=limma::topTable(fit2_PRE_Prot, coef="Timepiont", adjust="BH",number=Inf,  p.value = 0.05)

Prot_VO2max_PRE=topTable(fit2_PRE_Prot, coef="VO2max_rel", adjust="BH",number=Inf,  p.value = 1)
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


setwd("E:/Macsue - Hard drive/PhD Victoria University/Papers to work on/Methylation Meta-Analyses/New analyses - VO2max and aerobic only/Final Codes")
write.table(Timepoint_prot,
            file="Proteomics_Timepoint_prot_atul.tbl",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")  

#Run Bacon

library(bacon)

#List all files ending in "tbl" in the folder
directory=c("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise")

setwd=directory
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
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise")
Meta_Proteomics=read.delim("METAANALYSIS1.tbl")

#Select only genes present on 3 or more studies
Meta_Proteomics_subset=Meta_Proteomics%>%
  filter(HetDf >=0)
Meta_Proteomics_subset=Meta_Proteomics_subset%>%
  mutate(FDR=p.adjust(Meta_Proteomics_subset$P.value, method = "fdr"))%>%
  mutate(t=Effect/StdErr)

signif_meta_Proteomics=Meta_Proteomics_subset[Meta_Proteomics_subset$FDR<0.05,]
signif_meta_Proteomics2=signif_meta_Proteomics[,c(1,4,5,6,12)]

write.csv(signif_meta_Proteomics2, "Metanalyses Proteomics significant results.csv")

sum(signif_meta_Proteomics$Effect >0)*100/nrow(signif_meta_Proteomics)
sum(signif_meta_Proteomics$Effect <0)*100/nrow(signif_meta_Proteomics)

#Build volcano plot with all the results
####PLOT 1: Volcano plot
library(ggplot2)
library(ggpubr)

#Add in different color the points that are significant
library(dplyr)
Gene=Meta_Proteomics_subset$MarkerName
all_results_P = mutate(Meta_Proteomics_subset, sig=ifelse(Meta_Proteomics_subset$FDR <0.05, "FDR<0.05", "Not Sig"))
all_results_P=data.frame(Gene,all_results_P)
all_results_P=all_results[order(all_results_P$FDR),]

all_results_P=mutate(all_results_P, diff=">0.05")
# if Estimate > 0 and pvalue < 0.05, set as "UP" 
all_results$diff[all_results$Effect > 0.0 & all_results$FDR <0.05] <- "Overexpressed"
# if log2Foldchange < -0.0 and pvalue < 0.1, set as "DOWN"
all_results$diff[all_results$Effect < -0.0 & all_results$FDR <0.05] <- "Downexpressed"

library(ggrepel)

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise")

tiff('Volcano_signif Proteomics.tiff', width =7, height = 5, units = 'in', res=600)
ggplot() +
  geom_point(data = all_results,
             aes(Effect,
                 -log10(`P.value`),
                 col=diff), size=0.85)+
  scale_color_manual(values=c("grey","#1BD2F7", "#1AEDD8"))+
  labs(y="-log10(p-value)")+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

####Transcriptomics innactivity####
setwd("E:/Macsue - Hard drive/PhD Victoria University/Papers to work on/Methylation Meta-Analyses/New analyses - VO2max and aerobic only/Inactivity - LITER study")

transcrip_innactivity=read_excel("Differentially expressed genes with disuse from 2021.xlsx", col_names=TRUE,
                                 col_types=c("numeric","text", "numeric","numeric","numeric","numeric","numeric","numeric","numeric"))

Significant_transcrip_innactivity=transcrip_innactivity[transcrip_innactivity$`Stouffer P-value: FDR` <0.05,]
signif_trans_inact=unique(Significant_transcrip_innactivity$`Gene symbol`)

transcrip_innactivity2=transcrip_innactivity%>%
  mutate(t=`LFC: mean`/`LFC: SE`)%>%
  drop_na()

sum(Significant_transcrip_innactivity$`LFC: mean`>0)*100/nrow(Significant_transcrip_innactivity)
sum(Significant_transcrip_innactivity$`LFC: mean`<0)*100/nrow(Significant_transcrip_innactivity)

####Use Mitch to contrast results####
library("mitch")
library(tidyverse)



proteome <- Prot_timepoint_PRE
proteome$gene <- rownames(proteome)
head(proteome)
str(proteome)

transcriptome <- Extrameta
colnames(transcriptome)[6] <- "gene"
head(transcriptome)
str(transcriptome)

transcriptome_inactivity <- transcrip_innactivity2
colnames(transcriptome_inactivity)[2] <- "gene"
head(transcriptome_inactivity)
str(transcriptome_inactivity)

#methylation results
meth_A1 <- meth_aerobic%>%
  mutate(t=Effect.size/SE)%>%
  select("CpG","Effect.size","SE","P.value","FDR", "t","Annotated.gene.s.")
dim(meth_A1)
methg_A <- meth_A1$`Annotated.gene.s.`
methg_spl_A <- strsplit(methg_A,";")
meth_A1 <- meth_A1[which(lapply(methg_spl_A,length)>0),]
head(meth_A1)
dim(meth_A1)
x <- apply(meth_A1,1,function(x) {
  GENES <- as.character(x[length(x)])
  GENES <- unlist(strsplit(GENES,";"))
  z <- lapply(GENES,function(G) {
    y <- x[1:(length(x)-1)] 
    y$Gene <- G
    y
  })
  do.call(rbind,z)
})
meth_fix_A <- do.call(rbind,x)
meth_fix_A2 <- as.data.frame(apply(meth_fix_A[,c(2:6)],2,as.numeric))
meth_fix_A2$gene <- unlist(meth_fix_A[,"Gene"])
head(meth_fix_A2)
dim(meth_fix_A2)

#add genes to methylation results - resistance
meth_R1 <- meth_resistance%>%
  mutate(t=Effect.size/SE)%>%
  select("CpG","Effect.size","SE","P.value","FDR", "t","Annotated.gene.s.")
dim(meth_R1)
methg_R <- meth_R1$`Annotated.gene.s.`
methg_spl_R <- strsplit(methg_R,";")
meth_R1 <- meth_R1[which(lapply(methg_spl_R,length)>0),]
head(meth_R1)
dim(meth_R1)
x <- apply(meth_R1,1,function(x) {
  GENES <- as.character(x[length(x)])
  GENES <- unlist(strsplit(GENES,";"))
  z <- lapply(GENES,function(G) {
    y <- x[1:(length(x)-1)] 
    y$Gene <- G
    y
  })
  do.call(rbind,z)
})
meth_fix_R <- do.call(rbind,x)
meth_fix_R2 <- as.data.frame(apply(meth_fix_R[,c(2:6)],2,as.numeric))
meth_fix_R2$gene <- unlist(meth_fix_R[,"Gene"])
head(meth_fix_R2)
dim(meth_fix_R2)


#Now import with mitch
meth_A1 <- meth_fix_A2[,c("gene","t")]
meth_A1 <- aggregate(. ~ gene, meth_A1, sum)
dim(meth_A1)

meth_R1 <- meth_fix_R2[,c("gene","t")]
meth_R1 <- aggregate(. ~ gene, meth_R1, sum)
dim(meth_R1)

l <- list("DNAm aerobic"=meth_A1, "mRNA"=transcriptome, "mRNA inactivity"=transcriptome_inactivity, "proteome"=proteome)
m <- mitch_import(x = l,geneIDcol = "gene",DEtype = "limma")

#Fetch gene sets 
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip", destfile="ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
genesets <- gmt_import("ReactomePathways.gmt")

#Calculate enrichment priority as significance
res <- mitch_calc(x=m,genesets = genesets ,priority = "significance",minsetsize = 5, resrows = 25)
head(res$enrichment_result,25)
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise")
mitch_report(res = res, outfile = "myreport_significance_Aerobic_innactivity.html")
mitch_plots(res,outfile="mycharts_significance_Aerobic_innactivity.pdf")

#save results as table
res=res$enrichment_result
res_sig=res[res$p.adjustMANOVA<0.05 ,]
write.csv(res_sig, "mitch_aerobic_inactivity_signif.csv")

#Reproduce plot using results table
res1 = mutate(res, sig=ifelse(res$p.adjustMANOVA<0.05, "FDR<0.05", "Not Sig"))

library("ggrepel")

tiff('Effect size versus significance_random_inactivity.tiff', width =8, height = 6, units = 'in', res=600)
pvsef=ggplot(res1, aes(x=`s.dist`, y=-log(`p.adjustMANOVA`)))+
  geom_point(aes(col=sig))+
  scale_color_manual(values=c("#1BD2F7", "gray"), name = " ")+
  theme_minimal()+
  labs(y="-log(p.adjustedMANOVA)", x="Effect size")
pvsef
dev.off()

#Reproduce HeatMap
d=ncol(res$input_profile)

if (d==1) {
  formatted<-t(as.data.frame(res$analysis_metrics[ c(1,2,3,4,5 ) ]))
} else if (d==2) {
  unformatted<-t(as.data.frame(res$analysis_metrics[ c(2,3,4,5,11,12 )  ]))
  formatted<-unformatted
  formatted[1:4]<-as.character(round(as.numeric(unformatted[1:4]) , digits=0))
  formatted[5:6]<-as.character(round(as.numeric(unformatted[5:6]) , digits=5))
} else if (d>2) {
  formatted<-t(as.data.frame(res$analysis_metrics[ c(2,3,4,5 )  ]))
}

colnames(formatted)="Profile metrics"

formatted %>%
  kbl(caption="Profiling data metrics") %>%
  kable_styling("hover", full_width = FALSE)

hmapx<-head(res$enrichment_result[,4:(4+d-1)],14)
rownames(hmapx)<-head(res$enrichment_result$set,14)
colnames(hmapx)<-gsub("^s.","",colnames(hmapx))
my_palette <- colorRampPalette(c("#CC99FF", "#1BD2F7","white","#1AEDD8"))(n = 10)
hmapx=as.matrix(hmapx)
colnames(hmapx)=c("DNAm Aerobic", "mRNA exercise", "mRNA innactivity","Proteome")
hmapx=hmapx[,c("DNAm Aerobic", "mRNA exercise", "Proteome","mRNA innactivity")]

tiff('Heatplot integration Aerobic and innactivity.tiff',
     width =200,
     height = 250,
     units = 'mm',
     res=600)
heatmap.2(hmapx,scale="none",Rowv = T, Colv = FALSE ,margin=c(10, 25),
          cexRow=0.8,trace="none",cexCol=1,col=my_palette,
          dendrogram="none")
dev.off()

####Enrichment with cluster profiler####
# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(pheatmap)
library(clusterProfiler) # for PEA analysis
library('org.Hs.eg.db')
library(DOSE)
library(enrichplot) # for visualisations
library(ggupset) # for visualisations

## Function: Adjacency matrix to list ####
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

# Get the genes that are present in your dataframe
genes_in_proteome <- all_results_P$gene
# Read in the .gmt file
pwl2 <- read.gmt("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/ReactomePathways.gmt") 
# Subset to the genes that are present in our dataset
pwl2_R <- pwl2[pwl2$gene %in% genes_in_proteome,] 
# Save the filtered background gene set
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise")
filename <- 'proteome_path_R.RDS'
saveRDS(pwl2_R, filename)

# Remove non-significant genes
prot_path <- all_results_P[all_results_P$FDR < 0.05, ]
# Substitute names so they are annotated nicely in the heatmap later
# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
prot_results_list <- split(prot_path, prot_path$diff)

## Run ClusterProfiler -----------------------------------------------
# Settings
name_of_comparison <- 'Proteome' # for our filename
background_genes <- 'reactome' # for our filename
bg_genes_prot <- readRDS("proteome_path_R.RDS") # read in the background genes
padj_cutoff <- 0.05 # p-adjusted threshold, used to filter out pathways
genecount_cutoff <- 5 # minimum number of genes in the pathway, used to filter out pathways

res_prot <- lapply(names(prot_results_list),
                    function(x) enricher(gene = prot_results_list[[x]]$gene,
                                         TERM2GENE = bg_genes_prot))
names(res_prot) <- names(prot_results_list)
#Convert the list of enrichResults for each sample_pattern to a dataframe with the pathways
res_df_prot <- lapply(names(res_prot), function(x) rbind(res_prot[[x]]@result))
names(res_df_prot) <- names(res_prot)
res_df_prot <- do.call(rbind, res_df_prot)
head(res_df_prot)

res_df_prot <- res_df_prot %>% mutate(minuslog10padj = -log10(p.adjust),
                                        diffexpressed = gsub('\\.REACTOME.*$', '', rownames(res_df_prot)))

# Subset to those pathways that have p adj < cutoff and gene count > cutoff (you can also do this in the enricher function)
target_pws <- unique(res_df_prot$ID[res_df_prot$p.adjust < padj_cutoff ]) # select only target pathways have p adjusted < 0.05 and at least 6 genes
res_df_prot <- res_df_prot[res_df_prot$ID %in% target_pws, ]

write.csv(res_df_prot, 'proteome_resclusterp.csv', row.names = FALSE)

res_df_down_prot <- res_df_prot %>% filter(grepl("Downexpressed", diffexpressed)) %>% 
  dplyr::select(!c('minuslog10padj', 'diffexpressed')) 
rownames(res_df_down_prot) <- res_df_down_prot$ID

# For visualisation purposes, let's shorten the pathway names
res_df_down_prot$Description <- gsub('(H|h)iv', 'HIV', 
                                    gsub('pd 1', 'PD-1',
                                         gsub('ecm', 'ECM', 
                                              gsub('(I|i)nterleukin', 'IL', 
                                                   gsub('(R|r)na', 'RNA', 
                                                        gsub('(D|d)na', 'DNA',
                                                             gsub(' i ', ' I ', 
                                                                  gsub('(A|a)tp ', 'ATP ', 
                                                                       gsub('(N|n)adh ', 'NADH ', 
                                                                            gsub('(N|n)ad ', 'NAD ',
                                                                                 gsub('t cell', 'T cell',
                                                                                      gsub('b cell', 'B cell',
                                                                                           gsub('built from .*', ' (...)',
                                                                                                gsub('mhc', 'MHC',
                                                                                                     gsub('mhc class i', 'MHC I', 
                                                                                                          gsub('mhc class ii', 'MHC II', 
                                                                                                               stringr::str_to_sentence(
                                                                                                                 gsub('_', ' ',  
                                                                                                                      gsub('REACTOME_', '', res_df_down_prot$Description)))))))))))))))))))

enrichres_down_prot <- new("enrichResult",
                          readable = FALSE,
                          result = res_df_down_prot,
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2,
                          organism = "human",
                          ontology = "UNKNOWN",
                          gene = all_results_P$gene,
                          keytype = "UNKNOWN",
                          universe = unique(bg_genes_prot$gene),
                          gene2Symbol = character(0),
                          geneSets = bg_genes_prot)
class(enrichres_down_prot)


# Dotplot
tiff('dot plot down pathways proteome.tiff',
     width =200,
     height = 300,
     units = 'mm',
     res=300)
dotplot(enrichres_down_prot, showCategory = 20) + ggtitle("Top pathways from down expressed DEPs")
dev.off()

tiff('dot plot up pathways proteome.tiff',
     width =200,
     height = 300,
     units = 'mm',
     res=300)
dotplot(enrichres_up_prot, showCategory = 20) + ggtitle("Top pathways from up expressed DEPs")
dev.off()

####intersect DMRs with sig proteins and mRNAs####
#VennDiagram of DNA meth and Proteins
library("ggvenn")
library('ggplot2')

x <- list("DMGs Aerobic" = genesZscore,
          "DMGs Resistance" = DNAm_R$Gene,
          "DEPs" = signif_meta_Proteomics$MarkerName,
          "DEGs exercise"=Significant_Extrameta$gene,
          "DEGs innactivity"=Significant_transcrip_innactivity$`Gene symbol`)
value=c(DNAm_A$Gene,DNAm_R$Gene, signif_meta_Proteomics$MarkerName, Significant_Extrameta$gene, Significant_transcrip_innactivity$`Gene symbol`)
value=unique(value)

data_venn <- data.frame(value = value,    # Create example data frame
                        DMRs_A = FALSE,
                        DMRs_R = FALSE,
                        Proteins = FALSE,
                        mRNA = FALSE,
                        mRNA_inn= FALSE)
data_venn$DMRs_R <- data_venn$value %in% x$`DMGs Resistance`
data_venn$DMRs_A <- data_venn$value %in% x$`DMGs Aerobic`
data_venn$Proteins <- data_venn$value %in% x$DEPs
data_venn$mRNA <- data_venn$value %in% x$`DEGs exercise`
data_venn$mRNA_inn <- data_venn$value %in% x$`DEGs innactivity`

ggplot(data_venn, aes(A=DMRs_A, B=Proteins, C=mRNA, D=mRNA_inn)) +
  geom_venn(fill_color = c("#CC99FF", "#1BD2F7","#FF99FF","#1AEDD8", "#FFFF99"), 
            stroke_size = 1, set_names = c("DMRs Aerobic", "Proteins", "mRNA exercise", "mRNA innactivity"), set_name_size = 5, fill_alpha = 0.3, text_size = 5, 
            show_percentage = FALSE, text_color = "#333333")+
  theme_classic()+
  theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank())

tiff('VennDiagram Aerobic.tiff', width =10, height = 8, units = 'in', res=600)
ggplot(data_venn, aes(A=DMRs_A, B=Proteins, C=mRNA, D=mRNA_inn)) +
  geom_venn(fill_color = c("#CC99FF", "#1BD2F7","#FF99FF","#1AEDD8", "#FFFF99"), 
            stroke_size = 1, set_names = c("DMRs Aerobic", "Proteins", "mRNA exercise", "mRNA innactivity"), set_name_size = 5, fill_alpha = 0.3, text_size = 5, 
            show_percentage = FALSE, text_color = "#333333")+
  theme_classic()+
  theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank())
dev.off()

ggplot(data_venn, aes(A=DMRs_R, B=Proteins, C=mRNA, D=mRNA_inn)) +
  geom_venn(fill_color = c("#CC99FF", "#1BD2F7","#FF99FF","#1AEDD8", "#FFFF99"), 
            stroke_size = 1, set_names = c("DMRs Resistance", "Proteins", "mRNA exercise", "mRNA innactivity"), set_name_size = 5, fill_alpha = 0.3, text_size = 5, 
            show_percentage = FALSE, text_color = "#333333")+
  theme_classic()+
  theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank())

tiff('VennDiagram Resistance.tiff', width =10, height = 8, units = 'in', res=600)
ggplot(data_venn, aes(A=DMRs_R, B=Proteins, C=mRNA, D=mRNA_inn)) +
  geom_venn(fill_color = c("#CC99FF", "#1BD2F7","#FF99FF","#1AEDD8", "#FFFF99"), 
            stroke_size = 1, set_names = c("DMRs Resistance", "Proteins", "mRNA exercise", "mRNA innactivity"), set_name_size = 5, fill_alpha = 0.3, text_size = 5, 
            show_percentage = FALSE, text_color = "#333333")+
  theme_classic()+
  theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank())
dev.off()
