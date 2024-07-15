setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation")

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
#Load data
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/E-MTAB-11282")
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
sd(pheno$`VO2 max Pre`)
mean(pheno$`VO2 max Pre`)

#EWAS
design=model.matrix(~Age+
                      Gender*`VO2 max Pre`+
                      Timepoint+
                      `Lean or Obese`,
                    pheno)
library(limma)
corfit <- duplicateCorrelation(M,
                               design,
                               block=pheno$ID)
cor = corfit$consensus.correlation #0.150
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

coef = "`VO2 max Pre`"
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
  labs(title="Distribution of raw p-values for VO2max DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "E-MTAB-11282_VO2max",
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
tiff('E-MTAB-11282_pvalhist_DMPs.tiff',
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

coef = "GenderMale"
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
  labs(title="Distribution of raw p-values for Sex DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "E-MTAB-11282_Sex",
               directory = getwd())

coef = "GenderMale:`VO2 max Pre`"
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
  labs(title="Distribution of raw p-values for Sex|VO2max DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "E-MTAB-11282_Sex_VO2max",
               directory = getwd())

#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="E-MTAB-11282_M_res.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")


####GSE223786####
#Load data
#Get the GEO dataset
library(tidyverse)
library(GEOquery)
da=read.delim("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/GSE223786/GSE223786_GEO_submission_sparta_matrix_processed.txt")
rownames(da)=da$ID_REF
B=da[,-c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47)]

pheno=read_xlsx("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/GSE223786/Individualised data for GEO samples.xlsx")
colnames(pheno)=c("Paticipant ID","Age","Weight", "Body Fat","Height","VO2max", "VO2max_l", "Timepoint","Diet")

library(minfi)
M <- logit2(B)
library(tidyverse)
nrow(pheno)
mean(pheno$Age)
sd(pheno$Age)
range(pheno$Age)
sd(pheno$`VO2max_Pre`)
mean(pheno$`VO2max_Pre`)

#EWAS
design=model.matrix(~Age+
                      VO2max_Pre+ Diet+
                      Timepoint,
                    pheno)
library(limma)
corfit <- duplicateCorrelation(M,
                               design,
                               block=pheno$`Participant ID`)
cor = corfit$consensus.correlation #0.14218
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

coef = "VO2max_Pre"
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
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/GSE223786")
tiff('GSE223786_pvalhist_DMPs.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for VO2max DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "GSE223786_VO2max",
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
tiff('GSE223786_pvalhist_DMPs_AGE.tiff',
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
               dataset_label = "GSE223786_Age",
               directory = getwd())

#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="GSE223786_M_res.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")

####Sharples GSE TBC####
#Load data
#Get the GEO dataset
library(tidyverse)
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Sharples GSE TBC")
B=read.delim("NORMALISED_ALL PROBES_matrix_processed_for GEO.txt")
rownames(B)=B$ID_REF
B=B[,-c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41)]

pheno=read_xlsx("Sample ID for GEO.xlsx")
colnames(pheno)=c("Array Barcode no","Sample number","Participant","Sex","Condition","Age","Weight","Height","VO2max","VO2max_l")

library(minfi)
M <- logit2(B)
library(tidyverse)
nrow(pheno)
mean(pheno$Age)
sd(pheno$Age)
range(pheno$Age)
sd(pheno$`VO2max_Pre`)
mean(pheno$`VO2max_Pre`)

#EWAS
design=model.matrix(~Age+
                      VO2max_Pre+ Sex+
                      Condition,
                    pheno)
library(limma)
corfit <- duplicateCorrelation(M,
                               design,
                               block=pheno$Participant)
cor = corfit$consensus.correlation #0.0552
fit1_M <- lmFit(M,
                design,
                block=pheno$Participant,
                correlation=cor)
fit2_M <- eBayes(fit1_M)
fit1_B <- lmFit(B,
                design,
                block=pheno$Participant,
                correlation=cor)
fit2_B <- eBayes(fit1_B)

coef = "VO2max_Pre"
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
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Sharples GSE TBC")
tiff('SharplesTBC_pvalhist_DMPs.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for VO2max DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "SharplesTBC_VO2max",
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
tiff('SharplesTBC_pvalhist_DMPs_AGE.tiff',
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
               dataset_label = "SharplesTBC_Age",
               directory = getwd())

#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="SharplesTBC_M_res.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")



####GSE244359#### 
#not completed
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/GSE244359")

B <- read.table("GSE244359_series_matrix.txt", comment.char = "!")
rownames(B)=B$V1
colnames(B)=B[1,]
B=B[-1,-1]
pheno <- read_xlsx("Physiological_data_Macsue.xlsx")





####Gene SMART 1####
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Gene SMART")

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
mean(pheno$VO2max_rel_pre, na.rm=T)
sd(pheno$VO2max_rel_pre, na.rm=T)

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
                      Sex*VO2max_rel_pre+
                      Timepoint+
                      Batch+Intervention,
                    pheno)
library(limma)
corfit <- duplicateCorrelation(M,
                               design,
                               block=pheno$ID)
cor = corfit$consensus.correlation #0.231
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

coef = "VO2max_rel_pre"
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
tiff('Gene_SMART_pvalhist_DMPs_VO2max.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for VO2max DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "Gene SMART_VO2max",
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

coef = "SexM"
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
tiff('Gene_SMART_pvalhist_DMPs_Sex.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Sex DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "Gene_SMART_Sex",
               directory = getwd())

coef = "SexM:VO2max_rel"
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
tiff('Gene_SMART_pvalhist_DMPs_Sex_VO2max.tiff',
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
               dataset_label = "Gene_SMART_Sex_VO2max",
               directory = getwd())


#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="Gene_SMART_M_res.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")

####EXACT####
#Load data
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/EXACT")

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
mean(pheno$VO2max_pre)
sd(pheno$VO2max_pre)


#EWAS
design=model.matrix(~Age+
                      Timepoint+
                      VO2max_pre,
                    pheno)
library(limma)
corfit <- duplicateCorrelation(M,
                               design,
                               block=pheno$ID)
cor = corfit$consensus.correlation #0.155
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

coef = "VO2max_pre"
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
tiff('EXACT_pvalhist_DMPs.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for age DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "EXACT",
               directory = getwd())

#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="EXACT_M_res.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")


####GSE213029####
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/GSE213029")

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
mean(pheno$VO2max_pre)
sd(pheno$VO2max_pre)

#EWAS
design=model.matrix(~Age+
                      Condition+
                      Time+
                      VO2max_pre,
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

coef = "VO2max_pre"
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
tiff('GSE213029_pvalhist_DMPs.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for VO2max DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "GSE213029_VO2max",
               directory = getwd())

#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="GSE213029_M_res.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")

####GSE60655####
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/GSE60655")

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
                      Sex*VO2max_pre+
                      Timepoint,
                    pheno)
library(limma)
corfit <- duplicateCorrelation(M,
                               design,
                               block=pheno$ID)
cor = corfit$consensus.correlation #0.216
#cor = 0.216
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

coef = "VO2max_pre"
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
tiff('GSE60655_pvalhist_DMPs_VO2max.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for VO2max DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "GSE60655_VO2max",
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

coef = "SexMale"
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
tiff('GSE60655_pvalhist_DMPs_Sex.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Sex DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "GSE60655_Sex",
               directory = getwd())

coef = "SexMale:VO2max_pre"
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
tiff('GSE60655_pvalhist_DMPs_Sex_VO2max.tiff',
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
               dataset_label = "GSE60655_Sex_VO2max",
               directory = getwd())



#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="GSE60655_M_res.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")


####epiH####
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/epiH")

B <- read.table("EpiH Normalized beta values after batch correction.txt")
library(minfi)
M <- logit2(B)
library(tidyverse)
library(readxl)
pheno <- read_xlsx("Phenotype data - EPI-H samples -For Andrew.xlsx")

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
                      `Condition`,
                    pheno)
library(limma)
corfit <- duplicateCorrelation(M,
                               design,
                               block=pheno$ID_rep)
cor = corfit$consensus.correlation #0.1
#cor = 0.1
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

coef = "`Fitness level (ml/min/kg)`"

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
tiff('epiH_pvalhist_DMPs_VO2max.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for VO2max DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "epiH_VO2max",
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

coef = "SexM"
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
tiff('EpiH_pvalhist_DMPs_Sex.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Sex DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "EpiH_Sex",
               directory = getwd())

coef = "SexM:`Fitness level (ml/min/kg)`"
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
tiff('EpiH_pvalhist_DMPs_Sex_VO2max.tiff',
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
               dataset_label = "EpiH_Sex_VO2max",
               directory = getwd())


#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="epiH_M_res.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")

####trainSMART####
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/trainSMART")

B <- read.table("beta trainsmart.txt")
library(minfi)
M <- logit2(B)
library(tidyverse)
library(readxl)
pheno <- read_xlsx("phenotype-train-smart.xlsx")%>%
  rename(Timepoint="Pre=T0_Post=T2")

#reorder beta and M table to match pheno ID
new_order=c(pheno$Sample_ID)
B=B[,new_order]
M=M[,new_order]

nrow(pheno)
mean(pheno$Age)
sd(pheno$Age)
range(pheno$Age)
sum(pheno$Sex=="M")/nrow(pheno)
mean(pheno$VO2_pre)
sd(pheno$VO2_pre)

#EWAS
design=model.matrix(~Age+
                      Sex*VO2_pre+
                      Timepoint,
                    pheno)
library(limma)
corfit <- duplicateCorrelation(M,
                               design,
                               block=pheno$CODE)
cor = corfit$consensus.correlation #0.1294
#cor = 0.1294
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

coef = "VO2_pre"

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
tiff('trainSMART_pvalhist_DMPs_VO2max.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for VO2max DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "trainSMART_VO2max",
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

coef = "SexM"
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
tiff('trainSMART_pvalhist_DMPs_Sex.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Sex DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "trainSMART_Sex",
               directory = getwd())

coef = "SexM:VO2_pre"
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
tiff('trainSMART_pvalhist_DMPs_Sex_VO2max.tiff',
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
               dataset_label = "trainSMART_Sex_VO2max",
               directory = getwd())



#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="trainSMART_M_res.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")

####Wellderly####
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Wellderly")

B <- read.table("beta wellderly.txt")
library(minfi)
M <- logit2(B)
library(tidyverse)
library(readxl)
pheno <- read_xlsx("phenotype-wellderly.xlsx")

#reorder beta and M table to match pheno ID
new_order=c(pheno$`Sample ID`)
B=B[,new_order]
M=M[,new_order]

nrow(pheno)
mean(pheno$Age)
sd(pheno$Age)
range(pheno$Age)
sum(pheno$Sex=="M")/nrow(pheno)
mean(pheno$V02)
sd(pheno$V02)

#EWAS
design=model.matrix(~Age+
                      Sex*V02,
                    pheno)
library(limma)
fit1_M <- lmFit(M,
                design)
fit2_M <- eBayes(fit1_M)
fit1_B <- lmFit(B,
                design)
fit2_B <- eBayes(fit1_B)

coef = "V02"

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
tiff('Wellderly_pvalhist_DMPs_VO2max.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for VO2max DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "Wellderly_VO2max",
               directory = getwd())


#Save tables for ageing
coef = "Age"

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
tiff('Wellderly_pvalhist_DMPs_ageing.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Ageing DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "Wellderly_Age",
               directory = getwd())


coef = "SexM"
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
tiff('Wellderly_pvalhist_DMPs_Sex.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Sex DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "Wellderly_Sex",
               directory = getwd())

coef = "SexM:V02"
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
tiff('Wellderly_pvalhist_DMPs_Sex_VO2max.tiff',
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
               dataset_label = "Wellderly_Sex_VO2max",
               directory = getwd())


#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="Wellderly_Age_M_res.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")

####Gene SMART EPICV2####
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Gene SMART V2")

B <- read.table("beta genesmart_950K.txt")
library(minfi)
M <- logit2(B)
library(tidyverse)
library(readxl)
pheno <- read_xlsx("phenotype-genesmart-2020.xlsx")

nrow(pheno)
mean(pheno$Age)
sd(pheno$Age)
range(pheno$Age)
sum(pheno$Sex=="M")/nrow(pheno)
pheno$V02=as.numeric(pheno$V02)
mean(pheno$V02,na.rm=TRUE)
sd(pheno$V02,na.rm=TRUE)

pheno=pheno[1:10,]
B=B[,1:10]
M=M[,1:10]

#EWAS
design=model.matrix(~Age+
                      Sex*V02,
                    pheno)
library(limma)
fit1_M <- lmFit(M,
                design)
fit2_M <- eBayes(fit1_M)
fit1_B <- lmFit(B,
                design)
fit2_B <- eBayes(fit1_B)

coef = "V02"

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
tiff('GeneSMARTV2_pvalhist_DMPs.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for VO2max DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "GeneSMART_EPICV2_VO2max",
               directory = getwd())

#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="GeneSMART_M_res.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")

#Save tables for ageing
coef = "Age"

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
tiff('GeneSMARTV2_pvalhist_DMPs_ageing.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Ageing DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "GeneSMART_EPICV2_Age",
               directory = getwd())

coef = "SexM"
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
tiff('Gene_SMARTV2_pvalhist_DMPs_Sex.tiff',
     width =5,
     height = 3,
     units = 'in',
     res = 200)
ggplot(results,
       mapping = aes(x = `P.Value`))+
  labs(title="Distribution of raw p-values for Sex DMPs",
       x="p-value")+
  geom_histogram(color="darkblue",
                 fill="lightblue")
dev.off()

create_summary(toptable = results,
               dataset_label = "Gene_SMARTV2_Sex",
               directory = getwd())

coef = "SexM:V02"
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
tiff('Gene_SMARTV2_pvalhist_DMPs_Sex_VO2max.tiff',
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
               dataset_label = "Gene_SMARTV2_Sex_VO2max",
               directory = getwd())


#Save residuals
ResidualsMatrix <- residuals(fit2_M, M)
write.table(signif(ResidualsMatrix,digits = 4),
            file="GeneSMART_Age_M_res.txt",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep="\t")


####Bacon - inflation correction####
#List all files ending in "tbl" in the folder
directory=c("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)")
setwd(directory)
library(bacon)
library(tidyverse)
files2 <- c("GSE223786_VO2max.tbl","SharplesTBC_VO2max.tbl")
#List tbl files
files <- list.files()[grep(".tbl",list.files())]

for (f in files2)
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

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)")
E_MTAB <- read_tsv("E-MTAB-11282_VO2max_corrected.TBL")
E_MTAB2 <- left_join(E_MTAB, annotation, join_by(CPG==probeID))
E_MTAB2=E_MTAB2[!duplicated(E_MTAB2$Unique_ID),]
E_MTAB2=E_MTAB2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(E_MTAB2, #save the file
            file="E-MTAB-11282_VO2max_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

GSE223786=read_tsv("GSE223786_VO2max_corrected.tbl")			
GSE223786_2 <- left_join(GSE223786, annotation, join_by(CPG==probeID))
GSE223786_2=GSE223786_2[!duplicated(GSE223786_2$Unique_ID),]
GSE223786_2=GSE223786_2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(GSE223786_2, #save the file
            file="GSE223786_VO2max_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

SharplesTBC=read_tsv("SharplesTBC_VO2max_corrected.tbl")			
SharplesTBC_2 <- left_join(SharplesTBC, annotation, join_by(CPG==probeID))
SharplesTBC_2=SharplesTBC_2[!duplicated(SharplesTBC_2$Unique_ID),]
SharplesTBC_2=SharplesTBC_2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(SharplesTBC_2, #save the file
            file="SharplesTBC_VO2max_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")			

EXACT <- read_tsv("EXACT_corrected.TBL")
EXACT2 <- left_join(EXACT, annotation, join_by(CPG==probeID))
EXACT2=EXACT2[!duplicated(EXACT2$Unique_ID),]
EXACT2=EXACT2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(EXACT2, #save the file
            file="EXACT_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

FTC <- read_tsv("FTC_vo2max_ml_kg_min_190224_2_corrected.TBL")
FTC2 <- left_join(FTC, annotation, join_by(CPG==probeID))
FTC2=FTC2[!duplicated(FTC2$Unique_ID),]
FTC2=FTC2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(FTC2, #save the file
            file="FTC_vo2max_ml_kg_min_190224_2_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

GeneSMART <- read_tsv("Gene SMART_VO2max_corrected.TBL")
GeneSMART2 <- left_join(GeneSMART, annotation, join_by(CPG==probeID))
GeneSMART2=GeneSMART2[!duplicated(GeneSMART2$Unique_ID),]
GeneSMART2=GeneSMART2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(GeneSMART2, #save the file
            file="Gene SMART_VO2max_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

GSE60655 <- read_tsv("GSE60655_VO2max_corrected.TBL")
GSE60655_2 <- left_join(GSE60655, annotation, join_by(CPG==probeID))
GSE60655_2=GSE60655_2[!duplicated(GSE60655_2$Unique_ID),]
GSE60655_2=GSE60655_2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(GSE60655_2, #save the file
            file="GSE60655_VO2max_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

GSE213029 <- read_tsv("GSE213029_VO2max_corrected.TBL")
GSE213029_2 <- left_join(GSE213029, annotation, join_by(CPG==probeID))
GSE213029_2=GSE213029_2[!duplicated(GSE213029_2$Unique_ID),]
GSE213029_2=GSE213029_2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(GSE213029_2, #save the file
            file="GSE213029_VO2max_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

annotation2=unite(annotation, "probeID", "probeID":"EPICv2")

GeneSMARTv2 <- read_tsv("GeneSMART_EPICV2_VO2max_corrected.TBL")
GeneSMART2v2_2 <- left_join(GeneSMARTv2, annotation2, join_by(CPG==probeID))
GeneSMART2v2_2=GeneSMART2v2_2[!duplicated(GeneSMART2v2_2$Unique_ID),]
GeneSMART2v2_2=GeneSMART2v2_2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(GeneSMART2v2_2, #save the file
            file="GeneSMART_EPICV2_VO2max_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

trainSMART <- read_tsv("trainSMART_VO2max_corrected.TBL")
trainSMART2 <- left_join(trainSMART, annotation2, join_by(CPG==probeID))
trainSMART2=trainSMART2[!duplicated(trainSMART2$Unique_ID),]
trainSMART2=trainSMART2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(trainSMART2, #save the file
            file="trainSMART_VO2max_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

Wellderly <- read_tsv("Wellderly_VO2max_corrected.TBL")
Wellderly2 <- left_join(Wellderly, annotation2, join_by(CPG==probeID))
Wellderly2=Wellderly2[!duplicated(Wellderly2$Unique_ID),]
Wellderly2=Wellderly2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(Wellderly2, #save the file
            file="Wellderly_VO2max_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

EpiH <- read_tsv("epiH_VO2max_corrected.TBL")
EpiH2 <- left_join(EpiH, annotation2, join_by(CPG==probeID))
EpiH2=EpiH2[!duplicated(EpiH2$Unique_ID),]
EpiH2=EpiH2[,c("Unique_ID","ALLELE1","ALLELE2","TESTSTAT","PVALUE","EFFECTSIZE", "SE","EFFECTSIZE_CORR","SE_CORR","PVALUE_CORR","EFFECTSIZE_CORR100","SE_CORR100")]
write.table(EpiH2, #save the file
            file="epiH_VO2max_corrected.TBL",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")


####Load METAL results and annotate####
library(tidyverse)
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max data")
meta_res <- read_tsv("METAANALYSIS_VO2.TBL")
#Change name to probeID
meta_res <- dplyr::rename(meta_res,
                          Unique_ID = MarkerName)

#Add CHR and position to meta analysis results
meta_res_tot <- left_join(x = meta_res,
                          y = annotation,
                          by = "Unique_ID")
meta_res_tot<-meta_res_tot[!duplicated(meta_res_tot$Unique_ID), ]

max(meta_res_tot$HetDf)

#Calculate nb of CpGs included in analysis for each nb of studies
nb_CpGs <- c()
for (i in 1:max(meta_res_tot$HetDf))
{
  nb_CpGs <- c(nb_CpGs,nrow(meta_res_tot %>% filter(HetDf>=(i-1))))
}

cumsum <- tibble(`Number of included studies` = 1:max(meta_res_tot$HetDf),
                 `Number of CpGs present in all studies` = nb_CpGs)

ggplot(cumsum,
       mapping = aes(x = `Number of included studies`,
                     y = `Number of CpGs present in all studies`))+
  geom_point()+
  geom_line()+
  theme_classic()

#View results including in at least 5 studies
meta_res_robust <- meta_res_tot %>% filter(HetDf>=4)
padj <- p.adjust(meta_res_robust$`P-value`,method="BH")
meta_res_robust <- meta_res_robust %>%
  mutate(`Adjusted p-value` = padj)
meta_res_robust <- meta_res_robust %>%
  mutate(sig = ifelse(`Adjusted p-value`<0.05,
                      "Sig",
                      "Not Sig"))
nrow(meta_res_robust) #740727
nDMPs <- nrow(meta_res_robust %>% filter(sig=="Sig")) #DMPs = 13195


#####Arrange
sigdir <- meta_res_robust$sig
sigdir[meta_res_robust$sig=="Sig"&
         meta_res_robust$Effect<0]="Hypo"
sigdir[meta_res_robust$sig=="Sig"&
         meta_res_robust$Effect>0]="Hyper"
meta_res_robust <- meta_res_robust %>%
  mutate(sigdir = sigdir)
meta_res_robust <- meta_res_robust %>%
  arrange(desc(sigdir))

nrow(meta_res_robust %>% filter(sigdir=="Hypo"))/nDMPs #% hypo DMPs = 54% hypo

#Select relevant columns only and rename them
meta_res_robust <- meta_res_robust %>%
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
meta_res_robust <- meta_res_robust %>%
  mutate(genesUniq = replace_na(genesUniq,""),
         GeneHancer_interaction = replace_na(GeneHancer_interaction,""))

#Replace CpG island names with better words
meta_res_robust$CGIposition[grep("Shore",meta_res_robust$CGIposition)]="Shore"
meta_res_robust$CGIposition[grep("Shelf",meta_res_robust$CGIposition)]="Shelf"
meta_res_robust$CGIposition[is.na(meta_res_robust$CGIposition)]="Open sea"

#Change HetDf by number of included studies
meta_res_robust$HetDf <- meta_res_robust$HetDf +1
#Rename columns
colnames(meta_res_robust) = c("CpG",
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

nDMPs <- nrow(meta_res_robust %>% filter(Significance=="Sig")) #3100

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max data")
write.table(meta_res_robust,
            file="METAL_muscle_VO2max.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

#### Tables & Figures ####
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Tables & Figures")
library(tidyverse)
#Load METAL results that have been formatted to annotate each CpG
to_write <- meta_res_robust
to_order <- tibble(CpG = to_write$CpG)

to_write <- to_write %>%
  dplyr::rename(`% DNAm change per VO2max` = `Effect size`,
                `P-value VO2max` = `P-value`,
                `FDR VO2max` = FDR)%>%
  arrange(`P-value VO2max`)

DMPs <- subset(to_write,
               `FDR VO2max`<0.05)

#Volcano plot for EWAS RE meta-analysis of VO2max

tiff('Volcano plot_EWAS RE meta-analysis of VO2max.tiff',
     width =100,
     height = 80,
     units = 'mm',
     res=300)
ggplot() +
  geom_point(data = to_write,
             aes(`% DNAm change per VO2max`,
                 -log10(`P-value VO2max`),
                 col=Direction),
             size=0.5)+
  scale_color_manual(values=c("#1BD2F7", "#1AEDD8","grey"))+
  labs(y="-log10(p-value)")+
  theme_classic()+
  theme(legend.position = "none")
dev.off()


#Show histogram of effect sizes for DMPs
tiff('Distribution of effect size in DMPs_EWAS RE meta-analysis of VO2max.tiff',
     width =100,
     height = 65,
     units = 'mm',
     res=300)
ggplot(data = DMPs,
       aes(x=`% DNAm change per VO2max`,
           fill=Direction)) +
  geom_histogram(colour = "black")+
  scale_fill_manual(values=c("red","blue"))+
  lims(x=c(min(to_write$`% DNAm change per VO2max`),
           max(to_write$`% DNAm change per VO2max`)))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

#Forest plot for top hypo and hyper DMP
library(readxl)
dataset_summary <- read_excel('//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Tables & Figures/DNAm datasets.xlsx',
                              na = c("","NA"))%>%
  dplyr::rename(prop_males = `% Male`)
dataset_summary=dataset_summary[1:10,]

#Create empty list
L.names <- meta_res_robust$`Unique ID`
L <- vector("list", length(L.names))
names(L) <- L.names

#Initiate list
tib <- tibble(Unique_ID = meta_res_robust$`Unique ID` ,
              CpG= meta_res_robust$CpG,
              Dataset = rep("Meta-analysis",nrow(meta_res_robust)),
              ES = meta_res_robust$`Effect size`,
              SE = meta_res_robust$`SE`,
              PVAL = meta_res_robust$`P-value`,
              FDR = meta_res_robust$FDR,
              Type = rep("Meta-analysis",nrow(meta_res_robust)))
tib <- tib %>%
  mutate_if(is.numeric,
            signif,
            digits = 2)

#Initiate the list from meta-analysis
L <- split(tib, seq(nrow(tib))) #splits each row of the tibble into a list component
names(L) <- meta_res_robust$`Unique ID`

#List tbl files
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max data")
files <- list.files(pattern="_corrected.tbl",
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
  summary <- tibble(Unique_ID = file$Unique_ID,
                    Dataset = rep(s,nrow(file)),
                    ES = file$EFFECTSIZE_CORR100,
                    SE = file$SE_CORR100,
                    PVAL = file$PVALUE_CORR,
                    FDR = file$FDR,
                    Type = rep("Individual study",nrow(file)))
  summary <- summary %>%
    mutate_if(is.numeric,
              signif,
              digits = 2)
  
  summary <- summary %>%
    filter(Unique_ID %in% names(L))
  subL <- split(summary, seq(nrow(summary)))
  names(subL) = summary$Unique_ID
  
  #Merge the pieces of the two lists that are in common
  L2 <- L[names(subL)]
  L3 <- L[setdiff(names(L),names(subL))]
  listinter <- Map(bind_rows,
                   L2,
                   subL)
  L <- c(listinter,
         L3)
}

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Tables & Figures")

saveRDS(L,
        file = "METAL_muscle_VO2max.rds")


L <- read_rds("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Tables & Figures/METAL_muscle_VO2max.rds")

top_hypo <- DMPs %>%
  filter(`% DNAm change per VO2max` < 0) %>%
  dplyr::slice(1)%>%
  pull(`Unique ID`) #CpG "cg17999652"
tib <- left_join(L[[top_hypo]],
                 dataset_summary)%>%
  arrange(-`n`)


top_hyper <- DMPs %>%
  filter(`% DNAm change per VO2max` > 0) %>%
  dplyr::slice(1)%>%
  pull(`Unique ID`) #CpG "17174613"


cg <- L[["chr18_79871624_79871626"]]
cg$CpG

cg$Dataset=c("Meta-analysis","E-MTAB-11282","EpiH","EXACT", "Finnish Twin Cohort (FTC)","Gene SMART","Gene SMART2","CAUSE","GSE60655","train SMART","Wellderly")
Datasets$Dataset=c("Gene SMART","Finnish Twin Cohort (FTC)","EpiH",  "train SMART", "EXACT","GSE60655","CAUSE","E-MTAB-11282","Wellderly","Gene SMART2", "Meta-analysis" )

cg <- left_join(cg, Datasets, by = "Dataset")


cg$N <- as.numeric(cg$n)
N_meta <- sum(cg$N, na.rm = TRUE)
cg$N[is.na(cg$N)]<-sum(cg$N, na.rm = T)

#cg$Dataset=c("Meta-analysis", "Gene SMART","Finnish Twin Cohort (FTC)", "EXACT", "CAUSE", "E-MTAB-11282","GSE60655","EpiH","train SMART", "Gene SMART2","Welderly")

cg_Meta <- filter(cg, Dataset == "Meta-analysis")
cg <- cg %>% filter(Dataset != "Meta-analysis")
?tiff
tiff("cg17174613_forest.tiff",
     height = 11,
     width = 7,
     unit = 'in',
     res = 200)
forest(cg$ES,
       sei = cg$SE,
       slab = cg$Dataset,
       xlab = "% DNAm change per year",
       psize = 1,
       header = "Dataset",
       ilab = cg$N,
       ilab.xpos = -1.6,
       order = -(cg$N),
       cex = 0.6)
addpoly(cg_Meta$ES,
        sei = cg_Meta$SE,
        cex = 0.6,
        mlab = "Meta-analysis",
        col = "Purple",
        border = "Purple")
text(-1.6, 53,
     "N",
     cex = 0.8)
dev.off()



cg <- L[["chr15_51366359_51366361"]]
cg$CpG

Datasets <- tib %>% select(Dataset,n)

cg <- left_join(cg, Datasets, by = "Dataset")


cg$N <- as.numeric(cg$n)
N_meta <- sum(cg$N, na.rm = TRUE)
cg$N[is.na(cg$N)]<-sum(cg$N, na.rm = T)

cg$Dataset=c("Meta-analysis", "Gene SMART","Finnish Twin Cohort (FTC)", "EXACT", "CAUSE", "E-MTAB-11282","GSE60655","EpiH",
             "train SMART", "Gene SMART2","Welderly")

cg_Meta <- filter(cg, Dataset == "Meta-analysis")
cg <- cg %>% filter(Dataset != "Meta-analysis")
?tiff
tiff("cg17999652_forest.tiff",
     height = 11,
     width = 7,
     unit = 'in',
     res = 200)
forest(cg$ES,
       sei = cg$SE,
       slab = cg$Dataset,
       xlab = "% DNAm change per year",
       psize = 1,
       header = "Dataset",
       ilab = cg$N,
       ilab.xpos = -1.6,
       order = -(cg$N),
       cex = 0.6)
addpoly(cg_Meta$ES,
        sei = cg_Meta$SE,
        cex = 0.6,
        mlab = "Meta-analysis",
        col = "Purple",
        border = "Purple")
text(-1.6, 53,
     "N",
     cex = 0.8)
dev.off()






#CpG <- "cg22833204"

tib <- left_join(L[[top_hyper]],
                 dataset_summary)%>%
  arrange(-`n`)

tiff('C:/Users/e5103562/OneDrive - Victoria University/Forest plot_sex moderator STAT1.tiff',
     width =150,
     height = 180,
     units = 'mm',
     res=300)
forest(meta_RE(tib),
       header="Dataset ID",
       xlab = "% DNAm change per VO2max",
       digits=c(2L,4L),
       slab = tib %>%
         pull(Dataset))
)
dev.off()

#Plot ES as a function of % males
library(ggrepel)
tiff('C:/Users/e5103562/OneDrive - Victoria University/ES vs prop males_sex moderator STAT1.tiff',
     width =150,
     height = 100,
     units = 'mm',
     res=300)
ggplot(tib,
       mapping = aes(x=prop_males,
                     y=ES,
                     col = prop_males))+
  geom_point(size=3)+
  geom_smooth(method="lm")+
  labs(x = "% of males in dataset",
       y = "% DNAm change per year of age")+
  geom_text_repel(aes(label = Dataset),
                  col = "black",
                  #size = 2,
                  #box.padding = unit(0.45, "lines"),
                  #point.padding = unit(0.45, "lines"),
                  max.overlaps = 30)+
  theme_bw()

dev.off()

#Only keep age-related DMPs
to_write <- DMPs




