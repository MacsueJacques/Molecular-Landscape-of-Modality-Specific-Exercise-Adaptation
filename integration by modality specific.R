setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise/Transcriptomics")
library(tidyverse)

#Load extrameta results for exercise
#Extrameta
load("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise/Transcriptomics/meta_analysis_results.RData")

#Subset the results from chronic, muscle
chronic_muscle <- all_meta_analysis_res$`longterm,muscle`
#Have a look at number of genes in muscle for chronic training
length(chronic_muscle) #14,309 genes present
#Convert to gene symbol
library(org.Hs.eg.db)
library(annotate)
names(chronic_muscle) <- getSYMBOL(names(chronic_muscle),
                                   data='org.Hs.eg')
#Remove NA gene symbols from results
chronic_muscle <- chronic_muscle[-which(is.na(names(chronic_muscle)))]

#Extract model type (i.e. simple model or moderators), effect size, standard error and p-value for each gene.
#For each gene, the first entry is the best fit
my_fn <- function(l)
{
  tib <- tibble(mod_type = sub(".*:", "", names(l)[1]),
                ES = l[[1]]$coeffs["intrcpt","beta"],
                SE = l[[1]]$coeffs["intrcpt","se"],
                p = l[[1]]$coeffs["intrcpt","pval"],
                I2 = l[[1]]$I2)
  return(tib)
}

results <- lapply(chronic_muscle,
                  my_fn)

#Unlist and bind into a single table
all_results <- results %>%
  purrr::reduce(bind_rows)

#Add gene name
all_results <- all_results %>%
  mutate(Gene = names(chronic_muscle),
         ES = log(2^ES))%>%
  dplyr::rename(`log2FC after training` = ES,
                `SE training` = SE,
                `P-value training` = p)

#Save results
write.table(all_results,
            file="Extrameta all genes results.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

#Investigate whether the age-related genes changed with training tended to be blunted in older age
library(tidyverse)
#Load RE meta-analysis results from metafor for age
age <- read_delim('C:/Users/e5103562/OneDrive - Victoria University/TWAS meta analysis of age/metafor_muscle_age_RE.txt')
age <- age %>%
  dplyr::rename(`log2FC per year of age` = ES_RE,
                SE = SE_RE,
                `P-value age` = pval_RE,
                `FDR age` = FDR)%>%
  dplyr::mutate(`log2FC per year of age` = `log2FC per year of age`/100,
                SE = SE/100)%>%
  arrange(`P-value age`)
DEGs <- subset(age,`FDR age` < 0.005)

#Load ExtraMETA results for training
extrameta <- read_tsv("Extrameta all genes results.txt")%>%
  #dplyr::filter(Gene %in% DEGs$Gene)%>%
  dplyr::rename(`log2FC after training` = `logFC after training`,
                `Moderators?` = mod_type)%>%
  mutate(`Moderators?` = str_replace(`Moderators?`,
                                     "base_model",
                                     "None")) %>%
  mutate(`Moderators?` = str_replace(`Moderators?`,
                                     "prop_males",
                                     "Sex")) %>%
  mutate(`Moderators?` = str_replace(`Moderators?`,
                                     "avg_age",
                                     "Age")) %>%
  mutate(`Moderators?` = str_replace(`Moderators?`,
                                     "training",
                                     "Exercise modality")) %>%
  mutate(`Moderators?` = str_replace(`Moderators?`,
                                     "time",
                                     "Length of training")) %>%
  mutate(`FDR training` = p.adjust(`P-value training`))%>%
  dplyr::select(Gene,
                `log2FC after training`,
                `SE training`,
                `P-value training`,
                `FDR training`,
                `Moderators?`)

extrameta_resistance <- extrameta%>%
  filter(str_detect(`Moderators?`, "Exercise modality"))

extrameta_aerobic <- extrameta%>%
  filter(!str_detect(`Moderators?`, "Exercise modality"))

extrameta <- extrameta %>%
  filter(`Moderators?`!="None")


extrameta_resistance_sig=extrameta_resistance%>%
  filter(`FDR training`<0.05)

extrameta_aerobic_sig=extrameta_aerobic%>%
  filter(`FDR training`<0.05)

# Library
library(ggplot2)
library(dplyr)
library(hrbrthemes)

#Then load meth results for downstream analyses
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
library(readxl)
transcrip_innactivity=read_excel("Differentially expressed genes with disuse from 2021.xlsx", col_names=TRUE,
                                 col_types=c("numeric","text", "numeric","numeric","numeric","numeric","numeric","numeric","numeric"))

Significant_transcrip_innactivity=transcrip_innactivity[transcrip_innactivity$`Stouffer P-value: FDR` <0.05,]
signif_trans_inact=unique(Significant_transcrip_innactivity$`Gene symbol`)

transcrip_innactivity2=transcrip_innactivity%>%
  mutate(t=`LFC: mean`/`LFC: SE`)%>%
  drop_na()

sum(Significant_transcrip_innactivity$`LFC: mean`>0)*100/nrow(Significant_transcrip_innactivity)
sum(Significant_transcrip_innactivity$`LFC: mean`<0)*100/nrow(Significant_transcrip_innactivity)

#Load methylation tables
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise/Methylation/Metal exercise")

meth_aerobic=read.delim("METAL_muscle_Aerobic.txt")
meth_resistance=read.delim("METAL_muscle_resistance.txt")

####Use Mitch to contrast results####
library("mitch")
library(tidyverse)

proteome <- Meta_Proteomics
proteome$gene <- proteome$MarkerName
head(proteome)
str(proteome)

transcriptome_A <- extrameta_aerobic
colnames(transcriptome_A)[1] <- "gene"
head(transcriptome_A)
str(transcriptome_A)

transcriptome_R <- extrameta_resistance
colnames(transcriptome_R)[1] <- "gene"
head(transcriptome_R)
str(transcriptome_R)

transcriptome_inactivity <- transcrip_innactivity2
colnames(transcriptome_inactivity)[2] <- "gene"
head(transcriptome_inactivity)
str(transcriptome_inactivity)

#methylation results
meth_A1 <- meth_aerobic%>%
  mutate(t=Effect.size/SE)%>%
  dplyr::select("CpG","Effect.size","SE","P.value","FDR", "t","Annotated.gene.s.")
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
  dplyr::select("CpG","Effect.size","SE","P.value","FDR", "t","Annotated.gene.s.")
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

proteome=proteome%>%
  mutate(t=proteome$Effect/proteome$StdErr)%>%
  dplyr::select("gene","t")

transcriptome_all=extrameta%>%
  mutate(t=extrameta$`log2FC after training`/extrameta$`SE training`)%>%
  rename(gene="Gene")%>%
  dplyr::select("gene","t")

transcriptome_A=transcriptome_A%>%
  mutate(t=transcriptome_A$`log2FC after training`/transcriptome_A$`SE training`)%>%
  dplyr::select("gene","t")

transcriptome_R=transcriptome_R%>%
  mutate(t=transcriptome_R$`log2FC after training`/transcriptome_R$`SE training`)%>%
  dplyr::select("gene","t")

#transcriptome

l_rna <- list("mRNA aerobic"=transcriptome_A, "mRNA resistance"=transcriptome_R)
m_rna <- mitch_import(x = l_rna,geneIDcol = "gene",DEtype = "limma")

res_rna<- mitch_calc(x=m_meth,genesets = genesets ,priority = "significance",minsetsize = 5, resrows = 25)
head(res_meth$enrichment_result,25)

#Reproduce HeatMap
d=ncol(res_meth$input_profile)

if (d==1) {
  formatted<-t(as.data.frame(res_meth$analysis_metrics[ c(1,2,3,4,5 ) ]))
} else if (d==2) {
  unformatted<-t(as.data.frame(res_meth$analysis_metrics[ c(2,3,4,5,11,12 )  ]))
  formatted<-unformatted
  formatted[1:4]<-as.character(round(as.numeric(unformatted[1:4]) , digits=0))
  formatted[5:6]<-as.character(round(as.numeric(unformatted[5:6]) , digits=5))
} else if (d>2) {
  formatted<-t(as.data.frame(res_meth$analysis_metrics[ c(2,3,4,5 )  ]))
}

colnames(formatted)="Profile metrics"

formatted %>%
  kbl(caption="Profiling data metrics") %>%
  kable_styling("hover", full_width = FALSE)

hmapx<-head(res_meth$enrichment_result[,4:(4+d-1)],25)
rownames(hmapx)<-head(res_meth$enrichment_result$set,25)
colnames(hmapx)<-gsub("^s.","",colnames(hmapx))
my_palette <- colorRampPalette(c("#CC99FF", "white","#1AEDD8"))(n = 10)
hmapx=as.matrix(hmapx)
colnames(hmapx)=c("DNAm Arobic","DNAm Resistance")
hmapx=hmapx[,c("DNAm Resistance", "mRNA Resistance","mRNA innactivity")]

tiff('Heatplot integration methylation.tiff',
     width =200,
     height = 250,
     units = 'mm',
     res=600)
heatmap.2(hmapx,scale="none",Rowv = T, Colv = FALSE ,margin=c(10, 25),
          cexRow=0.8,trace="none",cexCol=1,col=my_palette,
          dendrogram="none")
dev.off()

l_meth <- list("DNAm aerobic"=meth_A1, "DNAm resistance"=meth_R1)
m_meth <- mitch_import(x = l_meth,geneIDcol = "gene",DEtype = "limma")

res_meth<- mitch_calc(x=m_meth,genesets = genesets ,priority = "significance",minsetsize = 5, resrows = 25)
head(res_meth$enrichment_result,25)

#Reproduce HeatMap
d=ncol(res_meth$input_profile)

if (d==1) {
  formatted<-t(as.data.frame(res_meth$analysis_metrics[ c(1,2,3,4,5 ) ]))
} else if (d==2) {
  unformatted<-t(as.data.frame(res_meth$analysis_metrics[ c(2,3,4,5,11,12 )  ]))
  formatted<-unformatted
  formatted[1:4]<-as.character(round(as.numeric(unformatted[1:4]) , digits=0))
  formatted[5:6]<-as.character(round(as.numeric(unformatted[5:6]) , digits=5))
} else if (d>2) {
  formatted<-t(as.data.frame(res_meth$analysis_metrics[ c(2,3,4,5 )  ]))
}

colnames(formatted)="Profile metrics"

formatted %>%
  kbl(caption="Profiling data metrics") %>%
  kable_styling("hover", full_width = FALSE)

hmapx<-head(res_meth$enrichment_result[,4:(4+d-1)],25)
rownames(hmapx)<-head(res_meth$enrichment_result$set,25)
colnames(hmapx)<-gsub("^s.","",colnames(hmapx))
my_palette <- colorRampPalette(c("#CC99FF", "white","#1AEDD8"))(n = 10)
hmapx=as.matrix(hmapx)
colnames(hmapx)=c("DNAm Arobic","DNAm Resistance")
hmapx=hmapx[,c("DNAm Resistance", "mRNA Resistance","mRNA innactivity")]

tiff('Heatplot integration methylation.tiff',
     width =200,
     height = 250,
     units = 'mm',
     res=600)
heatmap.2(hmapx,scale="none",Rowv = T, Colv = FALSE ,margin=c(10, 25),
          cexRow=0.8,trace="none",cexCol=1,col=my_palette,
          dendrogram="none")
dev.off()

l_A <- list("DNAm aerobic"=meth_A1, "mRNA"=transcriptome_all, "mRNA inactivity"=transcriptome_inactivity, "proteins"=proteome)
m_A <- mitch_import(x = l_A,geneIDcol = "gene",DEtype = "limma")

#Fetch gene sets 
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip", destfile="ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
genesets <- gmt_import("ReactomePathways.gmt")

#Calculate enrichment priority as significance
res_A <- mitch_calc(x=m_A,genesets = genesets ,priority = "significance",minsetsize = 5, resrows = 25)
head(res_A$enrichment_result,25)

#Reproduce HeatMap
d=ncol(res_A$input_profile)

if (d==1) {
  formatted<-t(as.data.frame(res_A$analysis_metrics[ c(1,2,3,4,5 ) ]))
} else if (d==2) {
  unformatted<-t(as.data.frame(res_A$analysis_metrics[ c(2,3,4,5,11,12 )  ]))
  formatted<-unformatted
  formatted[1:4]<-as.character(round(as.numeric(unformatted[1:4]) , digits=0))
  formatted[5:6]<-as.character(round(as.numeric(unformatted[5:6]) , digits=5))
} else if (d>2) {
  formatted<-t(as.data.frame(res_A$analysis_metrics[ c(2,3,4,5 )  ]))
}

colnames(formatted)="Profile metrics"

formatted %>%
  kbl(caption="Profiling data metrics") %>%
  kable_styling("hover", full_width = FALSE)

hmapx<-head(res_A$enrichment_result[,4:(4+d-1)],25)
rownames(hmapx)<-head(res_A$enrichment_result$set,25)
colnames(hmapx)<-gsub("^s.","",colnames(hmapx))
my_palette <- colorRampPalette(c("#CC99FF", "white","#1AEDD8"))(n = 10)
hmapx=as.matrix(hmapx)
colnames(hmapx)=c("DNAm Arobic","mRNA", "mRNA innactivity", "Proteome")
hmapx=hmapx[,c("DNAm Arobic","mRNA","Proteome", "mRNA innactivity")]

tiff('Heatplot integration aerobic.tiff',
     width =200,
     height = 250,
     units = 'mm',
     res=600)
heatmap.2(hmapx,scale="none",Rowv = T, Colv = FALSE ,margin=c(10, 25),
          cexRow=0.8,trace="none",cexCol=1,col=my_palette,
          dendrogram="none")
dev.off()

l_R <- list("DNAm resistance"=meth_R1, "mRNA"=transcriptome_all, "mRNA inactivity"=transcriptome_inactivity)
m_R <- mitch_import(x = l_R,geneIDcol = "gene",DEtype = "limma")

res_R <- mitch_calc(x=m_R,genesets = genesets ,priority = "significance",minsetsize = 5, resrows = 25)
head(res_R$enrichment_result,25)


#Reproduce HeatMap
d=ncol(res_R$input_profile)

if (d==1) {
  formatted<-t(as.data.frame(res_R$analysis_metrics[ c(1,2,3,4,5 ) ]))
} else if (d==2) {
  unformatted<-t(as.data.frame(res_R$analysis_metrics[ c(2,3,4,5,11,12 )  ]))
  formatted<-unformatted
  formatted[1:4]<-as.character(round(as.numeric(unformatted[1:4]) , digits=0))
  formatted[5:6]<-as.character(round(as.numeric(unformatted[5:6]) , digits=5))
} else if (d>2) {
  formatted<-t(as.data.frame(res_R$analysis_metrics[ c(2,3,4,5 )  ]))
}

colnames(formatted)="Profile metrics"

formatted %>%
  kbl(caption="Profiling data metrics") %>%
  kable_styling("hover", full_width = FALSE)

hmapx<-head(res_R$enrichment_result[,4:(4+d-1)],25)
rownames(hmapx)<-head(res_R$enrichment_result$set,25)
colnames(hmapx)<-gsub("^s.","",colnames(hmapx))
my_palette <- colorRampPalette(c("#CC99FF","white","#1AEDD8"))(n = 10)
hmapx=as.matrix(hmapx)
colnames(hmapx)=c("DNAm Resistance","mRNA exercise", "mRNA innactivity")
hmapx=hmapx[,c("DNAm Resistance","mRNA exercise", "mRNA innactivity")]

tiff('Heatplot integration resistance.tiff',
     width =200,
     height = 250,
     units = 'mm',
     res=600)
heatmap.2(hmapx,scale="none",Rowv = T, Colv = FALSE ,margin=c(10, 25),
          cexRow=0.8,trace="none",cexCol=1,col=my_palette,
          dendrogram="none")
dev.off()


res<- mitch_calc()

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise")
mitch_report(res = res, outfile = "myreport_significance_Aerobic_innactivity.html")
mitch_plots(res,outfile="mycharts_significance_Aerobic_innactivity.pdf")

#save results as table
resA=res_A$enrichment_result
res_sig_A=resA[resA$p.adjustMANOVA<0.05 ,]
write.csv(res_sig_A, "mitch_aerobic_inactivity_signif_noprot.csv")

#save results as table
resR=res_R$enrichment_result
res_sig_R=resR[resR$p.adjustMANOVA<0.05 ,]
write.csv(res_sig, "mitch_resistance_inactivity_signif.csv")

#Reproduce plot using results table
res1 = mutate(res, sig=ifelse(res$p.adjustMANOVA<0.05, "FDR<0.05", "Not Sig"))

library("ggrepel")

tiff('Effect size versus significance_Aerobic_inactivity.tiff', width =8, height = 6, units = 'in', res=600)
pvsef=ggplot(res1, aes(x=`s.dist`, y=-log(`p.adjustMANOVA`)))+
  geom_point(aes(col=sig))+
  scale_color_manual(values=c("#1BD2F7", "gray"), name = " ")+
  theme_minimal()+
  labs(y="-log(p.adjustedMANOVA)", x="Effect size")
pvsef
dev.off()

#Reproduce HeatMap
d=ncol(res_R$input_profile)

if (d==1) {
  formatted<-t(as.data.frame(res_R$analysis_metrics[ c(1,2,3,4,5 ) ]))
} else if (d==2) {
  unformatted<-t(as.data.frame(res_R$analysis_metrics[ c(2,3,4,5,11,12 )  ]))
  formatted<-unformatted
  formatted[1:4]<-as.character(round(as.numeric(unformatted[1:4]) , digits=0))
  formatted[5:6]<-as.character(round(as.numeric(unformatted[5:6]) , digits=5))
} else if (d>2) {
  formatted<-t(as.data.frame(res_R$analysis_metrics[ c(2,3,4,5 )  ]))
}

colnames(formatted)="Profile metrics"

formatted %>%
  kbl(caption="Profiling data metrics") %>%
  kable_styling("hover", full_width = FALSE)

hmapx<-head(res_R$enrichment_result[,4:(4+d-1)],25)
rownames(hmapx)<-head(res_R$enrichment_result$set,25)
colnames(hmapx)<-gsub("^s.","",colnames(hmapx))
my_palette <- colorRampPalette(c("#CC99FF", "#1BD2F7","white","#1AEDD8"))(n = 5)
hmapx=as.matrix(hmapx)
colnames(hmapx)=c("DNAm Resistance", "mRNA Resistance", "mRNA innactivity")
hmapx=hmapx[,c("DNAm Resistance", "mRNA Resistance","mRNA innactivity")]

tiff('Heatplot integration Resistance and innactivity.tiff',
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

#Load annotation
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Annotation")
annotation <- read.delim("Annotation.txt")

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets Exercise/Methylation/Tables & Figures")
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

DNAm_R <- meth_resistance %>%
  filter(FDR<0.05)%>%
  dplyr::rename(Gene = `Annotated.gene.s.`)%>%
  separate_rows(Gene,
                sep=";")%>%
  dplyr::select(Gene,
                `Effect.size`)

DNAm_A <- meth_aerobic %>%
  filter(FDR<0.05)%>%
  dplyr::rename(Gene = `Annotated.gene.s.`)%>%
  separate_rows(Gene,
                sep=";")%>%
  dplyr::select(Gene,
                `Effect.size`)


x <- list("DMGs Aerobic" = genesZscore,
          "DMGs Resistance" = DNAm_R$Gene,
          "DEPs" = signif_meta_Proteomics$MarkerName,
          "DEGs Aerobic"=extrameta_aerobic_sig$Gene,
          "DEGs Resistance"=extrameta_resistance_sig$Gene,
          "DEGs innactivity"=Significant_transcrip_innactivity$`Gene symbol`)
value=c(DNAm_A$Gene,DNAm_R$Gene, signif_meta_Proteomics$MarkerName,extrameta_aerobic_sig$Gene, extrameta_resistance_sig$Gene, Significant_transcrip_innactivity$`Gene symbol`)
value=unique(value)

data_venn <- data.frame(value = value,    # Create example data frame
                        DMRs_A = FALSE,
                        DMRs_R = FALSE,
                        Proteins = FALSE,
                        mRNA_A = FALSE,
                        mRNA_R=FALSE,
                        mRNA_inn= FALSE)
data_venn$DMRs_R <- data_venn$value %in% x$`DMGs Resistance`
data_venn$DMRs_A <- data_venn$value %in% x$`DMGs Aerobic`
data_venn$Proteins <- data_venn$value %in% x$DEPs
data_venn$mRNA_A <- data_venn$value %in% x$`DEGs Aerobic`
data_venn$mRNA_R <- data_venn$value %in% x$`DEGs Resistance`
data_venn$mRNA_inn <- data_venn$value %in% x$`DEGs innactivity`

ggplot(data_venn, aes(A=DMRs_A, B=Proteins, C=mRNA_A, D=mRNA_inn)) +
  geom_venn(fill_color = c("#CC99FF", "#1BD2F7","#FF99FF","#1AEDD8", "#FFFF99"), 
            stroke_size = 1, set_names = c("DMRs Aerobic", "Proteins", "mRNA Aerobic", "mRNA innactivity"), set_name_size = 5, fill_alpha = 0.3, text_size = 5, 
            show_percentage = FALSE, text_color = "#333333")+
  theme_classic()+
  theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank())

tiff('VennDiagram Aerobic.tiff', width =10, height = 8, units = 'in', res=600)
ggplot(data_venn, aes(A=DMRs_A, B=Proteins, C=mRNA_A, D=mRNA_inn)) +
  geom_venn(fill_color = c("#CC99FF", "#1BD2F7","#FF99FF","#1AEDD8", "#FFFF99"), 
            stroke_size = 1, set_names = c("DMRs Aerobic", "Proteins", "mRNA Aerobic", "mRNA innactivity"), set_name_size = 5, fill_alpha = 0.3, text_size = 5, 
            show_percentage = FALSE, text_color = "#333333")+
  theme_classic()+
  theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank())
dev.off()

ggplot(data_venn, aes(A=DMRs_R, B=mRNA_R, D=mRNA_inn)) +
  geom_venn(fill_color = c("#CC99FF", "#1BD2F7","#1AEDD8"), 
            stroke_size = 1, set_names = c("DMRs Resistance", "mRNA Resistance", "mRNA innactivity"), set_name_size = 5, fill_alpha = 0.3, text_size = 5, 
            show_percentage = FALSE, text_color = "#333333")+
  theme_classic()+
  theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank())

tiff('VennDiagram Resistance.tiff', width =10, height = 8, units = 'in', res=600)
ggplot(data_venn, aes(A=DMRs_R, B=mRNA_R, D=mRNA_inn)) +
  geom_venn(fill_color = c("#CC99FF", "#1BD2F7","#1AEDD8"), 
            stroke_size = 1, set_names = c("DMRs Resistance", "mRNA Resistance", "mRNA innactivity"), set_name_size = 5, fill_alpha = 0.3, text_size = 5, 
            show_percentage = FALSE, text_color = "#333333")+
  theme_classic()+
  theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank())
dev.off()
