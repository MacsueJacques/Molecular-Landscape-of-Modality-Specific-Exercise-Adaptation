####Load METAL results and annotate####
library(tidyverse)
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)")
meta_res <- read_tsv("METAANALYSIS1.TBL")
#Change name to probeID
meta_res <- dplyr::rename(meta_res,
                          Unique_ID = MarkerName)

#Load annotation table
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Annotation")
annotation <- read.delim("Annotation.txt")


setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)")

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
nrow(meta_res_robust) #720489
nDMPs <- nrow(meta_res_robust %>% filter(sig=="Sig")) #DMPs = 9822


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

nrow(meta_res_robust %>% filter(sigdir=="Hypo"))/nDMPs #% hypo DMPs = 58% hypo

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
                Unique_ID,
                E107,E108)
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
                              "Unique ID",
                              "Chromatin state in male skeletal muscle", "Chromatin state in female skeletal muscle")

#remove the letters in front of Chrom States
meta_res_robust$`Chromatin state in male skeletal muscle` =ifelse(grepl("1_TssA", meta_res_robust$`Chromatin state in male skeletal muscle`), "Active TSS", meta_res_robust$`Chromatin state in male skeletal muscle`)
meta_res_robust$`Chromatin state in male skeletal muscle`=ifelse(grepl("2_TssAFlnk", meta_res_robust$`Chromatin state in male skeletal muscle`), "Flanking Active TSS", meta_res_robust$`Chromatin state in male skeletal muscle`)
meta_res_robust$`Chromatin state in male skeletal muscle`=ifelse(grepl("3_TxFlnk", meta_res_robust$`Chromatin state in male skeletal muscle`), "Trans at gene 5 and 3 primer", meta_res_robust$`Chromatin state in male skeletal muscle`)
meta_res_robust$`Chromatin state in male skeletal muscle`=ifelse(grepl("4_Tx", meta_res_robust$`Chromatin state in male skeletal muscle`), "Strong transcription", meta_res_robust$`Chromatin state in male skeletal muscle`)
meta_res_robust$`Chromatin state in male skeletal muscle`=ifelse(grepl("5_TxWk", meta_res_robust$`Chromatin state in male skeletal muscle`), "Weak transcription", meta_res_robust$`Chromatin state in male skeletal muscle`)
meta_res_robust$`Chromatin state in male skeletal muscle`=ifelse(grepl("6_EnhG", meta_res_robust$`Chromatin state in male skeletal muscle`), "Genic enhancers", meta_res_robust$`Chromatin state in male skeletal muscle`)
meta_res_robust$`Chromatin state in male skeletal muscle`=ifelse(grepl("7_Enh", meta_res_robust$`Chromatin state in male skeletal muscle`), "Enhancers", meta_res_robust$`Chromatin state in male skeletal muscle`)
meta_res_robust$`Chromatin state in male skeletal muscle`=ifelse(grepl("8_ZNF/Rpts", meta_res_robust$`Chromatin state in male skeletal muscle`), "ZNF Genes & Repeats", meta_res_robust$`Chromatin state in male skeletal muscle`)
meta_res_robust$`Chromatin state in male skeletal muscle`=ifelse(grepl("9_Het", meta_res_robust$`Chromatin state in male skeletal muscle`), "Heterochromatin", meta_res_robust$`Chromatin state in male skeletal muscle`)
meta_res_robust$`Chromatin state in male skeletal muscle`=ifelse(grepl("10_TssBiv", meta_res_robust$`Chromatin state in male skeletal muscle`), "Bivalent/Poised TSS", meta_res_robust$`Chromatin state in male skeletal muscle`)
meta_res_robust$`Chromatin state in male skeletal muscle`=ifelse(grepl("11_BivFlnk", meta_res_robust$`Chromatin state in male skeletal muscle`), "Flanking bivalent TSS/Enh", meta_res_robust$`Chromatin state in male skeletal muscle`)
meta_res_robust$`Chromatin state in male skeletal muscle`=ifelse(grepl("12_EnhBiv", meta_res_robust$`Chromatin state in male skeletal muscle`), "Bivalent enhancer", meta_res_robust$`Chromatin state in male skeletal muscle`)
meta_res_robust$`Chromatin state in male skeletal muscle`=ifelse(grepl("13_ReprPC", meta_res_robust$`Chromatin state in male skeletal muscle`), "Repressed PolyComb", meta_res_robust$`Chromatin state in male skeletal muscle`)
meta_res_robust$`Chromatin state in male skeletal muscle`=ifelse(grepl("14_ReprPCWk", meta_res_robust$`Chromatin state in male skeletal muscle`), "Weak Repressed PolyComb", meta_res_robust$`Chromatin state in male skeletal muscle`)
meta_res_robust$`Chromatin state in male skeletal muscle`=ifelse(grepl("15_Quies", meta_res_robust$`Chromatin state in male skeletal muscle`), "Quiescent/Low", meta_res_robust$`Chromatin state in male skeletal muscle`)

#remove the letters in front of Chrom States
meta_res_robust$`Chromatin state in female skeletal muscle`=ifelse(grepl("1_TssA", meta_res_robust$`Chromatin state in female skeletal muscle`), "Active TSS", meta_res_robust$`Chromatin state in female skeletal muscle`)
meta_res_robust$`Chromatin state in female skeletal muscle`=ifelse(grepl("2_TssAFlnk", meta_res_robust$`Chromatin state in female skeletal muscle`), "Flanking Active TSS", meta_res_robust$`Chromatin state in female skeletal muscle`)
meta_res_robust$`Chromatin state in female skeletal muscle`=ifelse(grepl("3_TxFlnk", meta_res_robust$`Chromatin state in female skeletal muscle`), "Trans at gene 5 and 3 primer", meta_res_robust$`Chromatin state in female skeletal muscle`)
meta_res_robust$`Chromatin state in female skeletal muscle`=ifelse(grepl("4_Tx", meta_res_robust$`Chromatin state in female skeletal muscle`), "Strong transcription", meta_res_robust$`Chromatin state in female skeletal muscle`)
meta_res_robust$`Chromatin state in female skeletal muscle`=ifelse(grepl("5_TxWk", meta_res_robust$`Chromatin state in female skeletal muscle`), "Weak transcription", meta_res_robust$`Chromatin state in female skeletal muscle`)
meta_res_robust$`Chromatin state in female skeletal muscle`=ifelse(grepl("6_EnhG", meta_res_robust$`Chromatin state in female skeletal muscle`), "Genic enhancers", meta_res_robust$`Chromatin state in female skeletal muscle`)
meta_res_robust$`Chromatin state in female skeletal muscle`=ifelse(grepl("7_Enh", meta_res_robust$`Chromatin state in female skeletal muscle`), "Enhancers", meta_res_robust$`Chromatin state in female skeletal muscle`)
meta_res_robust$`Chromatin state in female skeletal muscle`=ifelse(grepl("8_ZNF/Rpts", meta_res_robust$`Chromatin state in female skeletal muscle`), "ZNF Genes & Repeats", meta_res_robust$`Chromatin state in female skeletal muscle`)
meta_res_robust$`Chromatin state in female skeletal muscle`=ifelse(grepl("9_Het", meta_res_robust$`Chromatin state in female skeletal muscle`), "Heterochromatin", meta_res_robust$`Chromatin state in female skeletal muscle`)
meta_res_robust$`Chromatin state in female skeletal muscle`=ifelse(grepl("10_TssBiv", meta_res_robust$`Chromatin state in female skeletal muscle`), "Bivalent/Poised TSS", meta_res_robust$`Chromatin state in female skeletal muscle`)
meta_res_robust$`Chromatin state in female skeletal muscle`=ifelse(grepl("11_BivFlnk", meta_res_robust$`Chromatin state in female skeletal muscle`), "Flanking bivalent TSS/Enh", meta_res_robust$`Chromatin state in female skeletal muscle`)
meta_res_robust$`Chromatin state in female skeletal muscle`=ifelse(grepl("12_EnhBiv", meta_res_robust$`Chromatin state in female skeletal muscle`), "Bivalent enhancer", meta_res_robust$`Chromatin state in female skeletal muscle`)
meta_res_robust$`Chromatin state in female skeletal muscle`=ifelse(grepl("13_ReprPC", meta_res_robust$`Chromatin state in female skeletal muscle`), "Repressed PolyComb", meta_res_robust$`Chromatin state in female skeletal muscle`)
meta_res_robust$`Chromatin state in female skeletal muscle`=ifelse(grepl("14_ReprPCWk", meta_res_robust$`Chromatin state in female skeletal muscle`), "Weak Repressed PolyComb", meta_res_robust$`Chromatin state in female skeletal muscle`)
meta_res_robust$`Chromatin state in female skeletal muscle`=ifelse(grepl("15_Quies", meta_res_robust$`Chromatin state in female skeletal muscle`), "Quiescent/Low", meta_res_robust$`Chromatin state in female skeletal muscle`)

nDMPs <- nrow(meta_res_robust %>% filter(Significance=="Sig")) #9822
write.table(meta_res_robust,
            file="METAL_muscle_VO2max_withoutGSv2.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

meta_res_robust=read.delim("METAL_muscle_VO2max_withoutGSv2.txt")

#How many DMPs in Active TSS regions
DMPs_ActiveTSS=meta_res_robust%>%filter(Significance=="Sig")%>%filter(`Chromatin state in male skeletal muscle`=="Active TSS")
nrow(DMPs_ActiveTSS %>% filter(Direction=="Hypo"))/440 #0.4636364
DMPs_Flanking_ActiveTSS=meta_res_robust%>%filter(Significance=="Sig")%>%filter(`Chromatin state in male skeletal muscle`=="Flanking Active TSS")
nrow(DMPs_Flanking_ActiveTSS %>% filter(Direction=="Hypo"))/1496 # 0.2760695
DMPs_Enhancers=meta_res_robust%>%filter(Significance=="Sig")%>%filter(`Chromatin state in male skeletal muscle`=="Enhancers")
nrow(DMPs_Enhancers %>% filter(Direction=="Hypo"))/3082 #0.3848151

#### Tables & Figures ####
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures")
library(tidyverse)
#Load METAL results that have been formatted to annotate each CpG

to_write <- meta_res_robust
to_order <- tibble(CpG = to_write$CpG)

to_write <- to_write %>%
  dplyr::rename(`% DNAm change per VO2max` = `Effect.size`,
                `P-value VO2max` = `P.value`,
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

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures")

#Show histogram of effect sizes for DMPs
tiff('Distribution of effect size in DMPs.tiff',
     width =100,
     height = 65,
     units = 'mm',
     res=300)
ggplot(data = DMPs,
       aes(x=`% DNAm change per VO2max`,
           fill=Direction)) +
  geom_histogram(colour = "grey")+
  scale_fill_manual(values=c("#1BD2F7","#1AEDD8"))+
  lims(x=c(min(to_write$`% DNAm change per VO2max`),
           max(to_write$`% DNAm change per VO2max`)))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

####Forest plot for top hypo and hyper DMP####
####Prep####
library(readxl)
dataset_summary <- read_excel('//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures/DNAm datasets.xlsx',
                              na = c("","NA"))%>%
  dplyr::rename(prop_males = `% Male`)
dataset_summary=dataset_summary[1:9,]

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

dataset_summary$`Dataset ID`=c("Gene SMART","Finnish Twin Cohort", "EXACT", "CAUSE", "E-MTAB-11282","GSE60655","EpiH","train SMART","Wellderly", "GSE223786","Sharples")
colnames(dataset_summary)=c("Database","Dataset", "Muscle type", "n","Age (mean ± SD)", "Age range (min - max)",
                            "VO2max (mean ± SD)","Training Type","Array","prop_males","Phenotype / Disease status","Variables included in the linear model"  )

#List tbl files
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures/Tables Forest Plot")
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

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures/Tables Forest Plot")

saveRDS(L,
        file = "METAL_muscle_VO2max_10studies.rds")

####Forestplot - start here####
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures/Tables Forest Plot")

L <- read_rds("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures/Tables Forest Plot/METAL_muscle_VO2max.rds")

dataset_summary2=dataset_summary[,c(2,4,10)]
colnames(dataset_summary2)=c("Dataset", "n","prop_males")

#1
top_hypo1 <- DMPs %>%
  filter(`% DNAm change per VO2max` < 0) %>%
  dplyr::slice(1)%>%
  pull(`Unique.ID`) #CpG "cg21130958"
top_hypo1=L[[top_hypo1]]
top_hypo1$Dataset=sub("\\..*","",top_hypo1$Dataset)

tib_hypo1 <- left_join(top_hypo1,
                       dataset_summary2)%>%
  arrange(-`n`)
tib_hypo1 [10,11]=702
tib_hypo1 [10,17]=0.41

tib_hypo1=tib_hypo1 %>%
  dplyr::select("Unique_ID","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")
#2
top_hypo2 <- DMPs %>%
  filter(`% DNAm change per VO2max` < 0) %>%
  dplyr::slice(2)%>%
  pull(`Unique.ID`) #CpG "cg04356588"
top_hypo2=L[[top_hypo2]]
top_hypo2$Dataset=sub("\\..*","",top_hypo2$Dataset)

tib_hypo2 <- left_join(top_hypo2,
                       dataset_summary2)%>%
  arrange(-`n`)
tib_hypo2 [9,11]=668
tib_hypo2 [9,17]=0.45

tib_hypo2=tib_hypo2 %>%
  dplyr::select("Unique_ID","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")


#3
top_hypo3 <- DMPs %>%
  filter(`% DNAm change per VO2max` < 0) %>%
  dplyr::slice(3)%>%
  pull(`Unique.ID`) #CpG "cg12411159"
top_hypo3=L[[top_hypo3]]
top_hypo3$Dataset=sub("\\..*","",top_hypo3$Dataset)

tib_hypo3 <- left_join(top_hypo3,
                       dataset_summary2)%>%
  arrange(-`n`)
tib_hypo3 [6,11]=496
tib_hypo3 [6,17]=0.55

tib_hypo3=tib_hypo3 %>%
  dplyr::select("Unique_ID","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")

#4
top_hypo4 <- DMPs %>%
  filter(`% DNAm change per VO2max` < 0) %>%
  dplyr::slice(4)%>%
  pull(`Unique.ID`) #CpG "cg08440987"
top_hypo4=L[[top_hypo4]]
top_hypo4$Dataset=sub("\\..*","",top_hypo4$Dataset)

tib_hypo4 <- left_join(top_hypo4,
                       dataset_summary2)%>%
  arrange(-`n`)
tib_hypo4 [9,11]=668
tib_hypo4 [9,17]=0.45

tib_hypo4=tib_hypo4 %>%
  dplyr::select("Unique_ID","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")

#5
top_hypo5 <- DMPs %>%
  filter(`% DNAm change per VO2max` < 0) %>%
  dplyr::slice(5)%>%
  pull(`Unique.ID`) #CpG "cg12916671"
top_hypo5=L[[top_hypo5]]
top_hypo5$Dataset=sub("\\..*","",top_hypo5$Dataset)

tib_hypo5 <- left_join(top_hypo5,
                       dataset_summary2)%>%
  arrange(-`n`)
tib_hypo5 [10,11]=702
tib_hypo5 [10,17]=0.41

tib_hypo5=tib_hypo5 %>%
  dplyr::select("Unique_ID","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")


top5_hypo=rbind(tib_hypo1,tib_hypo2,tib_hypo3,tib_hypo4,tib_hypo5)

top5_hypo$Unique_ID[which(top5_hypo$Unique_ID=="chr8_74986002_74986004")]<- "CRISPLD1"
top5_hypo$Unique_ID[which(top5_hypo$Unique_ID=="chr5_68243748_68243750")]<- "PIK3R1"
top5_hypo$Unique_ID[which(top5_hypo$Unique_ID=="chr5_181039678_181039680")]<- "BTNL9"
top5_hypo$Unique_ID[which(top5_hypo$Unique_ID=="chr11_128609888_128609890")]<- "ETS1"

top5_hypo<-top5_hypo%>%
  mutate(Dataset=factor(Dataset, levels= c("Meta-analysis", "Gene SMART", "Finnish Twin Cohort", "EpiH", "train SMART",
                                           "EXACT", "GSE60655", "CAUSE", "E-MTAB-11282", "Wellderly")))

tiff("Top5 Hypo CpGs.tiff",
     height = 11,
     width = 10,
     unit = 'in',
     res = 200)
ggforestplot::forestplot(
  df = top5_hypo,
  name = Unique_ID,
  estimate = ES,
  se = SE,
  pvalue = FDR,
  psignif = 0.05,
  colour = Dataset,
  shape = Dataset,
  xlab = "Odds ratio for increased VO2max (95% CI) per 1−SD decrease in CpG methylation",
  title = "Top 5 Hypo-methylated CpGs associated with VO2max",
  logodds = F,
  size=3
)+
  ggplot2::scale_shape_manual(
    values = c(23L, 21L, 21L, 21L, 21L,21L,21L,21L,21L,21L),
    labels = c("Meta-analysis", "Gene SMART", "Finnish Twin Cohort", "EpiH", "train SMART",
               "EXACT", "GSE60655", "CAUSE", "E-MTAB-11282", "Wellderly")
  )+
  ggplot2::scale_colour_manual(values = c("#00112B","#002255","#003380","#0044AA","#0055D4","#0066ff","#3771c8", "#5599FF","#80B3FF", 
                                          "#AACCFF"))
dev.off()


#1
top_hyper1 <- DMPs %>%
  filter(`% DNAm change per VO2max` >0) %>%
  dplyr::slice(1)%>%
  pull(`Unique.ID`) #CpG "cg02787506"
top_hyper1=L[[top_hyper1]]
top_hyper1$Dataset=sub("\\..*","",top_hyper1$Dataset)

tib_hyper1 <- left_join(top_hyper1,
                        dataset_summary2)%>%
  arrange(-`n`)
#tib_hyper1 [10,11]=702
#tib_hyper1 [10,17]=0.41

tib_hyper1=tib_hyper1 %>%
  select("Unique_ID","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")
#2
top_hyper2 <- DMPs %>%
  filter(`% DNAm change per VO2max` >0) %>%
  dplyr::slice(2)%>%
  pull(`Unique.ID`) #CpG "cg00546403"
top_hyper2=L[[top_hyper2]]
top_hyper2$Dataset=sub("\\..*","",top_hyper2$Dataset)

tib_hyper2 <- left_join(top_hyper2,
                        dataset_summary2)%>%
  arrange(-`n`)
#tib_hyper2 [9,11]=668
#tib_hyper2 [9,17]=0.45

tib_hyper2=tib_hyper2 %>%
  select("Unique_ID","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")


#3
top_hyper3 <- DMPs %>%
  filter(`% DNAm change per VO2max` >0) %>%
  dplyr::slice(3)%>%
  pull(`Unique.ID`) #CpG "cg06380443"
top_hyper3=L[[top_hyper3]]
top_hyper3$Dataset=sub("\\..*","",top_hyper3$Dataset)

tib_hyper3 <- left_join(top_hyper3,
                        dataset_summary2)%>%
  arrange(-`n`)
#tib_hyper3 [6,11]=496
#tib_hyper3 [6,17]=0.55

tib_hyper3=tib_hyper3 %>%
  select("Unique_ID","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")

#4
top_hyper4 <- DMPs %>%
  filter(`% DNAm change per VO2max` >0) %>%
  dplyr::slice(4)%>%
  pull(`Unique.ID`) #CpG "cg06207917"
top_hyper4=L[[top_hyper4]]
top_hyper4$Dataset=sub("\\..*","",top_hyper4$Dataset)

tib_hyper4 <- left_join(top_hyper4,
                        dataset_summary2)%>%
  arrange(-`n`)
#tib_hyper4 [9,11]=668
#tib_hyper4 [9,17]=0.45

tib_hyper4=tib_hyper4 %>%
  select("Unique_ID","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")

#5
top_hyper5 <- DMPs %>%
  filter(`% DNAm change per VO2max` >0) %>%
  dplyr::slice(5)%>%
  pull(`Unique.ID`) #CpG "cg16124410"
top_hyper5=L[[top_hyper5]]
top_hyper5$Dataset=sub("\\..*","",top_hyper5$Dataset)

tib_hyper5 <- left_join(top_hyper5,
                        dataset_summary2)%>%
  arrange(-`n`)
#tib_hyper5 [10,11]=702
#tib_hyper5 [10,17]=0.41

tib_hyper5=tib_hyper5 %>%
  select("Unique_ID","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")


top5_hyper=rbind(tib_hyper1,tib_hyper2,tib_hyper3,tib_hyper4,tib_hyper5)

top5_hyper$Unique_ID[which(top5_hyper$Unique_ID=="chr10_97540766_97540768")]<- "UBTD1"
top5_hyper$Unique_ID[which(top5_hyper$Unique_ID=="chr7_131010972_131010974")]<- "LINC-PINT"
top5_hyper$Unique_ID[which(top5_hyper$Unique_ID=="chr5_142935391_142935393")]<- "ARHGAP26"
top5_hyper$Unique_ID[which(top5_hyper$Unique_ID=="chr1_205255830_205255832")]<- "TMCC2"
top5_hyper$Unique_ID[which(top5_hyper$Unique_ID=="chr16_47396362_47396364")]<- "ITFG1"

top5_hyper<-top5_hyper%>%
  mutate(Dataset=factor(Dataset, levels= c("Meta-analysis", "Gene SMART", "Finnish Twin Cohort", "EpiH", "train SMART",
                                           "EXACT", "GSE60655", "CAUSE", "E-MTAB-11282", "Wellderly")))
fullTopCPGs=rbind(top5_hyper, top5_hypo)
fullTopCPGs<-fullTopCPGs%>%
  mutate(Dataset=factor(Dataset, levels= c("Meta-analysis", "Gene SMART", "Finnish Twin Cohort", "EpiH", "train SMART",
                                           "EXACT", "GSE60655", "CAUSE", "E-MTAB-11282", "Wellderly")))

tiff("Top5 CpGs.tiff",
     height = 11,
     width = 10,
     unit = 'in',
     res = 200)
ggforestplot::forestplot(
  df = fullTopCPGs,
  name = Unique_ID,
  estimate = ES,
  se = SE,
  pvalue = FDR,
  psignif = 0.05,
  colour = Dataset,
  shape = Dataset,
  #xlab = "Odds ratio pre change in VO2max (95% CI) per 1−SD decrease in CpG methylation",
  title = "Top 5 hypo- and hyper- methylated CpGs associated with VO2max",
  logodds = F
)+
  ggplot2::scale_shape_manual(
    values = c(23L, 21L, 21L, 21L, 21L,21L,21L,21L,21L,21L),
    labels = c("Meta-analysis", "Gene SMART", "Finnish Twin Cohort", "EpiH", "train SMART",
               "EXACT", "GSE60655", "CAUSE", "E-MTAB-11282", "Wellderly")
  )+
  ggplot2::scale_colour_manual(values = c("#00112B","#002255","#003380","#0044AA","#0055D4","#0066ff","#3771c8", "#5599FF","#80B3FF", 
                                          "#AACCFF"))
dev.off()

####Downstream analysis####
####Identify differentially methylathed regions####
library(DMRcate)
library(IRanges)

CpGs <- DMPs$CpG # the ones with a p-value below 0.05
Unique_IDs <- DMPs$`Unique.ID`


library(tidyverse)
annotation_overlap_only=filter(annotation, Unique_ID %in% Unique_IDs)
annotation_overlap_only=annotation_overlap_only[!duplicated(annotation_overlap_only$Unique_ID),]


annotation_overlap_only=left_join(annotation_overlap_only, DMPs, by=join_by(Unique_ID==`Unique.ID`))
annotation_overlap_only=mutate(annotation_overlap_only, t=annotation_overlap_only$`% DNAm change per VO2max`/annotation_overlap_only$SE)

rownames(annotation_overlap_only)=annotation_overlap_only$Unique_ID

annotated <- GenomicRanges::GRanges(annotation_overlap_only[Unique_IDs,"CpG_chrm" ], #chromosome
                                    IRanges(annotation_overlap_only[Unique_IDs,"CpG_beg"],annotation_overlap_only[Unique_IDs,"CpG_end"]), #add location #put the same location twice so it recognises as a region
                                    stat = annotation_overlap_only[Unique_IDs,"t"], #t-statistic
                                    diff = annotation_overlap_only[Unique_IDs,"% DNAm change per VO2max"], #effect size from beta 
                                    ind.fdr = annotation_overlap_only[Unique_IDs,"FDR VO2max"], #adjusted p-value
                                    is.sig = annotation_overlap_only[Unique_IDs,"FDR VO2max"] < 0.05) #p-value threshold #logical operator TRUE of FALSE so it recognises later on

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
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures")
write.csv(DMRs,
          file="DMRs with VO2max PRE.csv")
DMRs=read.csv("DMRs with VO2max PRE.csv")

library(tidyverse)

top_DMRs=DMRs%>%
  dplyr::filter(Stouffer.score<0.005&
                  Harmonic.mean.of.the.individual.component.FDRs <0.005&
                  Fisher.multiple.comparison.statistic<0.005)

#Find how many unique genes in DMRs
genes=top_DMRs$Annotated.gene.s.
genes=c(str_split(genes, ";"))
genes=unlist(genes)
genes=c(str_split(genes, ","))
genes2=unlist(genes)
genesZscore=unique(genes2)
view(genesZscore) #983 unique genes in the DMRs

#Calculate % of hypo and hyper DMRs
sum(DMRs$Mean.effect.size.in.DMR >0)*100/nrow(DMRs) #37.1
sum(DMRs$Mean.effect.size.in.DMR <0)*100/nrow(DMRs) #62.9

sum(DMPs$`% DNAm change per VO2max`>0)*100/nrow(DMPs)
sum(DMPs$`% DNAm change per VO2max`<0)*100/nrow(DMPs)


####DMRs ideogram
library(biomaRt)
library(regioneR)
library(AnnotationHub)
ahub <- AnnotationHub()
ahub["AH5086"] #human DNA CpGs islands
cpgs <- ahub[["AH5086"]]
cpgs

epiFiles <- query(ahub, c("EpigenomeRoadMap","E107", "segmentations"))
epiFiles
metadata <- epiFiles[["AH46960"]]

metadata_promoters=metadata
metadata_enhancers=metadata

View(as.data.frame(metadata_promoters@elementMetadata@listData))

metadata_promoters@elementMetadata@listData$color_code =ifelse(grepl("#FFFFFF", metadata_promoters@elementMetadata@listData$color_code), "#FFFFFF00", metadata_promoters@elementMetadata@listData$color_code)
metadata_promoters@elementMetadata@listData$color_code =ifelse(grepl("#C0C0C0", metadata_promoters@elementMetadata@listData$color_code), "#C0C0C000", metadata_promoters@elementMetadata@listData$color_code)
metadata_promoters@elementMetadata@listData$color_code =ifelse(grepl("#808080", metadata_promoters@elementMetadata@listData$color_code), "#80808000", metadata_promoters@elementMetadata@listData$color_code)
metadata_promoters@elementMetadata@listData$color_code =ifelse(grepl("#CD5C5C", metadata_promoters@elementMetadata@listData$color_code), "#CD5C5C00", metadata_promoters@elementMetadata@listData$color_code)
metadata_promoters@elementMetadata@listData$color_code =ifelse(grepl("#006400", metadata_promoters@elementMetadata@listData$color_code), "#00640000", metadata_promoters@elementMetadata@listData$color_code)
metadata_promoters@elementMetadata@listData$color_code =ifelse(grepl("#BDB76B", metadata_promoters@elementMetadata@listData$color_code), "#BDB76B00", metadata_promoters@elementMetadata@listData$color_code)
metadata_promoters@elementMetadata@listData$color_code =ifelse(grepl("#FF4500", metadata_promoters@elementMetadata@listData$color_code), "#FF450000", metadata_promoters@elementMetadata@listData$color_code)
metadata_promoters@elementMetadata@listData$color_code =ifelse(grepl("#FFFF00", metadata_promoters@elementMetadata@listData$color_code), "#F6BE00", metadata_promoters@elementMetadata@listData$color_code)
metadata_promoters@elementMetadata@listData$color_code =ifelse(grepl("#008000", metadata_promoters@elementMetadata@listData$color_code), "#00800000", metadata_promoters@elementMetadata@listData$color_code)
metadata_promoters@elementMetadata@listData$color_code =ifelse(grepl("#C2E105", metadata_promoters@elementMetadata@listData$color_code), "#C2E10500", metadata_promoters@elementMetadata@listData$color_code)
metadata_promoters@elementMetadata@listData$color_code =ifelse(grepl("#66CDAA", metadata_promoters@elementMetadata@listData$color_code), "#BDB76B00", metadata_promoters@elementMetadata@listData$color_code)
metadata_promoters@elementMetadata@listData$color_code =ifelse(grepl("#8A91D0", metadata_promoters@elementMetadata@listData$color_code), "#BDB76B00", metadata_promoters@elementMetadata@listData$color_code)
metadata_promoters@elementMetadata@listData$color_code =ifelse(grepl("#E9967A", metadata_promoters@elementMetadata@listData$color_code), "#BDB76B00", metadata_promoters@elementMetadata@listData$color_code)
metadata_promoters@elementMetadata@listData$color_code =ifelse(grepl("#32CD32", metadata_promoters@elementMetadata@listData$color_code), "#BDB76B00", metadata_promoters@elementMetadata@listData$color_code)



View(as.data.frame(metadata_promoters@elementMetadata@listData))


gene.symbols <- genesZscore
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                         filters = 'hgnc_symbol', values =gene.symbols, mart = ensembl))
seqlevelsStyle(genes) <- "UCSC"


head(genes)

library(karyoploteR) 

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures/Tables Forest Plot")

tiff('New_DMRs 1 to 10.tiff',
     width =18,
     height = 25,
     units = 'in',
     res=600)
kp <- plotKaryotype(chromosomes=c(c("chr10","chr9","chr8","chr7","chr6","chr5","chr4","chr3","chr2","chr1")), plot.type = 2)
kpPlotRegions(kp, data=metadata_promoters, col=metadata_promoters$color_code, lwd=1)
kpPlotDensity(kp, data=metadata, data.panel=1, col="#9999FFCC",window.size = 1000000, r0=0, r1=0.6 )
kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol, data.panel=2, r1=0.3, cex=1)

dev.off()

tiff('New_DMRs 11 to 22.tiff',
     width =18,
     height = 26,
     units = 'in',
     res=600)
kp <- plotKaryotype(chromosomes=c(c("chr22","chr21","chr20","chr19","chr18","chr17","chr16","chr15","chr14","chr13","chr12","chr11")), plot.type = 2)
kpPlotRegions(kp, data=metadata_promoters, col=metadata_promoters$color_code, lwd=1)
kpPlotDensity(kp, data=metadata, data.panel=1, col="#9999FFCC",window.size = 1000000, r0=0, r1=0.6 )
kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol, data.panel=2, r1=0.3, cex=1)
dev.off()



#####Chrom states plot####
#First extract the regions hypo and hyper methylated
regions_hypo <- makeGRangesFromDataFrame(DMRs %>% filter(`Maximum effect size in DMR`<0),
                                         seqnames.field="Chromosome",
                                         start.field="Position start (hg38)",
                                         end.field="Position end (hg38)")
regions_hyper <- makeGRangesFromDataFrame(DMRs %>% filter(`Maximum effect size in DMR`>0),
                                          seqnames.field="Chromosome",
                                          start.field="Position start (hg38)",
                                          end.field="Position end (hg38)")
cpgs <- GenomicRanges::GRanges(seqnames = annotation_overlap_only$CpG_chrm ,
                               ranges = IRanges::IRanges(start = annotation_overlap_only$CpG_beg,
                                                         end = annotation_overlap_only$CpG_end),
                               name = annotation_overlap_only$probeID)
overlaps_hypo <- GenomicRanges::findOverlaps(cpgs,
                                             regions_hypo)
overlaps_hyper <- GenomicRanges::findOverlaps(cpgs,
                                              regions_hyper)
DMR.cpg.hypo <- cpgs$name[from(overlaps_hypo)]
DMR.cpg.hyper <- cpgs$name[from(overlaps_hyper)]

#Load background of CpGs
#meta_res_robust <- read_tsv("Meta-analysis results with robust CpGs new.txt")
#chromHMM <- read.delim("Chromatin states annotation.txt")
#rownames(chromHMM) <- paste(chromHMM$`STATE.NO.`,chromHMM$MNEMONIC,sep="_")
chromHMM <- chromHMM %>%
  mutate(color = rgb(RED, GREEN, BLUE, max = 255))
#Add a column to meta_res_robust corresponding to Direction for DMR
Direction.DMR <- rep("Not Sig" ,nrow(annotation_overlap_only))
names(Direction.DMR) <- annotation_overlap_only$probeID
Direction.DMR[DMR.cpg.hyper] <- "Hyper"
Direction.DMR[DMR.cpg.hypo] <- "Hypo"
res_robust <- bind_cols(annotation_overlap_only,
                        Direction.DMR = Direction.DMR)
#Remove CpGs with ambiguous classification
summary_fisher <- res_robust %>%
  filter(!is.na(`E107`)) %>%
  filter(!grepl(",",`E107`)) %>%
  group_by(Direction.DMR,`E107`) %>%
  summarise(n = n()) %>%
  pivot_wider(id_cols = `E107`,
              values_from = n,
              names_from = Direction.DMR) %>%
  as.data.frame()


rownames(summary_fisher) <- c("Bivalent/Poised TSS","Flanking Bivalent TSS/Enh","Bivalent Enhancer","Repressed PolyComb",
                              "Weak Repressed PolyComb","Quiescent/Low","Active TSS","Flanking Active TSS","Transcr. at gene 5' and 3'",
                              "Strong transcription","Weak transcription","Genic enhancers","Enhancers","Heterochromatin", "ZNF genes & repeats")

summary_fisher <- summary_fisher[,-1]
summary_fisher[is.na(summary_fisher)] = 0

chisqtest <- chisq.test(summary_fisher) #p-value < 2.2e-16
chisqtest
chisqtest$stdres
library(chisq.posthoc.test)
chisq.posthoc.test(summary_fisher,
                   method = "BH")

#Graph residuals
library(corrplot)
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
res <- chisqtest$residuals
res <- res[chromHMM$DESCRIPTION,c("Hypo","Hyper","Not Sig")]
colnames(res) <- c("hypo","hyper","non-DMR")

tiff('Chromatin states enrichment male residuals.tiff',
     width =7,
     height = 5,
     units = 'in',
     res=600)
c=corrplot(t(res),
           col = rev(col2(200)),
           is.cor = FALSE,
           cl.pos = "n",
           tl.col = "black",
           tl.pos = "lt")
c
dev.off()
summary_fisher_percent <- t(t(summary_fisher) / rowSums(t(summary_fisher)) * 100)
summary_fisher_percent <- bind_cols(as_tibble(summary_fisher_percent),
                                    `Chromatin state in male skeletal muscle` = rownames(summary_fisher_percent)) %>%
  pivot_longer(cols = c("Hypo","Not Sig","Hyper"),
               names_to = "Direction.DMR",
               values_to = "perc") %>%
  mutate(`Chromatin state in male skeletal muscle` = fct_relevel(`Chromatin state in male skeletal muscle`,chromHMM$DESCRIPTION),
         Direction.DMR = fct_recode(Direction.DMR,
                                    "hypo"="Hypo",
                                    "hyper"="Hyper",
                                    "non-DMR"="Not Sig"))
summary_fisher_percent$Direction.DMR <- fct_relevel(summary_fisher_percent$Direction.DMR,
                                                    "hypo",
                                                    "hyper",
                                                    "non-DMR")

tiff('Chromatin states enrichment male.tiff',
     width =10.5,
     height = 4,
     units = 'in',
     res=600)
b=ggplot(summary_fisher_percent, aes(x=Direction.DMR,
                                     y=perc,
                                     fill=`Chromatin state in male skeletal muscle`))+
  geom_bar(stat="identity",
           color="black",
           aes(alpha=as.numeric(as.factor(Direction.DMR))))+
  facet_grid(.~`Chromatin state in male skeletal muscle`,
             labeller = labeller(`Chromatin state in male skeletal muscle` = label_wrap_gen(10)))+
  scale_fill_manual(values=chromHMM$color)+
  scale_alpha(range = c(0.4, 1))+
  labs(x="",y="% of CpGs")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust=1),
        strip.text.x = element_text(angle = 90,
                                    size = 12),
        axis.title.x=element_blank(),
        legend.position="none",
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
b
dev.off()

#Enrichment of DMPs (or DMRs) in CpG islands
res_robust$CGIposition[grep("Shore",res_robust$CGIposition)]="Shore"
res_robust$CGIposition[grep("Shelf",res_robust$CGIposition)]="Shelf"
res_robust$CGIposition[is.na(res_robust$CGIposition)]="Open sea"

#Remove CpGs with ambiguous classification
summary_fisher_I <- res_robust %>%
  filter(!is.na(`CGIposition`)) %>%
  filter(!grepl(",",`CGIposition`)) %>%
  group_by(Direction.DMR,`CGIposition`) %>%
  summarise(n = n()) %>%
  pivot_wider(id_cols = `CGIposition`,
              values_from = n,
              names_from = Direction.DMR) %>%
  as.data.frame()
rownames(summary_fisher_I) <- summary_fisher_I$CGIposition
summary_fisher_I <- summary_fisher_I[,-1]
summary_fisher_I=summary_fisher_I[c("Island", "Shore", "Shelf", "Open sea"),]
colnames(summary_fisher_I)=c("hyper", "hypo", "non-DMR")
summary_fisher_I=summary_fisher_I[,c("hypo", "hyper", "non-DMR")]

chisqtest <- chisq.test(summary_fisher_I) #p-value < 2.2e-16
chisqtest
chisqtest$stdres
library(chisq.posthoc.test)
chisq.posthoc.test(summary_fisher_I,
                   method = "BH")

#Graph residuals
library(corrplot)
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
res <- chisqtest$residuals


tiff('CpG islands enrichment male residuals.tiff',
     width =3,
     height = 3,
     units = 'in',
     res=600)
c=corrplot(t(res),
           col = rev(col2(200)),
           is.cor = FALSE,
           cl.pos = "n",
           tl.col = "black",
           tl.pos = "lt")
c
dev.off()
summary_fisher_percent2 <- t(t(summary_fisher_I) / rowSums(t(summary_fisher_I)) * 100)

summary_fisher_percent2 <- bind_cols(as_tibble(summary_fisher_percent2),
                                     `CpG Island in male skeletal muscle` = rownames(summary_fisher_percent2)) %>%
  pivot_longer(cols = c("hypo",
                        "hyper",
                        "non-DMR"),
               names_to = "Direction.DMR",
               values_to = "perc")

summary_fisher_percent2$Direction.DMR <- fct_relevel(summary_fisher_percent2$Direction.DMR,
                                                     "hypo", "hyper",
                                                     "non-DMR")
summary_fisher_percent2$`CpG Island in male skeletal muscle` <- fct_relevel(summary_fisher_percent2$`CpG Island in male skeletal muscle`,
                                                                            "Island",
                                                                            "Shore",
                                                                            "Shelf",
                                                                            "Open sea")
tiff('CpG Islands enrichment male.tiff',
     width =3,
     height = 3,
     units = 'in',
     res=600)
b=ggplot(summary_fisher_percent2, aes(x=Direction.DMR,
                                      y=perc,
                                      fill=`CpG Island in male skeletal muscle`))+
  geom_bar(stat="identity",
           color="black",
           aes(alpha=as.numeric(as.factor(Direction.DMR))))+
  facet_grid(.~`CpG Island in male skeletal muscle`,
             labeller = labeller(`CpG Island in male skeletal muscle` = label_wrap_gen(10)))+
  scale_fill_manual(values=c("#663399", "#3399FF", "#FFCCFF","#FF6666"))+
  scale_alpha(range = c(0.4, 1))+
  labs(x="",y="% of CpGs")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust=1),
        strip.text.x = element_text(angle = 90,
                                    size = 12),
        axis.title.x=element_blank(),
        legend.position="none",
        panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
b
dev.off()


####Perform Gene Set Enrichment Analysis to identify any enriched pathways####
library(tidyverse)
#Load METAL results that have been formatted to annotate each CpG
to_write <- meta_res_robust
to_order <- tibble(CpG = to_write$CpG)

#Load gene sets
library(mitch)
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/MSigDB files/MSigDB/msigdb_v2022.1.Hs_GMTs")
genesets <- list.files(pattern = ".gmt",
                       recursive = TRUE)
cgp <- genesets[which(str_detect(genesets,
                                 fixed("cgp"))==TRUE)]
cp <- genesets[which(str_detect(genesets,
                                fixed("cp"))==TRUE)]
go <- genesets[which(str_detect(genesets,
                                fixed(".go."))==TRUE)]
hpo <- genesets[which(str_detect(genesets,
                                 fixed("hpo"))==TRUE)]
c8 <- genesets[which(str_detect(genesets,
                                fixed("c8"))==TRUE)]

cgp <- gmt_import(cgp[which(str_detect(cgp,
                                       "entrez")==TRUE)])
cp <- gmt_import(cp[which(str_detect(cp,
                                     "biocarta|kegg|pid|reactome|wiki|symbols",
                                     negate = TRUE)==TRUE)])
go <- gmt_import(go[which(str_detect(go,
                                     ".bp|.cc|.mf|symbols",
                                     negate = TRUE)==TRUE)])
hpo <- gmt_import("c5.hpo.v2022.1.Hs.entrez.gmt")
hpo_symbols<-gmt_import("c5.hpo.v2022.1.Hs.symbols.gmt")
c8 <- gmt_import(c8[which(str_detect(c8,
                                     "entrez")==TRUE)])
c8 <- c8[which(str_detect(names(c8),
                          "SKELETAL_MUSCLE")==TRUE)]
c8 <- c8[which(str_detect(names(c8),
                          "FETAL",
                          negate=TRUE)==TRUE)]



####Enrichment for each gene set####
library(missMethyl)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)


go_enrichment <- gometh(sig.cpg = to_write %>% filter(FDR<0.05) %>%pull(CpG),
                        all.cpg = to_write$CpG,
                        collection = "GO",
                        array.type = "EPIC",
                        plot.bias = FALSE,
                        prior.prob = TRUE,
                        equiv.cpg = TRUE,
                        fract.counts = TRUE,
                        sig.genes = TRUE) 
go_enrichment$ID <- rownames(go_enrichment)
go_enrichment <- go_enrichment %>%
  arrange(P.DE)
View(go_enrichment)
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures")
write.table(go_enrichment,
            file="GSEA_meta-analysis RE_GO.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

write_csv(go_enrichment,
          file="GSEA_meta-analysis RE_GO.csv")

#Graph without info about genes inside circles
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures")
library(tidyverse)
go_enrichment <- read_csv("GSEA_meta-analysis RE_GO.csv")%>%
  mutate(Description = str_extract(ID,
                                   "(?<=_)(.+)"))%>%
  mutate(Description = str_replace_all(Description,
                                       "_",
                                       " "))%>%
  mutate(Description = str_to_sentence(Description))

library(ggrepel)
library(pals)
p <- ggplot(data = go_enrichment,
            mapping = aes(x = 100*DE/N,
                          y = -log10(P.DE),
                          label = Description,
                          col = FDR<0.05
            ))+
  ylim(0,12)+
  geom_point()+
  scale_color_manual(values = c("Black","Red"))+
  geom_label_repel(data = go_enrichment %>% filter(FDR<0.05),
                   size = 2,
                   box.padding = unit(0.45, "lines"),
                   point.padding = unit(0.45, "lines"),
                   max.overlaps = 20)+
  labs(x = "% of genes in GO term that are\ndifferentially methylated",
       y = "Significance (-log10(pvalue))")+
  theme_bw()+
  theme(legend.position = "none")
p

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures")
tiff('GSEA_meta-analysis RE_GO.tiff',
     width =100,
     height = 100,
     units = 'mm',
     res=600)
p
dev.off()

#enrichment for hpo
hpo_enrichment <- gsameth(#sig.cpg = to_write %>% filter(`FDR age`<0.005) %>%pull(CpG),
  sig.cpg = DMPs$CpG,           
  all.cpg = to_write$CpG,
  collection = hpo,
  array.type = "EPIC",
  plot.bias = FALSE,
  prior.prob = TRUE,
  #anno = RSanno,
  equiv.cpg = TRUE,
  fract.counts = TRUE,
  sig.genes = TRUE) 
go_enrichment$ID <- rownames(go_enrichment)
go_enrichment <- go_enrichment %>%
  arrange(P.DE)
View(go_enrichment)
setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
write.table(c8_enrichment,
            file="GSEA_meta-analysis RE_C8.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")


####Enrichment analysis with mitch DNA meth only####
library("mitch")
library(tidyverse)

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures")
meth=read.delim("METAL_muscle_VO2max_withoutGSv2.txt")

meth1 <- meth%>%
  mutate(t=Effect.size/SE)%>%
  dplyr::select("CpG","Effect.size","SE","P.value","FDR", "t","Annotated.gene.s.")
dim(meth1)
methg <- meth1$`Annotated.gene.s.`
methg_spl <- strsplit(methg,";")
meth1 <- meth1[which(lapply(methg_spl,length)>0),]
head(meth1)
dim(meth1)
x <- apply(meth1,1,function(x) {
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
head(meth_fix2)
dim(meth_fix2)

meth1 <- meth_fix2[,c("gene","t")]
meth1 <- aggregate(. ~ gene, meth1, sum)
dim(meth1)
l <- list("methylation"=meth1)
m <- mitch_import(x = l,geneIDcol = "gene",DEtype = "limma")

download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip", destfile="ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
genesets <- gmt_import("ReactomePathways.gmt")

res <- mitch_calc(x=m,genesets = genesets ,priority = "significance",minsetsize = 5, resrows = 25)
head(res$enrichment_result,50)
mitch_report(res = res, outfile = "DNAm_significance_reactome.html")
mitch_plots(res,outfile="DNAm_significance_reactome.pdf")

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

write.table(res1, "results enrichment DNA methylation reactome.txt")

library("ggrepel")

tiff('pathway VO2max DNAmeth Reactome.tiff', width =8, height = 6, units = 'in', res=300)
pvsef=ggplot(res1, aes(x=`s.dist`, y=-log(`p.adjustANOVA`)))+
  geom_point(aes(col=sigdir))+
  scale_color_manual(values=c("#1BD2F7","gray","#1AEDD8"), name = " ")+
  theme_minimal()+
  labs(y="-log(p.adjustedANOVA)", x="Effect size")+
  geom_label_repel(data=filter(res1[1:10,]), aes(label=set), max.overlaps = 7, color=c("#1BD2F7","#1AEDD8", "#1AEDD8","#1BD2F7","#1BD2F7","#1BD2F7","#1AEDD8", "#1AEDD8",  "#1AEDD8","#1BD2F7" ))
pvsef
dev.off()

#Graph heatmap and circular plot
library(enrichplot)
library(mitch)


#Create an enrichResult object for the code to work
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures")
library(tidyverse)
reactome_enrichment <- res1%>%
  mutate(Description = str_extract(ID,
                                   "(?<=_)(.+)"))%>%
  mutate(Description = str_replace_all(Description,
                                       "_",
                                       " "))%>%
  mutate(Description = str_to_sentence(Description))

hpo_enrichment <- hpo_enrichment %>%
  mutate(SigGenesInSet = str_replace_all(SigGenesInSet,
                                         ",",
                                         "/"),
         qvalue=FDR,
         GeneRatio=paste(DE,N,sep="/"))%>%
  dplyr::rename(pvalue=P.DE,
                geneID=SigGenesInSet,
                `p.adjust`=FDR,
                Count=DE)%>%
  select(-c(N))
rownames(hpo_enrichment) <- hpo_enrichment$ID

library(qdapTools)
x <- new("enrichResult",
         result         = hpo_enrichment,
         pvalueCutoff   = 1,
         pAdjustMethod  = "BH",
         qvalueCutoff   = 0.005,
         gene           = unique(unlist(strsplit(hpo_enrichment$geneID,
                                                 ","))),
         universe       = list2df(hpo)$X1,
         geneSets       = hpo,
         organism       = "Homo sapiens",
         keytype        = "SYMBOL",
         ontology       = "HPO",
         readable       = FALSE
)


setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
tiff('GSEA_meta-analysis RE_HPO_circ.tiff',
     width =90,
     height = 130,
     units = 'mm',
     res=300)
cnetplot(x,
         showCategory=2,
         circular = TRUE,
         colorEdge = TRUE,
         cex_label_gene	= 0.5,
         show.legend = FALSE,
         cex_label_category	= 0.75) +
  theme(legend.position="none") 
dev.off()


#Cell type proportions as a dotplot
library(enrichplot)
library(mitch)
setwd("C:/Users/e5103562/OneDrive - Victoria University/")
genesets <- list.files(pattern = ".gmt",
                       recursive = TRUE)
c8 <- genesets[which(str_detect(genesets,
                                fixed("c8"))==TRUE)]
c8 <- gmt_import(c8[which(str_detect(c8,
                                     "entrez")==TRUE)])
c8 <- c8[which(str_detect(names(c8),
                          "SKELETAL_MUSCLE")==TRUE)]
c8 <- c8[which(str_detect(names(c8),
                          "FETAL",
                          negate=TRUE)==TRUE)]

#Create an enrichResult object for the code to work
setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
library(tidyverse)
c8_enrichment <- read_tsv("GSEA_meta-analysis RE_C8.txt")%>%
  mutate(Description = str_extract(ID,
                                   "(?<=_)(.+)"))%>%
  mutate(Description = str_replace_all(Description,
                                       "_",
                                       " "))%>%
  mutate(Description = str_replace_all(Description,
                                       "SKELETAL MUSCLE ",
                                       ""))%>%
  mutate(Description = str_to_sentence(Description))

c8_enrichment <- c8_enrichment %>%
  mutate(SigGenesInSet = str_replace_all(SigGenesInSet,
                                         ",",
                                         "/"),
         `Significance (-log10((p-value))` = -log10(P.DE),
         P.DE=1,
         qvalue=FDR,
         GeneRatio=paste(DE,N,sep="/"))%>%
  dplyr::rename(pvalue=P.DE,
                geneID=SigGenesInSet,
                `p.adjust`=FDR,
                Count=DE)%>%
  dplyr::select(-c(N))
rownames(c8_enrichment) <- c8_enrichment$ID

library(qdapTools)
library(clusterProfiler)
x <- new("enrichResult",
         result         = c8_enrichment,
         pvalueCutoff   = 1,
         pAdjustMethod  = "BH",
         qvalueCutoff   = 1,
         gene           = unique(unlist(strsplit(c8_enrichment$geneID,
                                                 ","))),
         universe       = list2df(c8)$X1,
         geneSets       = c8,
         organism       = "Homo sapiens",
         keytype        = "SYMBOL",
         ontology       = "c8",
         readable       = FALSE
)

setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
tiff('GSEA_meta-analysis RE_C8.tiff',
     width =100,
     height = 150,
     units = 'mm',
     res=600)
dotplot(x,
        color = "pvalue",
        x = "`Significance (-log10((p-value))`",
        showCategory=30)+
  xlim(0,4) +
  geom_vline(xintercept = -log10(0.00041),
             lty=2)+
  theme(legend.position="none") 
dev.off()


####Identify transcription factors for DNA methylation and fitness####
#Create a background with all CpGs
# +/- 200bp around the CpGs. Your foreground should be the regions around the CpGs in your DMRs and your background the regions around all the CpGs.

library(DMRcate)
library(IRanges)


#First extract the regions hypo and hyper methylated
regions_hypo <- makeGRangesFromDataFrame(DMRs %>% filter(`Maximum.effect.size.in.DMR`<0),
                                         seqnames.field="Chromosome",
                                         start.field="Position.start..hg38.",
                                         end.field="Position.end..hg38.")
regions_hyper <- makeGRangesFromDataFrame(DMRs %>% filter(`Maximum.effect.size.in.DMR`>0),
                                          seqnames.field="Chromosome",
                                          start.field="Position.start..hg38.",
                                          end.field="Position.end..hg38.")
cpgs <- GenomicRanges::GRanges(seqnames = annotation_overlap_only$CpG_chrm ,
                               ranges = IRanges::IRanges(start = annotation_overlap_only$CpG_beg,
                                                         end = annotation_overlap_only$CpG_end),
                               name = annotation_overlap_only$probeID)
overlaps_hypo <- GenomicRanges::findOverlaps(cpgs,
                                             regions_hypo)
overlaps_hyper <- GenomicRanges::findOverlaps(cpgs,
                                              regions_hyper)
DMR.cpg.hypo <- cpgs$name[from(overlaps_hypo)]
DMR.cpg.hyper <- cpgs$name[from(overlaps_hyper)]


#select the anno of all CpGs used in the analysis
all_cpgs=annotation[annotation$probeID %in% meta_res_robust$CpG,]

all_cpgs$CpG_beg200=c(all_cpgs$CpG_beg-200)
all_cpgs$CpG_end200=c(all_cpgs$CpG_end+200)
#turn this into a Granges object with the correct column names
userUniverse=makeGRangesFromDataFrame(all_cpgs,
                                      keep.extra.columns=FALSE,
                                      ignore.strand=FALSE,
                                      seqinfo=NULL,
                                      seqnames.field=c("CpG_chrm"),
                                      start.field="CpG_beg200",
                                      end.field="CpG_end200",
                                      strand.field="probe_strand",
                                      starts.in.df.are.0based=FALSE)

#convert to df
userUniverse_df2=as.data.frame(userUniverse)
#save as bed file
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures")
write.table(userUniverse_df2, file="background3.bed", quote=F, sep="\t", row.names=F, col.names=F)

#Do the same for significant ones now - Hypo
CpGs2 <- DMR.cpg.hypo
meta_res_robust2=meta_res_robust[meta_res_robust$CpG %in% CpGs2, ]


annotation_overlap_only=annotation[is.element(annotation$probeID,intersect(annotation$probeID,meta_res_robust2$CpG)),]
annotation_overlap_only=annotation_overlap_only[!duplicated(annotation_overlap_only$probeID),]

library(tidyverse)
meta_res_robust2=mutate(meta_res_robust2, t=meta_res_robust2$`Effect size`/meta_res_robust2$SE)
annotation_overlap_only_2=left_join(annotation_overlap_only, meta_res_robust2, by=join_by(probeID==CpG))

#We need to create a "CpGannotated" object to be used in dmrcate
annotated_only <- GRanges(as.character(annotation_overlap_only_2$CpG_chrm), #chromosome
                          IRanges(start=c(annotation_overlap_only_2$CpG_beg),end=c(annotation_overlap_only_2$CpG_end)), #position on chromosome
                          stat = annotation_overlap_only_2$t, #t-statistic
                          diff = annotation_overlap_only_2$`Effect size`, #effect size
                          ind.fdr = annotation_overlap_only_2$`FDR`, #adjusted p-value
                          is.sig = annotation_overlap_only_2$`FDR` < 0.05) #p-value threshold
names(annotated_only) <- annotation_overlap_only_2$probeID
annotated_only <- sort(annotated_only) # i think this is the object I put into the line to extract the DMRs 

#repeat for the CpGs in my DMRs (*must do CpGs in DMRS as opposed to DMR regions, wont work that way)
hypoDMRcpgs <- sapply(1:length(resultsRanges), function (x) names(subsetByOverlaps(annotated_only, resultsRanges[x]))) #run dmrcate, filter for DMRs<0.005
hypoDMRcpgs_unlisted=unlist(hypoDMRcpgs) #1115 hypo-cpgs in DMRs
hypoDMRcpgs_anno=annotation[annotation$probeID %in% hypoDMRcpgs_unlisted,]

hypoDMRcpgs_anno$CpG_beg200=c(hypoDMRcpgs_anno$CpG_beg-200)
hypoDMRcpgs_anno$CpG_end200=c(hypoDMRcpgs_anno$CpG_end+200)

hypoUserSets2=makeGRangesFromDataFrame(hypoDMRcpgs_anno,
                                       keep.extra.columns=FALSE,
                                       ignore.strand=FALSE,
                                       seqinfo=NULL,
                                       seqnames.field=c("CpG_chrm"),
                                       start.field="CpG_beg200",
                                       end.field="CpG_end200",
                                       strand.field="probe_strand",
                                       starts.in.df.are.0based=FALSE)
hypoUserSets_df2=as.data.frame(hypoUserSets2)
write.table(hypoUserSets_df2, file="hypoCpGsinDMRs2.bed", quote=F, sep="\t", row.names=F, col.names=F)

#Do the same for significant ones now - hyper
CpGs2 <- DMR.cpg.hyper
meta_res_robust2=meta_res_robust[meta_res_robust$CpG %in% CpGs2, ]


annotation_overlap_only=annotation[is.element(annotation$probeID,intersect(annotation$probeID,meta_res_robust2$CpG)),]
annotation_overlap_only=annotation_overlap_only[!duplicated(annotation_overlap_only$probeID),]

library(tidyverse)
meta_res_robust2=mutate(meta_res_robust2, t=meta_res_robust2$`Effect size`/meta_res_robust2$SE)
annotation_overlap_only_2=left_join(annotation_overlap_only, meta_res_robust2, by=join_by(probeID==CpG))


#We need to create a "CpGannotated" object to be used in dmrcate
annotated_only <- GRanges(as.character(annotation_overlap_only_2$CpG_chrm), #chromosome
                          IRanges(start=c(annotation_overlap_only_2$CpG_beg),end=c(annotation_overlap_only_2$CpG_end)), #position on chromosome
                          stat = annotation_overlap_only_2$t, #t-statistic
                          diff = annotation_overlap_only_2$`Effect size`, #effect size
                          ind.fdr = annotation_overlap_only_2$`FDR`, #adjusted p-value
                          is.sig = annotation_overlap_only_2$`FDR` < 0.05) #p-value threshold
names(annotated_only) <- annotation_overlap_only_2$probeID
annotated_only <- sort(annotated_only) # i think this is the object I put into the line to extract the DMRs 

#repeat for the CpGs in my DMRs (*must do CpGs in DMRS as opposed to DMR regions, wont work that way)
hyperDMRcpgs <- sapply(1:length(resultsRanges), function (x) names(subsetByOverlaps(annotated_only, resultsRanges[x]))) #run dmrcate, filter for DMRs<0.005
hyperDMRcpgs_unlisted=unlist(hyperDMRcpgs) #1115 hyper-cpgs in DMRs
hyperDMRcpgs_anno=annotation[annotation$probeID %in% hyperDMRcpgs_unlisted,]

hyperDMRcpgs_anno$CpG_beg200=c(hyperDMRcpgs_anno$CpG_beg-200)
hyperDMRcpgs_anno$CpG_end200=c(hyperDMRcpgs_anno$CpG_end+200)

hyperUserSets2=makeGRangesFromDataFrame(hyperDMRcpgs_anno,
                                        keep.extra.columns=FALSE,
                                        ignore.strand=FALSE,
                                        seqinfo=NULL,
                                        seqnames.field=c("CpG_chrm"),
                                        start.field="CpG_beg200",
                                        end.field="CpG_end200",
                                        strand.field="probe_strand",
                                        starts.in.df.are.0based=FALSE)
hyperUserSets_df2=as.data.frame(hyperUserSets2)
write.table(hyperUserSets_df2, file="hyperCpGsinDMRs2.bed", quote=F, sep="\t", row.names=F, col.names=F)

#Load TF results
library(tidyverse)

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures/UniBind_hypermethylation_VO2max")
TF_baseline_hyper=read.delim("allEnrichments.tsv")
TF_baseline_hyper$pvalue=10^(-TF_baseline_hyper$pValueLog)
TF_baseline_hyper$FDR=p.adjust(TF_baseline_hyper$pvalue)
TF_baseline_hyper_sig=filter(TF_baseline_hyper, TF_baseline_hyper$FDR<0.05)
unique(TF_baseline_hyper$collection)

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures/UniBind_hypomethylation_VO2max")
TF_baseline_hypo=read.delim("allEnrichments.tsv")
TF_baseline_hypo$pvalue=10^(-TF_baseline_hypo$pValueLog)
TF_baseline_hypo$FDR=p.adjust(TF_baseline_hypo$pvalue)
TF_baseline_hypo_sig=filter(TF_baseline_hypo, TF_baseline_hypo$FDR<0.05)
unique(TF_baseline_hypo_sig$collection)

####recreate swarm plot for TFs - hyper####
library(tidyverse)
library(ggbeeswarm)

#write.csv(TF_baseline_hyper_sig, "TF_baseline_hyper_sig.csv")
unique(TF_baseline_hyper_sig$collection)

TF_baseline_hyper = mutate(TF_baseline_hyper, sig=ifelse(TF_baseline_hyper$FDR <0.05, "FDR<0.05", "Not Sig"))
TF_baseline_hyper$collection2=TF_baseline_hyper$collection

TF_baseline_hyper$collection2=rep("NA")
TF_baseline_hyper$collection2[grep("SIX2", TF_baseline_hyper$collection)]="#99CCFF"
TF_baseline_hyper$collection2[grep("PPARG", TF_baseline_hyper$collection)]="#0066CC"
TF_baseline_hyper$collection2[grep("GRHL2", TF_baseline_hyper$collection)]="#99FF99"
TF_baseline_hyper$collection2[grep("NEUROG2", TF_baseline_hyper$collection)]="#009900"
TF_baseline_hyper$collection2[grep("ESR1", TF_baseline_hyper$collection)]="#CC99CC"

TF_baseline_hyper$collection2[grep("Not Sig", TF_baseline_hyper$sig)]="#CCCCCC"
TF_baseline_hyper$collection2[grep("NA", TF_baseline_hyper$collection2)]="#FFFFCC"

setwd("D:/Papers to work on/Methylation Meta-Analyses/New analyses - VO2max and aerobic only")
tiff('swarm baseline_hyper.tiff', width =8, height = 6, units = 'in', res=600)
beeswarm::beeswarm(TF_baseline_hyper$pValueLog,pch=19,
                   pwcol=TF_baseline_hyper$collection2)
legend("topright", legend=c("SIX2","PPARG","GRHL2","NEUROG2","ESR1","Other","Not Sig"), 
       col=c("#99CCFF","#0066CC","#99FF99", "#009900","#CC99CC", "#FFFFCC", "#CCCCCC"), pch=19)
dev.off()


####TFBS hyper-methylation####
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures")
TFBS_hyper=read.table("hypermeth_VO2max_TFBS.bed")

TFBS_hyper2=TFBS_hyper%>%
  dplyr::filter(V5 %in% c("SIX2" , "PPARG" , "GRHL2" , "NEUROG2" , "ESR1"))%>%
  dplyr::select("chromossome_name"=V1, "start"=V2, "end"=V3, "TF"=V5, "strain"=V8)%>%
  unique()

library("stringr")
TFBS_hyper3=str_replace_all(TFBS_hyper2$chromossome_name,"chr","" )
TFBS_hyper2$chromossome_name=TFBS_hyper3

library("biomaRt")
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")

all.genes <- getBM(
  attributes=c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position"),
  filters=c("chromosome_name", "start", "end"),
  values=list(TFBS_hyper2$chromossome_name , TFBS_hyper2$start ,TFBS_hyper2$end),
  mart=mart)


write_csv(all.genes, "TFBS for hypermethylation and VO2max.csv")

#59508 TFBS encountered for the 5 significant TFs

all.genes$hgnc_symbol[all.genes$hgnc_symbol==""]<-NA

all.genes2<-all.genes%>%
  drop_na()

uniquegenes=unique(all.genes2$hgnc_symbol)

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Transcriptome/Tables and Figures")
trans=read.delim("METAL_muscle_VO2max_Sarah.txt")

mRNA <- trans%>%
  dplyr::filter(FDR<0.05)

merged_hyper=inner_join(mRNA, all.genes2, by=join_by("Gene.name"=="hgnc_symbol"))

merged_hyper2=left_join(mRNA, all.genes2, by=join_by("Gene.name"=="hgnc_symbol"))

merged_hyper2=merged_hyper2%>%
  unite("location", "chromosome_name","start_position", sep="_")%>%
  unite("location", "location","end_position", sep="_")

TFBS_hyper2=TFBS_hyper2%>%
  unite("location", "chromossome_name","start", sep="_")%>%
  unite("location", "location","end", sep="_")

merged_hyper2=left_join(merged_hyper2, TFBS_hyper2, by="location")
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Transcriptome/Tables and Figures")
write_csv(merged_hyper2, "TFBS for HyperDMRS genes.csv")

####Check how likely intersected genes are simply by chance
#hypergeo test for enrichment
CpGs_DMRs=unique(c(DMR.cpg.hyper, DMR.cpg.hypo))

library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

## Ensembl IDs of interest
genenames <- mRNA$Gene.name
## Run biomaRt query
entrezs <- getBM(filters = "external_gene_name",
                 attributes = c("external_gene_name", "entrezgene_id"),
                 values = genenames,
                 mart = mart)

library("missMethyl")
gsa_hypergeo_enrich_chapman=gsameth(sig.cpg=CpGs_DMRs,
                                    all.cpg = meta_res_robust$CpG,
                                    collection = entrezs$entrezgene_id,
                                    array.type = "EPIC",
                                    plot.bias = TRUE,
                                    prior.prob = TRUE,
                                    anno = NULL,
                                    equiv.cpg = TRUE,
                                    fract.counts = TRUE,
                                    sig.genes = TRUE) 
gsa_hypergeo_enrich_chapman


SigGenesInSet=unlist(strsplit(gsa_hypergeo_enrich_chapman[1,"SigGenesInSet"], ","))

#FDR=0.00031


####TFBS hypo-methylation####
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures")
TFBS_hypo=read.table("hypometh_VO2max_TFBS.bed")

hypoTF=unique(TF_baseline_hypo_sig$collection)

TFBS_hypo2=TFBS_hypo%>%
  dplyr::filter(V5 %in% c("FLI1","ETS1","ERG","RELA","NR2F2","RUNX1","SPI1","GATA2","NR3C1","ELF1","MYC","MYCN","JUN","ETV6","BCL6"))%>%
  dplyr::select("chromossome_name"=V1, "start"=V2, "end"=V3, "TF"=V5, "strain"=V8)%>%
  unique()

library("stringr")
TFBS_hypo3=str_replace_all(TFBS_hypo2$chromossome_name,"chr","" )
TFBS_hypo2$chromossome_name=TFBS_hypo3

#3900 TFBS encountered for the 15 significant TFs

library("biomaRt")
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")

all.genes_hypo <- getBM(
  attributes=c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position"),
  filters=c("chromosome_name", "start", "end"),
  values=list(TFBS_hypo2$chromossome_name , TFBS_hypo2$start ,TFBS_hypo2$end),
  mart=mart)


write_csv(all.genes_hypo, "TFBS for hypomethylation and VO2max.csv")

#38409 TFBS encountered for the 15 significant TFs

all.genes_hypo$hgnc_symbol[all.genes_hypo$hgnc_symbol==""]<-NA

all.genes2_hypo<-all.genes_hypo%>%
  drop_na()

uniquegenes=unique(all.genes2_hypo$hgnc_symbol)

merged_hypo=inner_join(mRNA, all.genes2_hypo, by=join_by("Gene.name"=="hgnc_symbol"))

merged_hypo2=left_join(mRNA, all.genes2_hypo, by=join_by("Gene.name"=="hgnc_symbol"))

