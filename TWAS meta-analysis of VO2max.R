#Set working directory
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Transcriptome")
directory=c("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Transcriptome")
#List tbl files
files <- list.files(pattern=".tbl",
           recursive = TRUE)

#Add gene codes to Gene SMART data
setwd("E:/Macsue - Hard drive/PhD Victoria University/Papers to work on/Methylation Meta-Analyses/Transcriptomics")
trans=read_csv("Metanalyses Transcriptomics.csv")

Genes=data.frame(Genes=trans$ID, ENSEMBL=trans$MarkerName)

setwd(directory)
GS=read.delim("GeneSMART_mRNA_VO2max_corrected.tbl")
GS=left_join(GS, Genes, by="ENSEMBL")

write.table(GS, "GeneSMART_mRNA_VO2max_corrected.tbl",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

####Bacon process####
library(bacon)
library(tidyverse)

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
files <- list.files(pattern="corrected.tbl",
                    recursive = TRUE)
file <- "METAL_commands.txt"
write.table(paste("SCHEME STDERR",
                  "MARKER GENE",
                  "EFFECT EFFECTSIZE_CORR",
                  "SEPARATOR TAB",
                  "STDERR SE_CORR",
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



####Load METAL results####
library(tidyverse)
meta_res <- read_tsv("Y:/Ageing-and-exercise/TWAS meta-analysis of VO2max/METAL_muscle_VO2max.TBL")%>%
    dplyr::rename(`Gene name` = MarkerName)
max(meta_res$HetDf)

#Calculate nb of genes included in analysis for each nb of studies
nb_genes <- c()
for (i in 1:max(meta_res$HetDf))
{
    nb_genes <- c(nb_genes,nrow(meta_res %>% dplyr::filter(HetDf>=(i-1))))
}

cumsum <- tibble(`Number of included studies` = 1:max(meta_res$HetDf),
                 `Number of genes present in all studies` = nb_genes)

ggplot(cumsum,
       mapping = aes(x = `Number of included studies`,
                     y = `Number of genes present in all studies`))+
    geom_point()+
    geom_line()+
    theme_classic()

#View results including at least 2 studies
meta_res_robust <- meta_res %>% dplyr::filter(HetDf>=1)
padj <- p.adjust(meta_res_robust$`P-value`,method="BH")
meta_res_robust <- meta_res_robust %>%
    mutate(`Adjusted p-value` = padj)
meta_res_robust <- meta_res_robust %>%
    mutate(sig = ifelse(`Adjusted p-value`<0.05,
                        "Sig",
                        "Not Sig"))
nrow(meta_res_robust) #17405

#####Arrange
sigdir <- meta_res_robust$sig
sigdir[meta_res_robust$sig=="Sig"&
           meta_res_robust$Effect<0]="Under"
sigdir[meta_res_robust$sig=="Sig"&
           meta_res_robust$Effect>0]="Over"
meta_res_robust <- meta_res_robust %>%
    mutate(sigdir = sigdir)
meta_res_robust <- meta_res_robust %>%
    arrange(desc(sigdir))

#Select relevant columns only and rename them
meta_res_robust <- meta_res_robust %>%
    dplyr::select(`Gene name`,
                  Effect,
                  StdErr,
                  `P-value`,
                  `Adjusted p-value`,
                  HetDf,
                  HetISq,
                  HetPVal,
                  sig,
                  sigdir)%>%
    dplyr::rename(`log2FC per unit of VO2max` = Effect,
                  SE = StdErr,
                  FDR = `Adjusted p-value`,
                  `Number of studies` = HetDf,
                  `Heterogeneity index (I2)` = HetISq,
                  `Heterogeneity p-value` = HetPVal,
                  Significance = sig,
                  Direction = sigdir)

write.table(meta_res_robust,
            file="METAL_muscle_VO2max_Sarah.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")


#### Tables and Figures ####
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Transcriptome/Tables and Figures")
library(tidyverse)
#Load METAL results that have been formatted to annotate each CpG

to_write <- meta_res_robust

to_write <- to_write %>%
  dplyr::rename(`P-value VO2max` = `P-value`,
                `FDR VO2max` = FDR)%>%
  arrange(`P-value VO2max`)

DEGs <- subset(to_write,
               `FDR VO2max`<0.05) #162

#Volcano plot for EWAS RE meta-analysis of VO2max

tiff('Volcano plot_TWAS RE meta-analysis of VO2max_FINAL.tiff',
     width =100,
     height = 80,
     units = 'mm',
     res=300)
ggplot() +
  geom_point(data = to_write,
             aes(`log2FC per unit of VO2max`,
                 -log10(`P-value VO2max`),
                 col=Direction),
             size=0.5)+
  scale_color_manual(values=c("grey","#1BD2F7", "#1AEDD8"))+
  labs(y="-log10(p-value)")+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Transcriptome/Tables and Figures")

#Show histogram of effect sizes for DMPs
tiff('Distribution of effect size in mRNA_FINAL.tiff',
     width =100,
     height = 65,
     units = 'mm',
     res=300)
ggplot(data = DEGs,
       aes(x=`log2FC per unit of VO2max`,
           fill=Direction)) +
  geom_histogram(colour = "grey")+
  scale_fill_manual(values=c("#1BD2F7","#1AEDD8"))+
  lims(x=c(min(to_write$`log2FC per unit of VO2max`),
           max(to_write$`log2FC per unit of VO2max`)))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

####Forest plot for top hypo and hyper DMP####
####Prep####
library(readxl)
dataset_summary <- read_excel('//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Transcriptome/Tables and Figures/TWAS datasets.xlsx',
                              na = c("","NA"))%>%
  dplyr::rename(prop_males = `% Male`)
dataset_summary=dataset_summary[c(1:4,6),]

#Create empty list
L.names <- meta_res_robust$`Gene name`
L <- vector("list", length(L.names))
names(L) <- L.names

#Initiate list
tib <- tibble(Genes = meta_res_robust$`Gene name` ,
              Dataset = rep("Meta-analysis",nrow(meta_res_robust)),
              ES = meta_res_robust$`log2FC per unit of VO2max`,
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
names(L) <- meta_res_robust$`Gene name`

dataset_summary$`Dataset ID`=c("GSE18732","HERITAGE","GSE44818", "The Malmö Prevention Study","MSAT")
colnames(dataset_summary)=c("Dataset", "Muscle type", "n","Age (mean ± SD)", "prop_males")

#List tbl files
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Transcriptome/Tables and Figures")
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
  summary <- tibble(Genes = file$GENE,
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

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Transcriptome/Tables and Figures/")

saveRDS(L,
        file = "METAL_muscle_VO2max.rds")

####Forestplot - start here####
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Transcriptome/Tables and Figures")

L <- read_rds("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Transcriptome/Tables and Figures/METAL_muscle_VO2max.rds")

#1
top_down1 <- DEGs %>%
  filter(`log2FC per unit of VO2max` < 0) %>%
  dplyr::slice(1)%>%
  pull(`Gene name`) #RPUSD4
top_down1=L[[top_down1]]
top_down1$Dataset=sub("\\..*","",top_down1$Dataset)
top_down1[4,2]="MSAT"

tib_down1 <- left_join(top_down1,
                       dataset_summary)


tib_down1=tib_down1 %>%
  dplyr::select("Gene name","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")%>%
  rename("log2FC"="ES")

#2
top_down2 <- DEGs %>%
  filter(`log2FC per unit of VO2max` < 0) %>%
  dplyr::slice(2)%>%
  pull(`Gene name`) #SH3RF2
top_down2=L[[top_down2]]
top_down2$Dataset=sub("\\..*","",top_down2$Dataset)
top_down2[4,2]="MSAT"

tib_down2 <- left_join(top_down2,
                       dataset_summary)

tib_down2=tib_down2 %>%
  dplyr::select("Gene name","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")%>%
  rename("log2FC"="ES")


#3
top_down3 <- DEGs %>%
  filter(`log2FC per unit of VO2max` < 0) %>%
  dplyr::slice(3)%>%
  pull(`Gene name`) #PRKAG3
top_down3=L[[top_down3]]
top_down3$Dataset=sub("\\..*","",top_down3$Dataset)
top_down3[4,2]="MSAT"

tib_down3 <- left_join(top_down3,
                       dataset_summary)%>%
  arrange(-`n`)

tib_down3=tib_down3 %>%
  dplyr::select("Gene name","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")%>%
  rename("log2FC"="ES")

#4
top_down4 <- DEGs %>%
  filter(`log2FC per unit of VO2max` < 0) %>%
  dplyr::slice(4)%>%
  pull(`Gene name`) #PILRB
top_down4=L[[top_down4]]
top_down4$Dataset=sub("\\..*","",top_down4$Dataset)
top_down4[4,2]="MSAT"
top_down4[2,2]="The Malmö Prevention Study"

tib_down4 <- left_join(top_down4,
                       dataset_summary)%>%
  arrange(-`n`)

tib_down4=tib_down4 %>%
  dplyr::select("Gene name","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")%>%
  rename("log2FC"="ES")

#5
top_down5 <- DEGs %>%
  filter(`log2FC per unit of VO2max` < 0) %>%
  dplyr::slice(5)%>%
  pull(`Gene name`) #NOXA1
top_down5=L[[top_down5]]
top_down5$Dataset=sub("\\..*","",top_down5$Dataset)
top_down5[3,2]="MSAT"

tib_down5 <- left_join(top_down5,
                       dataset_summary)%>%
  arrange(-`n`)

tib_down5=tib_down5 %>%
  dplyr::select("Gene name","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")%>%
  rename("log2FC"="ES")


top5_down=rbind(tib_down1,tib_down2,tib_down3,tib_down4,tib_down5)

top5_down<-top5_down%>%
  mutate(Dataset=factor(Dataset, levels= c("Meta-analysis", "GSE18732", "HERITAGE","GSE44818","The Malmö Prevention Study", "MSAT")))%>%
  rename("Gene"="Gene name")

library(ggforestplot)

tiff("Top5 down mRNAs.tiff",
     height = 11,
     width = 10,
     unit = 'in',
     res = 200)
ggforestplot::forestplot(
  df = top5_down,
  name = Gene,
  estimate = log2FC,
  se = SE,
  pvalue = FDR,
  psignif = 0.05,
  colour = Dataset,
  shape = Dataset,
  #xlab = "Odds ratio for increased VO2max (95% CI) per 1−SD decrease in log2FC mRNA expression",
  title = "Top 5 down-expressed mRNAs associated with VO2max",
  logodds = F,
  size=3
)+
  ggplot2::scale_shape_manual(
    values = c(23L, 21L, 21L, 21L, 21L,21L,21L,21L,21L,21L),
    labels = c("Meta-analysis", "GSE18732", "HERITAGE","GSE44818","The Malmö Prevention Study","MSAT")
  )+
  ggplot2::scale_colour_manual(values = c("#FF9A97","#ACFF97","#1AEDD8","#1BD2F7", "#9A97FF","#DB90FF"))
dev.off()


#1
top_up1 <- DEGs %>%
  filter(`log2FC per unit of VO2max` >0) %>%
  dplyr::slice(1)%>%
  pull(`Gene name`) #LPL
top_up1=L[[top_up1]]
top_up1$Dataset=sub("\\..*","",top_up1$Dataset)
top_up1[5,2]="MSAT"
top_up1[2,2]="The Malmö Prevention Study"

tib_up1 <- left_join(top_up1,
                        dataset_summary)%>%
  arrange(-`n`)

tib_up1=tib_up1 %>%
  dplyr::select("Gene name","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")%>%
  rename("log2FC"="ES")

#2
top_up2 <- DEGs %>%
  filter(`log2FC per unit of VO2max` >0) %>%
  dplyr::slice(2)%>%
  pull(`Gene name`) #KIFBP
top_up2=L[[top_up2]]
top_up2$Dataset=sub("\\..*","",top_up2$Dataset)
top_up2[4,2]="MSAT"
top_up2[2,2]="The Malmö Prevention Study"

tib_up2 <- left_join(top_up2,
                        dataset_summary)%>%
  arrange(-`n`)

tib_up2=tib_up2 %>%
  dplyr::select("Gene name","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")%>%
  rename("log2FC"="ES")

#3
top_up3 <- DEGs %>%
  filter(`log2FC per unit of VO2max` >0) %>%
  dplyr::slice(3)%>%
  pull(`Gene name`) #ADGRF5
top_up3=L[[top_up3]]
top_up3$Dataset=sub("\\..*","",top_up3$Dataset)
top_up3[4,2]="MSAT"
top_up3[2,2]="The Malmö Prevention Study"
tib_up3 <- left_join(top_up3,
                        dataset_summary)%>%
  arrange(-`n`)

tib_up3=tib_up3 %>%
  dplyr::select("Gene name","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")%>%
  rename("log2FC"="ES")

#4
top_up4 <- DEGs %>%
  filter(`log2FC per unit of VO2max` >0) %>%
  dplyr::slice(4)%>%
  pull(`Gene name`) #NCEH1
top_up4=L[[top_up4]]
top_up4$Dataset=sub("\\..*","",top_up4$Dataset)
top_up4[4,2]="MSAT"

tib_up4 <- left_join(top_up4,
                        dataset_summary)%>%
  arrange(-`n`)

tib_up4=tib_up4 %>%
  dplyr::select("Gene name","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")%>%
  rename("log2FC"="ES")

#5
top_up5 <- DEGs %>%
  filter(`log2FC per unit of VO2max` >0) %>%
  dplyr::slice(5)%>%
  pull(`Gene name`) #LDHB
top_up5=L[[top_up5]]
top_up5$Dataset=sub("\\..*","",top_up5$Dataset)
top_up5[5,2]="MSAT"
top_up5[2,2]="The Malmö Prevention Study"

tib_up5 <- left_join(top_up5,
                        dataset_summary)%>%
  arrange(-`n`)

tib_up5=tib_up5 %>%
  dplyr::select("Gene name","Dataset", "ES", "SE", "PVAL","FDR","n","prop_males","Type")%>%
  rename("log2FC"="ES")

top5_up=rbind(tib_up1,tib_up2,tib_up3,tib_up4,tib_up5)

top5_up<-top5_up%>%
  mutate(Dataset=factor(Dataset, levels= c("Meta-analysis", "GSE18732", "GSE44818","The Malmö Prevention Study","MSAT")))

fullTopCPGs=rbind(top5_up, top5_down)
fullTopCPGs<-fullTopCPGs%>%
  mutate(Dataset=factor(Dataset, levels= c("Meta-analysis", "GSE18732", "HERITAGE","GSE44818","The Malmö Prevention Study","MSAT")))%>%
  rename("Genes"="Gene name")

tiff("Top5 mRNAs.tiff",
     height = 11,
     width = 10,
     unit = 'in',
     res = 200)
ggforestplot::forestplot(
  df = fullTopCPGs,
  name = Genes,
  estimate = log2FC,
  se = SE,
  pvalue = FDR,
  psignif = 0.05,
  colour = Dataset,
  shape = Dataset,
  #xlab = "Odds ratio per change in VO2max (95% CI) per 1−SD decrease in log2FC mRNA expression",
  title = "Top 5 down- and up- expressed mRNA associated with VO2max",
  logodds = F
)+
  ggplot2::scale_shape_manual(
    values = c(23L, 21L, 21L, 21L, 21L,21L,21L,21L,21L,21L),
    labels = c("Meta-analysis", "GSE18732", "GSE44818", "The Malmö Prevention Study","MSAT")
  )+
  ggplot2::scale_colour_manual(values = c("#FF9A97","#1AEDD8","#1BD2F7", "#9A97FF","#DB90FF", 
                                           "#987C89"))
dev.off()

####pathway analysis####
library("mitch")
library(tidyverse)

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Transcriptome/Tables and Figures")
trans=read.delim("METAL_muscle_VO2max_Sarah.txt")

Meta_Trans_mitch=trans%>%
  mutate(t=log2FC.per.unit.of.VO2max/SE)%>%
  select("log2FC.per.unit.of.VO2max", "SE", t, P.value, FDR, Gene.name)%>%
  drop_na()
rownames(Meta_Trans_mitch)=Meta_Trans_mitch$Gene.name

trans <- Meta_Trans_mitch
colnames(trans)[6] <- "gene"
head(trans)
str(trans)

download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip", destfile="ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
genesets <- gmt_import("ReactomePathways.gmt")

trans <- trans[,c("gene","t")]
dim(trans)

l <- list("rna"=trans)
m <- mitch_import(x = l,geneIDcol = "gene",DEtype = "limma")


res <- mitch_calc(x=m,genesets = genesets ,priority = "significance",minsetsize = 5, resrows = 25)
head(res$enrichment_result,50)
mitch_report(res = res, outfile = "trans_significance_reactome.html")
mitch_plots(res,outfile="trans_significance_reactome.pdf")

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

write.table(res1, "results enrichment transcriptome reactome.txt")

library("ggrepel")

tiff('pathway VO2max transcriptome Reactome.tiff', width =10, height = 12, units = 'in', res=300)
pvsef=ggplot(res1, aes(x=`s.dist`, y=-log(`p.adjustANOVA`)))+
  geom_point(aes(col=sigdir))+
  scale_color_manual(values=c("#1BD2F7","gray","#1AEDD8"), name = " ")+
  theme_minimal()+
  labs(y="-log(p.adjustedANOVA)", x="s score")+
  geom_label_repel(data=filter(res1[1:10,]), aes(label=set), max.overlaps = 7, color=c("#1BD2F7","#1AEDD8", "#1BD2F7","#1BD2F7","#1BD2F7","#1AEDD8","#1BD2F7", "#1AEDD8",  "#1AEDD8","#1BD2F7" ))
pvsef
dev.off()

