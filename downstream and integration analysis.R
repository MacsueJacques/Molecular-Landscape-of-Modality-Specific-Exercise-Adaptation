#### Integration ####
library(tidyverse)

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)/Tables and Figures")
meth=read.delim("METAL_muscle_VO2max_withoutGSv2.txt")

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Transcriptome/Tables and Figures")
trans=read.delim("METAL_muscle_VO2max_Sarah.txt")

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Proteome/Metal Proteins/Tables and Figures")
prot=read_tsv("METAANALYSIS1.TBL")

prot <- prot %>% 
  mutate(FDR=p.adjust(prot$`P-value`, method = "fdr"))%>%
  mutate(diffexpressed = case_when(
    Effect > 0 & FDR < 0.05 ~ 'UP',
    Effect < 0 & FDR < 0.05 ~ 'DOWN',
    FDR > 0.05 ~ 'Not Sig'
  ))

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

#####################Methylation
meth1 <- meth%>%
  mutate(t=Effect.size/SE)%>%
  dplyr::select("CpG","Effect.size","SE","P.value","FDR", "t","Direction","Annotated.gene.s.")
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
meth_fix2$Direction <- unlist(meth_fix[,"Direction"])
head(meth_fix2)
dim(meth_fix2)
meth1 <- meth_fix2[,c("gene","t")]
meth1 <- aggregate(. ~ gene, meth1, sum)
dim(meth1)

# Get the genes that are present in your dataframe
genes_in_meth <- meth_fix2$gene
# Read in the .gmt file
pwl2 <- read.gmt("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/ReactomePathways.gmt") 
# Subset to the genes that are present in our dataset
pwl2 <- pwl2[pwl2$gene %in% genes_in_meth,] 
# Save the filtered background gene set
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax")
filename <- 'meth_path.RDS'
saveRDS(pwl2, filename)

# Remove non-significant genes
meth_path <- meth_fix2[meth_fix2$Direction != 'Not Sig', ]
# Substitute names so they are annotated nicely in the heatmap later
# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
meth_results_list <- split(meth_path, meth_path$Direction)

## Run ClusterProfiler -----------------------------------------------
# Settings
name_of_comparison <- 'VO2max' # for our filename
background_genes <- 'reactome' # for our filename
bg_genes <- readRDS("meth_path.RDS") # read in the background genes
padj_cutoff <- 0.05 # p-adjusted threshold, used to filter out pathways
genecount_cutoff <- 5 # minimum number of genes in the pathway, used to filter out pathways

res <- lapply(names(meth_results_list),
              function(x) enricher(gene = meth_results_list[[x]]$gene,
                                   TERM2GENE = bg_genes))
names(res) <- names(meth_results_list)
#Convert the list of enrichResults for each sample_pattern to a dataframe with the pathways
res_df <- lapply(names(res), function(x) rbind(res[[x]]@result))
names(res_df) <- names(res)
res_df <- do.call(rbind, res_df)
head(res_df)

#Convert the enrichResults to a dataframe with the pathways
res_df <- lapply(names(res), function(x) rbind(res[[x]]@result))
names(res_df) <- names(res)
res_df <- do.call(rbind, res_df)
head(res_df)

res_df <- res_df %>% mutate(minuslog10padj = -log10(p.adjust),
                            diffexpressed = gsub('\\.REACTOME.*$', '', rownames(res_df)))

# Subset to those pathways that have p adj < cutoff and gene count > cutoff (you can also do this in the enricher function)
target_pws <- unique(res_df$ID[res_df$p.adjust < padj_cutoff & res_df$Count > genecount_cutoff]) # select only target pathways have p adjusted < 0.05 and at least 6 genes
res_df <- res_df[res_df$ID %in% target_pws, ]

write.csv(res_df, 'meth_resclusterp.csv', row.names = FALSE)

res_df_up <- res_df %>% filter(grepl("Hyper", diffexpressed)) %>% 
  dplyr::select(!c('minuslog10padj', 'diffexpressed')) 
rownames(res_df_up) <- res_df_up$ID

# For visualisation purposes, let's shorten the pathway names
res_df_up$Description <- gsub('(H|h)iv', 'HIV', 
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
                                                                                                                gsub('REACTOME_', '', res_df_up$Description)))))))))))))))))))

enrichres_up <- new("enrichResult",
                    readable = FALSE,
                    result = res_df_up,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.2,
                    organism = "human",
                    ontology = "UNKNOWN",
                    gene = meth_fix2$gene,
                    keytype = "UNKNOWN",
                    universe = unique(bg_genes$gene),
                    gene2Symbol = character(0),
                    geneSets = bg_genes)
class(enrichres_down)


# Dotplot
tiff('dot plot hypo pathways.tiff',
     width =200,
     height = 300,
     units = 'mm',
     res=300)
dotplot(enrichres_down, showCategory = 20) + ggtitle("Top pathways from hypo methylated DMRs")
dev.off()

tiff('dot plot hyper pathways.tiff',
     width =200,
     height = 300,
     units = 'mm',
     res=300)
dotplot(enrichres_up, showCategory = 20) + ggtitle("Top pathways from hyper methylated DMRs")
dev.off()


####################Transcriptome
# Get the genes that are present in your dataframe
genes_in_trans <- trans$Gene.name
pwl2 <- read.gmt("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/ReactomePathways.gmt") 
# Subset to the genes that are present in our dataset
pwl2 <- pwl2[pwl2$gene %in% genes_in_trans,] 
# Subset to the genes that are present in our dataset
pwl2 <- pwl2[pwl2$gene %in% genes_in_trans,] 
# Save the filtered background gene set
filename <- 'trans_path.RDS'
saveRDS(pwl2, filename)

# Remove non-significant genes
trans <- trans[trans$Direction != 'Not Sig', ]

# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
trans_results_list <- split(trans, trans$Direction)

## Run ClusterProfiler -----------------------------------------------
# Settings
name_of_comparison <- 'VO2max' # for our filename
background_genes <- 'reactome' # for our filename
bg_genes <- readRDS("trans_path.RDS") # read in the background genes
padj_cutoff <- 0.05 # p-adjusted threshold, used to filter out pathways
genecount_cutoff <- 5 # minimum number of genes in the pathway, used to filter out pathways

res <- lapply(names(trans_results_list),
              function(x) enricher(gene = trans_results_list[[x]]$Gene.name,
                                   TERM2GENE = bg_genes))
names(res) <- names(trans_results_list)
#Convert the list of enrichResults for each sample_pattern to a dataframe with the pathways
res_df <- lapply(names(res), function(x) rbind(res[[x]]@result))
names(res_df) <- names(res)
res_df <- do.call(rbind, res_df)
head(res_df)

#Convert the enrichResults to a dataframe with the pathways
res_df <- lapply(names(res), function(x) rbind(res[[x]]@result))
names(res_df) <- names(res)
res_df <- do.call(rbind, res_df)
head(res_df)

res_df <- res_df %>% mutate(minuslog10padj = -log10(p.adjust),
                            diffexpressed = gsub('\\.REACTOME.*$', '', rownames(res_df)))

# Subset to those pathways that have p adj < cutoff and gene count > cutoff (you can also do this in the enricher function)
target_pws <- unique(res_df$ID[res_df$p.adjust < padj_cutoff]) # select only target pathways have p adjusted < 0.05 and at least 6 genes
res_df <- res_df[res_df$ID %in% target_pws, ]

write.csv(res_df, 'trans_resclusterp.csv', row.names = FALSE)

res_df_up <- res_df %>% filter(grepl("Over", diffexpressed)) %>% 
  dplyr::select(!c('minuslog10padj', 'diffexpressed')) 
rownames(res_df_up) <- res_df_up$ID

# For visualisation purposes, let's shorten the pathway names
res_df_up$Description <- gsub('(H|h)iv', 'HIV', 
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
                                                                                                                  gsub('REACTOME_', '', res_df_up$Description)))))))))))))))))))

enrichres_up <- new("enrichResult",
                    readable = FALSE,
                    result = res_df_up,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.2,
                    organism = "human",
                    ontology = "UNKNOWN",
                    gene = trans$Gene.name,
                    keytype = "UNKNOWN",
                    universe = unique(bg_genes$gene),
                    gene2Symbol = character(0),
                    geneSets = bg_genes)
class(enrichres_down)


# Dotplot
tiff('dot plot hypo pathways.tiff',
     width =200,
     height = 300,
     units = 'mm',
     res=300)
dotplot(enrichres_down, showCategory = 4) + ggtitle("Top pathways from under expressed mRNAs")
dev.off()

tiff('dot plot Over mRNAs.tiff',
     width =300,
     height = 200,
     units = 'mm',
     res=300)
dotplot(enrichres_up, showCategory = 6) + ggtitle("Top pathways from over expressed mRNAs")
dev.off()


#Proteome
# Get the genes that are present in your dataframe
genes_in_prot <- prot$MarkerName
pwl2 <- read.gmt("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/ReactomePathways.gmt") 
# Subset to the genes that are present in our dataset
pwl2 <- pwl2[pwl2$gene %in% genes_in_prot,] 
# Subset to the genes that are present in our dataset
pwl2 <- pwl2[pwl2$gene %in% genes_in_prot,] 
# Save the filtered background gene set
filename <- 'prot_path.RDS'
saveRDS(pwl2, filename)

# Remove non-significant genes
prot <- prot[prot$diffexpressed != 'Not Sig', ]

# Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
prot_results_list <- split(prot, prot$diffexpressed)

## Run ClusterProfiler -----------------------------------------------
# Settings
name_of_comparison <- 'VO2max' # for our filename
background_genes <- 'reactome' # for our filename
bg_genes <- readRDS("prot_path.RDS") # read in the background genes
padj_cutoff <- 0.05 # p-adjusted threshold, used to filter out pathways
genecount_cutoff <- 5 # minimum number of genes in the pathway, used to filter out pathways

res <- lapply(names(prot_results_list),
              function(x) enricher(gene = prot_results_list[[x]]$MarkerName,
                                   TERM2GENE = bg_genes))
names(res) <- names(prot_results_list)
#Convert the list of enrichResults for each sample_pattern to a dataframe with the pathways
res_df <- lapply(names(res), function(x) rbind(res[[x]]@result))
names(res_df) <- names(res)
res_df <- do.call(rbind, res_df)
head(res_df)

res_df <- res_df %>% mutate(minuslog10padj = -log10(p.adjust),
                            diffexpressed = gsub('\\.REACTOME.*$', '', rownames(res_df)))

# Subset to those pathways that have p adj < cutoff and gene count > cutoff (you can also do this in the enricher function)
target_pws <- unique(res_df$ID[res_df$p.adjust < padj_cutoff]) # select only target pathways have p adjusted < 0.05 and at least 6 genes
res_df <- res_df[res_df$ID %in% target_pws, ]

write.csv(res_df, 'prot_resclusterp.csv', row.names = FALSE)

res_df_up <- res_df %>% filter(grepl("UP", diffexpressed)) %>% 
  dplyr::select(!c('minuslog10padj', 'diffexpressed')) 
rownames(res_df_up) <- res_df_up$ID

# For visualisation purposes, let's shorten the pathway names
res_df_up$Description <- gsub('(H|h)iv', 'HIV', 
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
                                                                                                                  gsub('REACTOME_', '', res_df_up$Description)))))))))))))))))))

enrichres_up <- new("enrichResult",
                      readable = FALSE,
                      result = res_df_up,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.2,
                      organism = "human",
                      ontology = "UNKNOWN",
                      gene = prot$MarkerName,
                      keytype = "UNKNOWN",
                      universe = unique(bg_genes$gene),
                      gene2Symbol = character(0),
                      geneSets = bg_genes)
class(enrichres_down)


# Dotplot
tiff('dot plot Under Proteins.tiff',
     width =200,
     height = 300,
     units = 'mm',
     res=300)
dotplot(enrichres_down, showCategory = 4) + ggtitle("Top pathways from under expressed proteins")
dev.off()

tiff('dot plot Over Proteins.tiff',
     width =300,
     height = 270,
     units = 'mm',
     res=300)
dotplot(enrichres_up, showCategory = 20) + ggtitle("Top pathways from over expressed Proteins")
dev.off()



#ven diagram
library(tidyverse)
mRNA <- trans%>%
  filter(FDR<0.05)%>%
  dplyr::select("Gene.name" ,
                `log2FC.per.unit.of.VO2max`)%>%
  dplyr::rename(Gene="Gene.name")

nrow(trans %>% filter(Direction=="Under"))/162


DNAm <- meth %>%
  filter(FDR<0.05)%>%
  dplyr::rename(Gene = `Annotated.gene.s.`)%>%
  separate_rows(Gene,
                sep=";")%>%
  dplyr::select(Gene,
                `Effect.size`)
Prot=prot%>%
  mutate(FDR=p.adjust(prot$`P-value`, method = "fdr"))%>%
  filter(FDR<0.05)%>%
  dplyr::select("MarkerName", "Effect")%>%
  dplyr::rename(Gene="MarkerName", EFFECTSIZE="Effect")

merged <- left_join(DNAm,
                    mRNA)%>%
  drop_na()

merged2<- left_join(DNAm, Prot)%>%
  drop_na()

merged3<- left_join(mRNA, Prot)%>%
  drop_na()

merged4<- left_join(merged2,merged3)%>%
  drop_na()

merged5<- left_join(merged4,merged)%>%
  drop_na()
merged5=merged5[-1,]
merged5=unique(merged5$Gene)

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/")
library(ggrepel)
tiff('Effect size DNAm vs mRNA for common genes.tiff',
     width =100,
     height = 80,
     units = 'mm',
     res=300)
ggplot(data = merged,
       aes(x = `Effect.size`,
           y = `log2FC.per.unit.of.VO2max`,
           label = Gene)) +
  xlab("% DNAm change per unit of VO2max")+
  ylab("mRNA log2FC per unit of VO2max")+
  geom_point(colour="#1AEDD8")+
  geom_label_repel(data = merged,
                   size = 2,
                   box.padding = unit(0.45, "lines"),
                   point.padding = unit(0.45, "lines"),
                   max.overlaps = 30)+
  geom_hline(yintercept = 0,
             lty="dashed")+
  geom_vline(xintercept = 0,
             lty="dashed")+
  theme_classic()
dev.off()

tiff('Effect size mRNA & Prot for common genes.tiff',
     width =100,
     height = 80,
     units = 'mm',
     res=300)
ggplot(data = merged3,
       aes(x = `EFFECTSIZE` ,
           y = `log2FC.per.unit.of.VO2max`,
           label = Gene)) +
  xlab("Protein logFC change per unit of VO2max")+
  ylab("mRNA log2FC per unit of VO2max")+
  geom_point(colour="#1BD2F7")+
  geom_hline(yintercept = 0,
             lty="dashed")+
  geom_vline(xintercept = 0,
             lty="dashed")+
  geom_label_repel(data = merged3,
                   size = 2,
                   box.padding = unit(0.45, "lines"),
                   point.padding = unit(0.45, "lines"),
                   max.overlaps = 30)+
  theme_classic()
dev.off()

tiff('Effect size DNAm & Prot & mRNA for common genes.tiff',
     width =100,
     height = 80,
     units = 'mm',
     res=300)
ggplot(data = merged5,
       aes(x = `EFFECTSIZE`,
           y = `Effect.size`,
           label = Gene)) +
  xlab("Protein logFC change per unit of VO2max")+
  ylab("% DNAm change per unit of VO2max")+
  geom_point(colour="#CC99FF")+
  geom_point(aes(x = `EFFECTSIZE` ,
                 y = `log2FC.per.unit.of.VO2max`), colour="#1BD2F7")+
  geom_label_repel(data = merged5,
                   size = 2,
                   box.padding = unit(0.45, "lines"),
                   point.padding = unit(0.45, "lines"),
                   max.overlaps = 10)+
  geom_hline(yintercept = 0,
             lty="dashed")+
  geom_vline(xintercept = 0,
             lty="dashed")+
  theme_classic()
dev.off()


tiff('Effect size DNAm & Prot for common genes.tiff',
     width =150,
     height = 100,
     units = 'mm',
     res=300)
ggplot(data = merged2,
       aes(y = `EFFECTSIZE` ,
           x = `Effect.size`,
           label = Gene)) +
  ylab("Protein logFC change per unit of VO2max")+
  xlab("% DNAm change per unit of VO2max")+
  geom_point(colour="#CC99FF")+
  geom_hline(yintercept = 0,
             lty="dashed")+
  geom_vline(xintercept = 0,
             lty="dashed")+
  geom_label_repel(data = merged2,
                   size = 2,
                   box.padding = unit(0.45, "lines"),
                   point.padding = unit(0.45, "lines"),
                   max.overlaps = 30)+
  theme_classic()
dev.off()


#Venn diagram
library("ggVennDiagram")
library("ggplot2")
library("sf")
library("ggvenn")

x <- list(DMGs = DNAm$Gene,
          DEGs = mRNA$Gene,
          DEPs = Prot$Gene)
value=c(DNAm$Gene,mRNA$Gene,Prot$Gene)
value=unique(value)

data_venn <- data.frame(value = value,    # Create example data frame
                        DMGs = FALSE,
                        DEGs = FALSE,
                        DEPs = FALSE)
data_venn$DMGs <- data_venn$value %in% x$DMGs
data_venn$DEGs <- data_venn$value %in% x$DEGs
data_venn$DEPs <- data_venn$value %in% x$DEPs

ggplot(data_venn, aes(A=DMGs, B=DEGs, C=DEPs)) +
  geom_venn(fill_color = c("#1AEDD8", "#1BD2F7", "#CC99FF"), 
            stroke_size = 1, set_name_size = 5, fill_alpha = 0.3, text_size = 5, 
            show_percentage = FALSE, text_color = "#333333")+
  theme_classic()+
  theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank())

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax")
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

merged
"NCEH1"   "ALDH6A1" "HSPA2"   "PARK7"   "CAB39"  

####Use Mitch to contrast results####
#change table for mitch
Meta_Trans_mitch=trans%>%
  mutate(t=log2FC.per.unit.of.VO2max/SE)%>%
  select("log2FC.per.unit.of.VO2max", "SE", t, P.value, FDR, Gene.name)%>%
  drop_na()
rownames(Meta_Trans_mitch)=Meta_Trans_mitch$Gene.name

library("mitch")
library(tidyverse)

prot1 <- prot
colnames(prot1)[1] <- "gene"
head(prot1)
str(prot1)

prot1=prot1%>%
  mutate(t=Effect/StdErr)

trans <- Meta_Trans_mitch
colnames(trans)[6] <- "gene"
head(trans)
str(trans)

#add genes to methylation results
meth1 <- meth%>%
  mutate(t=Effect.size/SE)%>%
  select("CpG","Effect.size","SE","P.value","FDR", "t","Annotated.gene.s.")
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

#Now import with mitch
prot1 <- prot1[,c("gene","t")]
dim(prot1)
trans <- trans[,c("gene","t")]
dim(trans)
meth1 <- meth_fix2[,c("gene","t")]
meth1 <- aggregate(. ~ gene, meth1, sum)
dim(meth1)
l <- list("meth"=meth1, "rna"=trans, "prot"=prot1)
m <- mitch_import(x = l,geneIDcol = "gene",DEtype = "limma")

#Fetch gene sets 
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip", destfile="ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
genesets <- gmt_import("ReactomePathways.gmt")

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/MSigDB files/MSigDB/msigdb_v2022.1.Hs_GMTs")
genesets2<- gmt_import("c5.all.v2023.2.Hs.symbols.gmt")

#Calculate enrichment priority as significance
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax")

res <- mitch_calc(x=m,genesets = genesets ,priority = "significance",minsetsize = 5, resrows = 25)
head(res$enrichment_result,50)
mitch_report(res = res, outfile = "myreport_significance_reactome.html")
mitch_plots(res,outfile="mycharts_significance_reactome.pdf")

#save results as table
res=res$enrichment_result
res_sig=res[res$p.adjustMANOVA<0.05 ,]

#Reproduce plot using results table
res1 = mutate(res, sig=ifelse(res$p.adjustMANOVA<0.05, "FDR<0.05", "Not Sig"))

library("ggrepel")

tiff('Effect size versus significance pathways VO2max reactome.tiff', width =8, height = 6, units = 'in', res=600)
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
colnames(hmapx)=c("DNAm", "mRNA", "Prot")


tiff('Heatplot integration.tiff',
     width =200,
     height = 250,
     units = 'mm',
     res=600)
heatmap.2(hmapx,scale="none",Rowv = T, Colv = FALSE ,margin=c(10, 25),
          cexRow=0.8,trace="none",cexCol=1,col=my_palette,
          dendrogram="none")
dev.off()