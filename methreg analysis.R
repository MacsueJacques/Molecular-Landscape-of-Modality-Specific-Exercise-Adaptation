#MethReg for TFs directionality
library(MethReg)

#First load all raw data for analysis DNA Methylation and Transcriptome
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/E-MTAB-11282")
B_E_MTAB <- read.table("E-MTAB-11282_normalized_batch_corrected.txt")

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Gene SMART")
#Load data
M_Gene_SMART <- read.table("GeneSMART_M.txt")
library(minfi)
B_Gene_SMART <- ilogit2(M_Gene_SMART)

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/EXACT")
B_EXACT <- read.table("EXACT_normalized_batch_corrected.txt")

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/GSE213029")
B_GSE213029 <- read.table("GSE213029_normalized_batch_corrected_beta.txt")

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/GSE60655")
B_GSE60655 <- read.table("GSE60655_B.txt")

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/GSE213029")
B_GSE213029 <- read.table("GSE213029_normalized_batch_corrected_beta.txt")

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/epiH")
B_EpiH <- read.table("EpiH Normalized beta values after batch correction.txt")

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/trainSMART")
B_trainSMART <- read.table("beta trainsmart.txt")

setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Wellderly")
B_Wellderly <- read.table("beta wellderly.txt")


#Load annotation
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Annotation")
annotation <- read.delim("Annotation.txt")
annotation$EPICv2_ID=paste(annotation$probeID, annotation$EPICv2, sep="_")

#Load results to use only significant CpGs
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Methylation/Metal VO2max (withouth GeneSMARTv2)")
meta_res_robust=read.delim(file="METAL_muscle_VO2max_withoutGSv2.txt")
meta_res_robust=meta_res_robust%>%
  dplyr::filter(Significance=="Sig")

annotation_filtered=annotation%>%
  dplyr::filter(Unique_ID %in% meta_res_robust$Unique.ID)

annotation3=annotation_filtered[,c("probeID", "EPICv2_ID")]

#Filter B tables to include only significant CpGs
B_E_MTAB=B_E_MTAB[annotation_filtered$probeID,]
B_E_MTAB=drop_na(B_E_MTAB)
B_E_MTAB=B_E_MTAB[!duplicated(B_E_MTAB),]
B_E_MTAB$probeID=rownames(B_E_MTAB)

B_EXACT=B_EXACT[annotation_filtered$probeID,]
B_EXACT=drop_na(B_EXACT)
B_EXACT=B_EXACT[!duplicated(B_EXACT),]
B_EXACT$probeID=rownames(B_EXACT)

B_Gene_SMART=B_Gene_SMART[annotation_filtered$probeID,]
B_Gene_SMART=drop_na(B_Gene_SMART)
B_Gene_SMART=B_Gene_SMART[!duplicated(B_Gene_SMART),]
B_Gene_SMART$probeID=rownames(B_Gene_SMART)

B_GSE213029=B_GSE213029[annotation_filtered$probeID,]
B_GSE213029=drop_na(B_GSE213029)
B_GSE213029=B_GSE213029[!duplicated(B_GSE213029),]
B_GSE213029$probeID=rownames(B_GSE213029)

B_GSE60655=B_GSE60655[annotation_filtered$probeID,]
B_GSE60655=drop_na(B_GSE60655)
B_GSE60655=B_GSE60655[!duplicated(B_GSE60655),]
B_GSE60655$probeID=rownames(B_GSE60655)

B_trainSMART=B_trainSMART[annotation_filtered$EPICv2_ID,]
B_trainSMART=drop_na(B_trainSMART)
B_trainSMART=B_trainSMART[!duplicated(B_trainSMART),]
B_trainSMART$EPICv2_ID=rownames(B_trainSMART)
B_trainSMART=full_join(B_trainSMART, annotation3, by="EPICv2_ID")
B_trainSMART=drop_na(B_trainSMART)

B_Wellderly=B_Wellderly[annotation_filtered$EPICv2_ID,]
B_Wellderly=drop_na(B_Wellderly)
B_Wellderly=B_Wellderly[!duplicated(B_Wellderly),]
B_Wellderly$EPICv2_ID=rownames(B_Wellderly)
B_Wellderly=full_join(B_Wellderly, annotation3, by="EPICv2_ID")
B_Wellderly=drop_na(B_Wellderly)

B_EpiH=B_EpiH[annotation3$EPICv2_ID,]
B_EpiH=drop_na(B_EpiH)
B_EpiH=B_EpiH[!duplicated(B_EpiH),]
B_EpiH$EPICv2_ID=rownames(B_EpiH)
B_EpiH=full_join(B_EpiH, annotation3, by="EPICv2_ID")
B_EpiH=drop_na(B_EpiH)


combined_B=full_join(B_E_MTAB, B_EXACT, by="probeID")
combined_B=full_join(combined_B, B_Gene_SMART, by="probeID")
combined_B=full_join(combined_B, B_GSE213029,by="probeID")
combined_B=full_join(combined_B, B_GSE60655,by="probeID")
combined_B=full_join(combined_B, B_EpiH, by="probeID")
combined_B=full_join(combined_B, B_trainSMART, by="probeID")
combined_B=full_join(combined_B, B_Wellderly, by="probeID")
combined_B=combined_B[!duplicated(combined_B$probeID),]
rownames(combined_B)=combined_B$probeID
combined_B=combined_B[,-c(27,456,525,550)]

dna.meth_se<-make_dnam_se(dnam = combined_B, genome = "hg38", arrayType="EPIC",
                           betaToM = FALSE, # transform beta to m-values 
                           verbose = FALSE)


#Gene expression data
setwd("//ad.monash.edu/home/User050/mjac0029/Desktop/Papers in progress/Methylation Multi OMICS skeletal muscle VO2 and exercise/Datasets VOmax/Transcriptome")

expr_avg_ID=rownames(expr_avg)
expr_avg=as.data.frame(expr_avg)
expr_avg$ID=expr_avg_ID

expr_avg_ID2=rownames(expr_avg2)
expr_avg2=as.data.frame(expr_avg2)
expr_avg2$ID=expr_avg_ID2

expr_avg_ID3=rownames(expr_avg3)
expr_avg3=as.data.frame(expr_avg3)
expr_avg3$ID=expr_avg_ID3

combined_exp=full_join(expr_avg, expr_avg2, by="ID")
combined_exp=full_join(combined_exp, expr_avg3, by="ID")

combined_exp=drop_na(combined_exp)
rownames(combined_exp)=combined_exp$ID

geneSymbols <-  combined_exp$ID

library(EnsDb.Hsapiens.v79)

geneIDs2 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= geneSymbols, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
colnames(geneIDs2)=c("ID", "ENSEMBL")
combined_exp=full_join(combined_exp, geneIDs2, by="ID")

rownames(combined_exp)= combined_exp$ENSEMBL
combined_exp=combined_exp[,-c(119, 274)]
combined_exp=as.matrix(combined_exp)
combined_exp=combined_exp[!grepl("LRG", rownames(combined_exp)),]

gene.exp_se=make_exp_se(exp=combined_exp, genome = "hg38",verbose = FALSE)


#create triplets
triplet.promoter <- create_triplet_distance_based(
  region = dna.meth_se,
  target.method = "genes.promoter.overlap",
  genome = "hg38",
  target.promoter.upstream.dist.tss = 2000,
  target.promoter.downstream.dist.tss = 2000,
  motif.search.window.size = 400,
  motif.search.p.cutoff  = 1e-08,
  cores = 1  
)

triplet.distal.window <- create_triplet_distance_based(
  region = dna.meth_se,
  genome = "hg38", 
  target.method = "window",
  target.window.size = 500 * 10^3,
  target.rm.promoter.regions.from.distal.linking = TRUE,
  motif.search.window.size = 500,
  motif.search.p.cutoff  = 1e-08,
  cores = 1
)

triplet.distal.nearby.genes <- create_triplet_distance_based(
  region = dna.meth_se,
  genome = "hg38", 
  target.method = "nearby.genes",
  target.num.flanking.genes = 5,
  target.window.size = 500 * 10^3,
  target.rm.promoter.regions.from.distal.linking = TRUE,
  motif.search.window.size = 400,
  motif.search.p.cutoff  = 1e-08,
  cores = 1  
)

#here need the same number of samples (re-run with Gene SMART)
results.interaction.model <- interaction_model(
  triplet = triplet.promoter, 
  dnam = dna.meth_se,
  exp = gene.exp_se,
  dnam.group.threshold = 0.1,
  sig.threshold = 0.05,
  fdr = T,
  stage.wise.analysis = FALSE,
  filter.correlated.tf.exp.dnam = F,
  filter.triplet.by.sig.term = T
)

#Repeat analysis with only gene SMART to see relationship with VO2max
load("E:/Macsue - Hard drive/PhD Victoria University/Papers to work on/Methylation Meta-Analyses/Transcriptomics/Transx_data.rda")
library(tidyverse)
names(y.Norm)
GE <- y.Norm$counts
ncol(GE)

GE_samples <- str_replace(colnames(GE),"P0","PO")
GE_samples <- str_replace(GE_samples,"P4","4WP")
Meth_samples <- colnames(B_Gene_SMART)
Meth_samples<- gsub('\\.', '_', Meth_samples)
colnames(B_Gene_SMART)=Meth_samples

samples_common <- intersect(GE_samples,Meth_samples)

#Filter tables to include onlu samples in common
colnames(GE)=GE_samples

GE_rownames=rownames(GE)
GE=GE[,samples_common]

B_Gene_SMART_com=B_Gene_SMART[,samples_common]

#Make summarized experiment
gene.exp_se=make_exp_se(exp=GE, genome = "hg38",verbose = FALSE)
dna.meth_se<-make_dnam_se(dnam = B_Gene_SMART_com, genome = "hg38", arrayType="EPIC",
                          betaToM = FALSE, # transform beta to m-values 
                          verbose = FALSE)

#create triplets
triplet.promoter_GS <- create_triplet_distance_based(
  region = dna.meth_se,
  target.method = "genes.promoter.overlap",
  genome = "hg38",
  target.promoter.upstream.dist.tss = 2000,
  target.promoter.downstream.dist.tss = 2000,
  motif.search.window.size = 400,
  motif.search.p.cutoff  = 1e-08,
  cores = 1  
)

triplet.distal.window_GS <- create_triplet_distance_based(
  region = dna.meth_se,
  genome = "hg38", 
  target.method = "window",
  target.window.size = 500 * 10^3,
  target.rm.promoter.regions.from.distal.linking = TRUE,
  motif.search.window.size = 500,
  motif.search.p.cutoff  = 1e-08,
  cores = 1
)

triplet.distal.nearby.genes_GS <- create_triplet_distance_based(
  region = dna.meth_se,
  genome = "hg38", 
  target.method = "nearby.genes",
  target.num.flanking.genes = 5,
  target.window.size = 500 * 10^3,
  target.rm.promoter.regions.from.distal.linking = TRUE,
  motif.search.window.size = 400,
  motif.search.p.cutoff  = 1e-08,
  cores = 1  
)

#here need the same number of samples (re-run with Gene SMART)
results.interaction.model <- interaction_model(
  triplet = triplet.promoter_GS, 
  dnam = dna.meth_se,
  exp = gene.exp_se,
  dnam.group.threshold = 0.1,
  sig.threshold = 0.05,
  fdr = T,
  stage.wise.analysis = FALSE,
  filter.correlated.tf.exp.dnam = F,
  filter.triplet.by.sig.term = T
)

#Visualise
plots <- plot_interaction_model(
  triplet.results = results.interaction.model[150,], 
  dnam = dna.meth_se, 
  exp = gene.exp_se,
  dnam.group.threshold = 0.25
)
