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

DEGs <- left_join(DEGs,
                  extrameta)

DEGs <- DEGs %>%
  filter(`Moderators?`!="None")

# Library
library(ggplot2)
library(dplyr)
library(hrbrthemes)

#Extract coeff for age
my_fn <- function(l)
{
  tib <- tibble(ES_exercise = l[[1]]$coeffs["intrcpt","beta"],
                ES_agemod = l[[1]]$coeffs["avg_age","beta"])
  return(tib)
}
names <- DEGs %>%
  filter(str_detect(`Moderators?`,
                    "Age"))%>%
  pull(Gene)
subchronic_muscle <- chronic_muscle[names]
results <- lapply(subchronic_muscle,
                  my_fn)
results <- bind_rows(results)
results <- results %>%
  mutate(Gene = names,
         `20 year-old` = abs(ES_exercise-25*ES_agemod),
         `70 year-old` = abs(ES_exercise+25*ES_agemod))
results <- left_join(results,
                     DEGs %>% select(-`log2FC after training`))%>%
  pivot_longer(cols = contains("year-old"),
               names_to = "Age",
               values_to = "log2FC after training") %>%
  mutate(direction = ifelse(sign(`log2FC per year of age`)==1,
                            "Up DEG",
                            "Down DEG"),
         Gene = factor(Gene))%>%
  mutate(Gene = fct_reorder(Gene,
                            ES_agemod))

# Plot
tiff('Difference in exercise response between young and old at age-related DEGs.tiff',
     width =80,
     height = 120,
     units = 'mm',
     res=300)
ggplot(results,
       aes(`log2FC after training`,
           Gene)) +
  geom_line(aes(group = Gene)) +
  geom_point(aes(color = `Age`),
             size = 2)+
  geom_vline(xintercept = 0,
             lty = 2)+
  scale_colour_manual(values=c("#E69F00", "#009E73"))+
  facet_grid(direction ~ .,
             scales="free_y",
             space="free_y")+
  theme(legend.position="top")
dev.off()


#Extract coeff for sex
my_fn <- function(l)
{
  tib <- tibble(ES_exercise = l[[1]]$coeffs["intrcpt","beta"],
                ES_sexmod = l[[1]]$coeffs["prop_males","beta"])
  return(tib)
}
names <- extrameta %>%
  filter(str_detect(`Moderators?`,
                    "Sex"))%>%
  pull(Gene)
subchronic_muscle <- chronic_muscle[names]
results <- lapply(subchronic_muscle,
                  my_fn)
results <- bind_rows(results)
results <- results %>%
  mutate(Gene = names,
         `Male` = ES_exercise-ES_sexmod/2,
         `Female` = ES_exercise+ES_sexmod/2)
results <- left_join(results,
                     DEGs %>% select(-`log2FC after training`))%>%
  pivot_longer(cols = contains("ale"),
               names_to = "Sex",
               values_to = "log2FC after training") %>%
  mutate(direction = ifelse(sign(`log2FC per year of age`)==1,
                            "Up DEG",
                            "Down DEG"))%>%
  arrange(direction,
          ES_exercise)

# Plot
tiff('Difference in exercise response between males and females at age-related DEGs.tiff',
     width =100,
     height = 130,
     units = 'mm',
     res=300)
ggplot(results,
       aes(`log2FC after training`,
           Gene)) +
  geom_line(aes(group = Gene)) +
  geom_point(aes(color = `Sex`),
             size = 2)+
  geom_vline(xintercept = 0,
             lty = 2)+
  scale_colour_manual(values=c("#CC79A7", "#56B4E9"))+
  facet_grid(direction ~ .,
             scales="free_y",
             space="free_y")
dev.off()

#Extract coeff for modality
my_fn <- function(l)
{
  tib <- tibble(ES_exercise = l[[1]]$coeffs["intrcpt","beta"],
                ES_type = l[[1]]$coeffs["trainingresistance","beta"])
  return(tib)
}
names <- DEGs %>%
  filter(str_detect(`Moderators?`,
                    "Exercise modality"))%>%
  pull(Gene)
subchronic_muscle <- chronic_muscle[names]
results <- lapply(subchronic_muscle,
                  my_fn)
results <- bind_rows(results)
results <- results %>%
  mutate(Gene = names,
         `Endurance` = ES_exercise-ES_type/2,
         `Resistance` = ES_exercise+ES_type/2)
results <- left_join(results,
                     DEGs %>% select(-`log2FC after training`))%>%
  pivot_longer(cols = contains("ance"),
               names_to = "Modality",
               values_to = "log2FC after training") %>%
  mutate(direction = ifelse(sign(`log2FC per year of age`)==1,
                            "Up DEG",
                            "Down DEG"))%>%
  arrange(direction,
          ES_type)

# Plot
tiff('Difference in exercise response between END and RES at age-related DEGs.tiff',
     width =100,
     height = 130,
     units = 'mm',
     res=300)
ggplot(results,
       aes(`log2FC after training`,
           Gene)) +
  geom_line(aes(group = Gene)) +
  geom_point(aes(color = `Modality`),
             size = 2)+
  geom_vline(xintercept = 0,
             lty = 2)+
  scale_colour_manual(values=c("#F0E442", "#D55E00"))+
  facet_grid(direction ~ .,
             scales="free_y",
             space="free_y")
dev.off()


c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
