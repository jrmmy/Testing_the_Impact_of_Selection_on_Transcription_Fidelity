# Find frequencies of specific mutations from clinvar
setwd("/Users/jeremydawe/Documents/School/GoutLab/Prion_Work/")

library(AnnotationHub)
library(GenomicRanges)
library(GenomicFeatures)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(ggplot2)

ensTxDb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.112.chr.gff3.gz", format="gff3")
myCDS <- cdsBy(ensTxDb)
cdsGR <- unlist(myCDS)

amyloids <- read.delim("amyloidSN.txt")
alzheimers <- read.delim("clinVar_Alzheimers_SN_P_LP.txt")
cardioVascular <- read.delim("clinVar_Cardiovascular_SN_P_LP.txt")
covid <- read.delim("clinVar_Covid_SN_P_LP.txt")
neurodegenerative <- read.delim("clinVar_Neurodegenerative_SN_P_LP.txt")
prion <- read.delim("clinVar_Prion_SN_P_LP.txt")

countMutations <- function(amyloids){
  spdi_split <- strsplit(amyloids$Canonical.SPDI, ":", fixed = TRUE)
  amyloids$ref <- sapply(spdi_split, `[`, 3)
  amyloids$var <- sapply(spdi_split, `[`, 4)
  
  # Make amyloids GRANGES
  amyloids$strand <- "*"
  grAmyloids <- makeGRangesFromDataFrame(amyloids, keep.extra.columns = TRUE, seqnames.field = "GRCh38Chromosome", start.field = "GRCh38Location", end.field = "GRCh38Location", na.rm = TRUE)
  
  # Select for cds positive and cds negative strands
  positiveCDS = cdsGR[strand(cdsGR) == "+"]
  negativeCDS = cdsGR[strand(cdsGR) == "-"]
  
  # Update grAmyloids with strand information
  positiveOverlaps <- findOverlaps(grAmyloids, positiveCDS)
  positiveInformation<- strand(positiveCDS)[subjectHits(positiveOverlaps)]
  strand(grAmyloids)[queryHits(positiveOverlaps)] <- positiveInformation
  
  negativeOverlaps <- findOverlaps(grAmyloids, negativeCDS)
  negativeInformation <- strand(negativeCDS)[subjectHits(negativeOverlaps)]
  strand(grAmyloids)[queryHits(negativeOverlaps)] <- negativeInformation
  
  validStrands <- strand(grAmyloids) == "+" | strand(grAmyloids) == "-"
  grAmyloids <- grAmyloids[validStrands]
  
  # Find Negative strands in grAmyloids and complement them
  negStrandIndices <- which(strand(grAmyloids)== "-")
  
  complementer <- function(base){
    baseComplement <- chartr("ATCG", "TAGC", base)
    return(baseComplement)
  }
  
  grAmyloids$ref[negStrandIndices] <- sapply(grAmyloids$ref[amyloidNegLocations], complementer)
  grAmyloids$var[negStrandIndices] <- sapply(grAmyloids$var[amyloidNegLocations], complementer)
  
  # Get stats on specific mutations
  bases <- c("A", "T", "C", "G")
  allMutations <- expand.grid(ref = bases, var = bases, stringsAsFactors = FALSE)
  allMutations <- allMutations[allMutations$ref != allMutations$var, ]
  allMutations <- paste(allMutations$ref, allMutations$var, sep = ">")
  
  refs <- as.character(grAmyloids$ref)
  var <- as.character(grAmyloids$var)
  
  mutationType <- paste(refs, var, sep = ">")
  mutationType <- factor(mutationType, levels = allMutations)
  
  mutCounts <- table(mutationType)
  mutProportions <- prop.table(mutCounts)
  
  return(mutProportions)
}

neurodegenerative <- neurodegenerative[neurodegenerative$GRCh38Location != "12340279 - 12340280",]

allDFs <- list(Amyloids = amyloids, Alzheimers = alzheimers, CardioVascular = cardioVascular, neurodegenerative = neurodegenerative, prion = prion)
allCounts <- lapply(allDFs, countMutations)

## Next step is to make this all sexy with all the different clinvar datasets, convert to proportions, calculate error bars and graph them all labeled

for (sample in names(allCounts)){
  allCounts[[sample]]$Sample <- Sample
}

combinedCounts <- do.call(rbind, allCounts)
combinedCounts <- as.data.frame(combinedCounts)

# Convert Data Frame

combinedCounts_tidy <- combinedCounts %>%
  rownames_to_column(var = "Sample") %>%  # Convert row names to a column
  pivot_longer(cols = -Sample, names_to = "MutationType", values_to = "Proportion")

ggplot(combinedCounts_tidy, aes(x = MutationType, y = Proportion, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +  # Side-by-side bars
  labs(
    title = "Mutation Proportions by Sample",
    x = "Mutation Type",
    y = "Proportion",
    fill = "Sample"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels


