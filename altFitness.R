setwd("/Users/jeremydawe/Documents/School/GoutLab/Prion_Work/")
siftDir <- "/Users/jeremydawe/Documents/School/GoutLab/Prion_Work/SIFTscores"

library(rtracklayer)
library(data.table)
library(dplyr)
library(phastCons100way.UCSC.hg38)
library(stringr)
library(GenomicScores)
library(GenomicRanges)

# Import lRes and add conservation score and sift score columns
lRes <- readRDS ("lRes_MIN_COV-200_HS-3.rds")

lRes <- lapply(lRes, function(sublist){
  calls <- sublist[[3]]
  candidates <- sublist[[2]]
  mcols(candidates)$ConservationScore <- NA
  mcols(candidates)$SIFTscore <- NA
  mcols(calls)$ConservationScore <- NA
  mcols(calls)$SIFTscore <- NA
  sublist[[3]] <- calls
  sublist[[2]]
  return(sublist)
})

# Annotate lRes with SIFT scores
lRes <- lapply(lRes, function(sublist) {
  calls <- sublist[[3]]
  candidates <- sublist[[2]]       
  if (!inherits(calls, "GRanges")) {
    warning("Third element is not a GRanges object")
    print(class(calls))
    return(sublist)
  }
  chrs <- as.character(unique(seqnames(calls)))
  siftScores <- rep(NA_real_, length(calls))
  candidateScores <- rep(NA_real_, length(candidates))

  for (chr in chrs){
    cat("Working on chromosome: ", chr, "\n")
    file <- file.path(siftDir, paste0(chr, ".filtered.gz"))
    if (!file.exists(file)) {
      warning("SIFTwarnings file not found for chromosome: ", chr)
      next
    }
    
    # Read in SIFT score
    sif <- fread(file)
    setnames(sif, old = "#Position", new = "Position", skip_absent = TRUE)
    sif[, Position := as.integer(Position)]
    sif <- sif[!is.na(SIFT_score)]
    
    # Subset calls for this chromosome
    chrIdx <- which(as.character(seqnames(calls)) == chr)
    chrCalls <- calls[chrIdx]
    
    candidatesIdx <- which(as.character(seqnames(candidates)) == chr)
    chrCandidates <- candidates[candidatesIdx]
    
    # Match and annotate
    matched <- match(start(chrCalls), sif$Position)
    siftScores[chrIdx] <- sif$SIFT_score[matched]
    
    matches <- match(start(chrCandidates), sif$Position)
    candidateScores[candidatesIdx] <- sif$SIFT_score[matches]
  }
  mcols(calls)$SIFTscore <- siftScores
  mcols(candidates)$SIFTscore <- candidateScores
  
  sublist[[3]] <- calls
  sublist[[2]] <- candidates
  return(sublist)
})

##################################
# Annotate with conservation score

phast <- phastCons100way.UCSC.hg38

# Get conservation scores
lRes <- lapply(lRes, function(sublist){
  calls <- sublist[[3]]
  candidates <- sublist[[2]]
  
  # Filter
  unwanted <- c("5S_rRNA", "rRNA")
  calls <- calls[!as.character(seqnames(calls)) %in% unwanted]
  candidates <- candidates[!as.character(seqnames(candidates)) %in% unwanted]
  
  # add conservation scores
  mcols(calls)$ConservationScore <- score(phast, calls)
  mcols(candidates)$ConservationScore <- score(phast, candidates)
  
  # Update
  sublist[[3]] <- calls
  sublist[[2]] <- candidates
  return(sublist)
})

#####################################
# Plot Error Spectrum for SIFT scores

# Function to annotate bin information
lRes <- lapply(lRes, function(sublist) {
  for (i in c(2,3)){
    gr <- sublist[[i]]
    
    bin <- rep(NA_integer_, length(gr))
    sift <- mcols(gr)$SIFTscore
    bin[!is.na(sift) & sift < 0.05] <- 1
    bin[is.na(sift) & sift ==1] <- 2
    mcols(gr)$Bin <- bin
    sublist[[i]] <- gr
  }
  return(sublist)
})

# Calculate transcription error rates
allSubs <- c()

for (sublist in lRes){
  calls <- sublist[[3]]
  base <- mcols(calls)$base
  counts <- mcols(calls)[, c("A", "T", "C", "G")]
  
  altBase <- apply(counts, 1, function(row){
    if (all(row == 0)) return(NA)
    names(row)[which.max(row)]
  })
  
}

################
# Annotate bin information based on both sift and conservation scores
annotater <- function(gr){
  sift <- mcols(gr)$SIFTscore
  cons <- mcols(gr)$ConservationScore
  
  Bin <- rep(NA_character_, length(gr))
  
  # Define conditions
  goodSift <- which(!is.na(sift) & sift <= 0.05)
  badSift <- which(!is.na(sift) & sift > 0.05)
  
  goodCons <- which(!is.na(cons) & cons ==1)
  badCons <- which(!is.na(cons) & cons != 1)
  
  # Good bin
  gIndices <- union( intersect(goodSift, goodCons), union(setdiff(goodSift, c(goodCons, badCons)), setdiff(goodCons, c(goodSift, badSift))))
  
  # Bad bin
  bIndices <- union(intersect(badSift, badCons), union(setdiff(badSift, c(goodCons, badCons)), setdiff(badCons,c(goodSift, badSift))))
  
  #Assign
  Bin[gIndices] <- "Good"
  Bin[bIndices] <- "Bad"
  
  mcols(gr)$Bin <- Bin
  return(gr)
}

for (i in seq_along(lRes)){
  candidates <- lRes[[i]][[2]]
  calls <- lRes[[i]][[3]]
  
  calls <- annotater(calls)
  candidates <- annotater(candidates)
  
  lRes[[i]][[2]] <- candidates
  lRes[[i]][[3]] <- calls
}

# Generate histogram


# Second step
ts <- read.csv("finalSpec-union.csv")

ts <- ts[which(ts$region != 2) , ]
ts <- ts[which(ts$region == 1 | ts$region == 4) , ]
ts$ConservationScore <- "Low"
ts$ConservationScore[which(ts$Bin == 4)] <- "High"
ts$ConservationScore <- factor(ts$ConservationScore, levels = c("Low", "High"))

# Plot a histogram using ggplot
pb <- ggplot(data=ts, aes(fill=ConservationScore, x=mut, y=rate) ) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar( aes(min = min, max = max), width=.2, position=position_dodge(.9) ) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Base substitution", y = "Error rate") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.9,0.85))
plot(pb) 