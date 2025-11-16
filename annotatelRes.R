setwd("/Users/jeremydawe/Documents/School/GoutLab/Prion_Work/")
library(rtracklayer)
library(phastCons100way.UCSC.hg38)
library(ggplot2)

#store conservation scores
phast <- phastCons100way.UCSC.hg38

lRes <- readRDS ("lRes_MIN_COV-200_HS-3.rds")

# Extract calls
calls <- lapply(lRes, function(x) x$calls)

allCalls <- unlist(GRangesList(calls))
grCalls <- reduce(allCalls)
grCalls <- unlist(tile(grCalls, width = 1))

# Get con scores
grCalls <- keepSeqlevels(grCalls, seqlevels(grCalls)[!(seqlevels(grCalls) %in% c("5S_rRNA", "rRNA"))], pruning.mode="coarse")
conGR <- gscores(phast, grCalls)
seqlevelsStyle(conGR) <- "NCBI"
hist(conGR$default, main = "Histogram of Conservation Scores", xlab = "Conservation Score", ylab = "Frequency", col = "purple", border = "black")

# Create four new GRanges for bins
conGR_0 <- subset(conGR, default == 0)
conGR_0_05 <- subset(conGR, default > 0 & default <= 0.5)
conGR_05_1 <- subset(conGR, default > 0.5 & default < 1)
conGR_1 <- subset(conGR, default == 1)

# Function to annotate lRes with bin information
binAssigner <- function(gr, conGR_0, conGR_0_05, conGR_05_1, conGR_1) {
  bin <- integer(length(gr)) # intialize bin w/ NA
  
  bin[queryHits(findOverlaps(gr, conGR_0))] <- 1
  bin[queryHits(findOverlaps(gr, conGR_0_05))] <- 2
  bin[queryHits(findOverlaps(gr, conGR_05_1))] <- 3
  bin[queryHits(findOverlaps(gr, conGR_1))] <- 4
  
  mcols(gr)$Bin <- bin
  return(gr)
}

# Apply function to each sample in lRes
for (i in seq_along(lRes)){
  lRes[[i]]$calls <- binAssigner(lRes[[i]]$calls, conGR_0, conGR_0_05, conGR_05_1, conGR_1)
}

# Consolidate all Calls GR into one object again
calls2 <- lapply(lRes, function(x) x$calls)
allCalls2 <- unlist(GRangesList(calls2))

length(allCalls2[which(allCalls2$Bin == 1),]) #Check size of all bins
length(allCalls2[which(allCalls2$Bin == 2),])
length(allCalls2[which(allCalls2$Bin == 3),])
length(allCalls2[which(allCalls2$Bin == 4),])


# GET TRANSCRIPTION ERROR RATES FOR EACH BIN  

ts <- read.csv("finalSpec.csv")

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
