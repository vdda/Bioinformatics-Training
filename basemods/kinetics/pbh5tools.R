library(pbh5)
library(Hmisc)
library(plyr)
library(Biostrings)
options(scipen=20)
pdfOutPath <- "cmph5_demo/"

# The cmph.h5 file output as part of the basemod workflow has alignments for every read
# and kinetic information for every base in those reads.
# You can examine all of this data in R using the accessor functions in R-pbh5. 
# This tutorial will walk you through how to apply some of these.

# One thing you might want to do is examine the metadata in the cmp.h5 file.  You can
# easily generate a dataframe with stats like, polymerase readlength, subread length, accuracy,
# CCS number of passes, etc.  Here's a simple example plotting the metadata from the two SMRT Cells
# of ecoli used in the basemod example.

get.cmph5.Dataframe <- function(path){   
  print(path)
  cmph5data <- PacBioCmpH5(path)
  accuracy <- getAccuracy(cmph5data) 
  readlength <- getTemplateSpan(cmph5data)  #length of the *reference* in the alignment region
  holeNumber <- getHoleNumbers(cmph5data)
  instrument <- getMachineName(cmph5data)
  movie <- getMovieName(cmph5data)
  movieStamp <- as.character(lapply(as.character(movie), function(x){strsplit(x,split='_c')[[1]][1]}))
  df <- data.frame(accuracy,readlength, holeNumber, instrument, movieStamp) 
  return(df)
}
cmph5path <- c('aligned_reads.cmp.h5')
cmph5Subreads <- ldply(cmph5path, get.cmph5.Dataframe)

# You can see all the cmp.h5 accessor functions here:
# https://github.com/PacificBiosciences/R-pbh5/blob/master/R/cmpH5Utils.R

# How many SMRT Cells were used in the base mod analysis? How many reads aligned from each cell?
table(cmph5Subreads$movieStamp)

# What was the thoughput in bases for each cell?
throughput <- ddply(cmph5Subreads, c("movieStamp"), function(df) round_any(sum(df$readlength)/1000000, 0.1))
names(throughput)[2] <- "Mapped_Yield_in_Mb"
throughput

# One of the cells was run with a size-selected 20kb insert library, and the other is
# from a standard 10 bk library prep.  Are there differences in the subread length or accuracy?
#
# Subread lengths
pdf(paste(pdfOutPath, 'ecoli.readlengths.pdf', sep=""), width=11, height=8.5, onefile=T)
p <- qplot(data=cmph5Subreads, x=readlength, geom="freqpoly", theme="bw", color=movieStamp) +
  labs(x="Mapped Subread Length", color="Movie") +
  theme(legend.position = c(0.8,0.8)) +
  coord_cartesian(xlim=c(-100,15000)) + scale_x_continuous(breaks=seq(0, 14000, 2000))
show(p)
#
p <- qplot(data=cmph5Subreads, x=movieStamp, y=readlength, geom="boxplot", theme="bw", color=movieStamp) +
  labs(y="Mapped Subread Length", x="Movie", color="Movie") +
  theme(legend.position = c(0.85,0.85))
show(p)
dev.off()

# Accuracy
pdf(paste(pdfOutPath, 'mapped_accuracy.pdf', sep=""), width=11, height=8.5, onefile=T)
p <- qplot(data=cmph5Subreads, x=accuracy, geom="freqpoly", theme="bw", color=movieStamp) +
  labs(x="Mapped Accuracy", color="Movie") +
  theme(legend.position = c(0.2,0.9)) +
  coord_cartesian(xlim=c(0.74,0.96)) + scale_x_continuous(breaks=seq(0.75, 0.95, 0.05))
show(p)
dev.off()

### The cmp.h5 file also has kinetic information, which we can also delve into
### using the utilities in R-pbh5. 

# Since there is kinetic information for every base incorporation event in every read,
# loading kinetic information for all the data in both cells would be very slow. 
# Let's instead load just the kinetic info from reads which align to a  window +/- 10 bases
# around a GATC motif.  The 'GATC' from 1569-1572 will do nicely.

ecoliPath <- 'ecoli_K12_MG1655.fasta'
ecoli <- read.DNAStringSet(ecoliPath)
gatcMatches <- vmatchPattern('GATC', ecoli) 
gatcStart <- start(gatcMatches)[[1]][6] - 10  #look at just one motif in the genome
gatcEnd <- end(gatcMatches)[[1]][6] + 10

cmph5 <- PacBioCmpH5(cmph5path)
gatcReads <- getByTemplatePosition(cmph5,idx=getReadsInRange(cmph5, 1, gatcStart, gatcEnd), f = getIPD)
names(gatcReads)
subset(gatcReads, position >1568 & position < 1573 & idx > 160)

# Each idx is a separate read.  There are 119 reads which align to this region. 
levels(factor(gatcReads$idx))

# Our dataframe contain information for the entire length of all 119 reads, not just the bases in the
# window we're interested in.  Let's fix that.

gatc <- subset(gatcReads, position >= gatcStart & position <= gatcEnd)

table(gatc[for_gatc$strand==0,]$ref, gatc[gatc$strand==0,]$position)
table(gatc[for_gatc$strand==1,]$ref, gatc[gatc$strand==1,]$position)

# Plot ipds for the aligned bases
# Remember that the altered ipd is on the strand *opposite* the modification event.

dna <- c(rgb(116, 118, 121,maxColorValue=255),# ltgrey
         rgb(49, 94,161,maxColorValue=255),# dkblue 
         rgb(156, 19,46,maxColorValue=255),#red 
         rgb(6,88, 20,maxColorValue=255),#green 
         rgb(242,175,50,maxColorValue=255)) # yellow

pdf(paste(pdfOutPath, 'kinetics.pdf', sep=""), width=11, height=8.5, onefile=T)
p <- qplot(data=gatc, x=as.factor(position), y=elt, geom="point", theme="bw", color=ref) +
  labs(y="IPD", x="Template Position", color="Base in Reference") +
  facet_grid(strand~.) +
  scale_colour_manual(values=dna) + 
  theme(legend.position = c(0.85,0.85))
show(p)
#
p <- qplot(data=gatc, x=as.factor(position), y=elt, geom="boxplot", theme="bw") +
  labs(y="IPD", x="Template Position") +
  facet_grid(strand~.) +
  theme(legend.position = c(0.85,0.85))
show(p)
dev.off()

# You can see from the plot that the 'T' at position 1571(0 strand), across from the 'A' in 'GATC',
# has a higher value that other T's in this context.  The same is true on the opposite strand.
# Note also that 1565 on the '1' strand, which opposes the 5th base 5' to G*A*TC on the '0' strand, also
# has outlier values.  This is often seen with 6mA.

# In this example, we are looking at IPDRatios of individual reads in our data.  We can also calculate the 
# value that appears in the gff file - the *mean* IPD Ratio of all the reads that include position 1571 
# on the '0' strand.

pos1570_meanIPDRatio <- ddply(subset(gatc, strand==1 & position==1570 & read=='T'),c("position", "strand","read"), summarize, meanIPDRatio=mean(elt))


# We can also plot the distribution of IPD Ratios at each position.  If a sample is uniformly
# methylated at this position, there should be just one mode at the modified position.

motif = subset(gatc, strand ==1 & position >= 1569 & position <= 1572) # look at just the GATC motif
pdf(paste(pdfOutPath, 'one_gatc_IPDRatio_distribution.pdf', sep=""), width=11, height=8.5, onefile=T)
d=qplot(data=motif, x=elt, geom='freqpoly', colour=read) +
  facet_wrap(~ position, ncol=2, scales="free_y") +
  labs(title="The Distribution of IPD Ratios at each postion", x="IPD Ratio", y="Density", colour="Base")
show(d)
dev.off()
