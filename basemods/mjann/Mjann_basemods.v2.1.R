# This script is for use with SMRTPortal 2.1 modifications.csv and modifications.gff output files.
# It reads gff3 compliant files.
# Meredith Ashby, Khai Luong, Jonas Korlach 02/2014

source("BaseModFunctions.v2.1.R")

#### USER ENTERED VARIABLES HERE
# NOTE: The modification detection should be done on the reference you wish to use for circos, or the
# coordinates will not match up.

refPath <- 'mjann/Methanocaldococcus_jannaschii_DSM2661.fasta'
gff <- 'mjann/mjann.modifications.gff'
csv <- 'mjann/mjann.modifications.csv'
outputPath <- 'mjann/mj.'

#### READ IN THE GFF FILE
hits <- readModificationsGff(gff)

#### IDENTIFY AND EXAMINE THE HIGH CONFIDENCE HITS

###PAGE 1
# Plotting the base modification coverage vs scores on separate axie by base let's us consider applying different cutoffs in different channels.

pdf(paste(outputPath, '1.ModQV_vs_Coverage_scatter_plot.pdf', sep=""), width=11, height=8.5, onefile=T)
p <- qplot(coverage, score, colour=CognateBase, size=I(0.5), data=subset(hits, CognateBase %in% c('A','T','G','C') & score > 20)) +
  facet_wrap(~CognateBase) +
  labs(title="Modification QV vs. Coverage", x="Per Strand Coverage", y="Modification QV", colour="Cognate Base") +
  theme(legend.position = c(0.1,0.9))
show(p)
dev.off()


##################### OPTIONAL ADDITIONAL FILTERING OF HITS FROM THE GFF FILE

###PAGE 2
# Iteratively adjust the slope for each base, determining which detections to include in motif finding depending on your level of coveraqe.
# Your selection of of the minimum QV can affect the % detected for each motif.  In SMRT Portal or with Motif Maker, you
# can selected hits using just the minumim QV; in R you can set a cutoff that takes into account the coverage dependence of the modificaiton QV 
# score.

cutoffs <- data.frame(CognateBase=c('A','C','G','T'), slope=c(60/150, 30/150, 110/150, 110/150), intercept=c(0,0,0,0), x=c(50, 50, 50, 50))
cutoffs

# change the values for slope, intercept, and x for each base. X is a vertical line on the plot and
# specifies a minimum coverage, independent of score. 'Slope' and 'intercept' define an additional
# minimum score which depends on coverage. For a straight line cut independent of coverage, use
# slope = 0 and intercept = desired minimum score.

# This plot shows how the cutoff choices made above will filter hits.

pdf(paste(outputPath, '2.DataFiltering.pdf', sep=""), width=11, height=8.5, onefile=T)
equations <- data.frame(CognateBase=c('A','C','G','T'), equation=formatEquationText(cutoffs)) # print formatting, for the plot.
p <- qplot(data=subset(hits, CognateBase %in% c('A','T','G','C')), x=coverage, y=score, colour=CognateBase, size=I(0.5)) + 
  facet_wrap(~CognateBase, scales="free_y") + 
  labs(title="Subsetting the Potential Hits to Limit False Positives", x="Per Strand Coverage", y="Modification QV", colour="Cognate Base") +
  theme(legend.position = c(0.85,0.3)) +
  geom_abline(data=cutoffs, aes(slope=slope, intercept=intercept)) +
  geom_vline(data=cutoffs, aes(xintercept=x)) +
  geom_text(data=equations, aes(label=equation, x=30, y=200), vjust=0, hjust=0)
show(p)
dev.off()

### Select and order a subset of hits using the above cutoffs.
workHits <- aboveThreshholdHits(hits, cutoffs)


# At this point you may wish to repeat motif finding using just your selected hits.  To use Motif Maker, you
# will need a new gff file, so let's write one.

gffWrite(workHits, paste(outputPath,"selected.modifications.gff",sep=''))

###### INPUT MOTIF INFORMATION FROM SMRTPORTAL OR MOTIFFINDER

# You can use the functions below to further characterize motifs of interest - modified or unmodified.
# Input the motifs you wish to examine further below to use R for additional analysis

motifs <- c()
positions <- c()

# Expanded the degenerate bases.

motifs <- c()
positions <- c()
motifLabels <- labelContexts(workHits$context, motifs, positions)
table(motifLabels)  # reports how many motifs are found within the 41 base contexts in the gff file.


#### GENOME ANNOTATION

# Once you've found the motifs of interest among the hits, calculate how many instances of each motif within
# the entire genome were detected as modified. Merge all-genome motifs (identified using the reference fasta)
# with the results of the gff/contexts search to see what fraction of motifs in the genome are modified. 

genomeAnnotations <- genomeAnnotation(refPath, motifs, positions) #find the positions of all motifs
mm <- merge(workHits, genomeAnnotations, all = TRUE)  
mm$motif[is.na(mm$motif)] <- 'no_motif'
mm$type[is.na(mm$type)] <- 'not_detected'
mm$CognateBase[is.na(mm$CognateBase)] <- 'native'
mm$type <- factor(mm$type, c("m6A","m4C","modified_base","not_detected"))
table(mm$type, mm$motif)

###PAGE 3
# This shows the hit filtering plot above, but now faceted by motif and colored by event type.
pdf(paste(outputPath, '3.ModQV_vs_Coverage_by_Motif.pdf', sep=""), width=11, height=8.5, onefile=T)
p <- qplot(data=mm, x=coverage, y=score, colour=type, size=I(0.5)) +
  labs(title="ModQV Score vs. Coverage by Motif", x="Per Strand Coverage", y="Modification QV", colour="Type") +
  facet_wrap(~motif, scales="free_y") + guides(size=FALSE)
show(p)
# 
p <- qplot(data=subset(mm, CognateBase!="native"), x=coverage, y=score, colour=CognateBase, size=I(0.5)) +
  labs(title="ModQV Score vs. Coverage by Motif", x="Per Strand Coverage", y="Modification QV", colour="Cognate Base") +
  facet_wrap(~motif, scales="free_y") + guides(size=FALSE)
show(p)
dev.off()

### LOOK AT WHAT *ISN'T* IN THE GFF FILE
# Plot the score vs coverage information for all motif positions found in the reference, including those that aren't hits.
# To do this, we have to get kinetic information on unmethylated motifs from the modifications.csv file.

allKinData <- csv.get(csv)
mga <- getUnmodifiedMotifKin(mm, allKinData)


###PAGE 4
# Compare the score distributions for motif cognate-base positions to scores in the 'no_motif' set
pdf(paste(outputPath, '4.ModificationQV_distribution_by_Motif.pdf', sep=""), width=11, height=8.5, onefile=T)
d=qplot(data=subset(mga, type!="not_detected"), x=score, colour=motif, geom='density') +
  theme(legend.position = c(0.85,0.75)) +
  labs(title="Score Distributions of Each Detected Motif", x="Motif Score", y="Density", colour="Motif")
show(d)
#
# now with facets (useful when there's lots of motifs)
d=qplot(data=mga, x=score, geom='freqpoly', colour=type) +
  facet_wrap(~ motif, ncol=2, scales="free_y") +
  labs(title="Score Distributions for Motifs by Type", x="Motif Score", y="Density", colour="Type")
show(d)
dev.off()

###PAGE 5
# Compare the IPD ratio distributions of the found motifs. Consulting IPDRatios is very helpful when looking at
# at differences in methylation between two differently treated samples, as is isn't coverage dependent.
pdf(paste(outputPath, '5.IPDratio_Distribution_by_Motif.pdf', sep=""), width=11, height=8.5, onefile=T)
d=qplot(data=subset(mga, type!='Not_Detected'), IPDRatio, colour=motif, geom='density') +
  theme(legend.position = c(0.85,0.75)) +
  labs(title="IPD Ratio Distributions by Motif", x="IPD Ratio", y="Density", colour="Motif")
show(d)
#
# now with facets (useful when there's lots of motifs)
d=qplot(data=mga, x=IPDRatio, colour=type, geom='freqpoly') +
  facet_wrap(~ motif, ncol=2, scales="free_y") +
  labs(x="IPD Ratio", y="Density", colour="Type")
show(d)
dev.off()

#  What can these plots tells us about the confidence we have in the motifs we've uncovered? Which might we strike out?

motifs <- c()
positions <- c()
motifLabels <- labelContexts(workHits$context, motifs, positions)
table(motifLabels)  # reports how many motifs are found within the 41 base contexts in the gff file.


## Since we've changed our filtering criteria, let's create a new summary of results, similar to the table SMRTPortal produces
nGenome <- ddply(mga, c('motif'), summarise, nGenome=length(motif))
nDetected <- ddply(subset(mga, type!='not_detected'), c('motif'), summarise, nDetected=length(motif))
meanScore <- ddply(subset(mga, type!='not_detected'), c('motif'), summarise, meanScore=mean(score))
meanIPDRatio <- ddply(subset(mga, type!='not_detected'), c('motif'), summarise, meanIPDRatio=mean(IPDRatio))
meanCoverage <- ddply(subset(mga, type!='not_detected'), c('motif'), summarise, meanCoverage=mean(coverage))

a <- merge(nGenome,nDetected)
b <- merge(meanScore, meanIPDRatio)
c <- merge(a,b)
summary <- merge(c, meanCoverage)
summary$fraction <- summary$nDetected/summary$nGenome
write.table(summary, file = paste(outputPath, 'MotifSummary.csv', sep=''), quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)


# R is a very useful tools for making nearly any type of plot.  The 'names' command reveals
# all the information held in the dataframe 'mga' that you may be interested in further examining.
# Check out the ggplot online manual for documentation and additional plotting options:
# docs.ggplot2.org/current/
names(mga)
