# This script is for use with SMRTPortal 2.0 modifications.csv and modifications.gff output files.
# It reads gff3 compliant files.
# Meredith Ashby, Khai Luong, Jonas Korlach 06/2013

source("/data/BaseModFunctions.R")

#### USER ENTERED VARIABLES HERE
# NOTE: The modification detection should be done on the reference you wish to use for circos, or the
# coordinates will not match up.

refPath <- 'ecoli_K12_MG1655.fasta'
gff <- 'modifications.gff'
csv <- 'modifications.csv' 
pdf_out_path <- 'analysis/ecoli.'
circosPrefix <- 'circos/ecoli.'
confModel <- 'oneMotif.conf'

#### READ IN THE GFF FILE
hits <- readModificationsGff(gff)
workHits <- hits # That took forever! Let's keep a record of the whole gff file by making a copy we can play with.

##  Have a look at what is in our dataframe after loading the gff file.
names(workHits)

#### EXAMINE THE HIGH CONFIDENCE HITS IN THE GFF FILE

###PAGE 1
# Distribution of base modification scores broken out by base.
pdf(paste(pdf_out_path, '1.ModQV_distribution.pdf', sep=""), width=11, height=8.5, onefile=T)
d <- qplot(score, colour=CognateBase, geom='freqpoly', data=subset(hits, CognateBase %in% c('A','T','G','C')), binwidth=1) +
  coord_cartesian(xlim=c(0,420)) +
  scale_x_continuous(breaks=seq(100,400,100)) +
  labs(x="Modification QV", y="Counts", title="Modification QV Distribution")
show(d)
dev.off()

###PAGE 2
# Scatter plot of base modification score vs coverage on one set of axes, colored by base
pdf(paste(pdf_out_path, '2.ModQV_vs_Coverage_scatter_plot.pdf', sep=""), width=11, height=8.5, onefile=T)
p <- qplot(coverage, score, colour=CognateBase, size=I(2.0), data=subset(hits, CognateBase %in% c('A','T','G','C') & score > 20))+
  labs(title="Modification QV vs. Coverage", x="Per Strand Coverage", y="Modification QV", colour="Cognate Base") +
  theme(legend.position = c(0.15,0.8))
show(p)
dev.off()


##################### OPTIONAL ADDITIONAL FILTERING OF HITS FROM THE GFF FILE

###PAGE 3
# Scatter plot, allowing USER TO CHOOSE THE CUTOFFS which will feed forward into downstream analysis. Iteratively adjust the slope for each base
# As previously mentioned, your selection of of the minimum QV can affect the % detected for each motif.  In SMRT Portal or with Motif Finder, you
# can selected hits using just the minumim QV; in R you can set a cutoff that takes into account the coverage dependence of the modificaiton QV 
# score.

cutoffs <- data.frame(CognateBase=c('A','C','G','T'), slope=c(110/150, 110/150, 110/150, 110/150), intercept=c(0,0,0,0), x=c(20, 20, 20, 20))
cutoffs

# change the values for slope, intercept, and x for each base. X is a vertical line on the plot and
# specifies a minimum coverage, independent of score. 'Slope' and 'intercept' define an additional
# minimum score which depends on coverage. For a straight line cut independent of coverage, use
# slope = 0 and intercept = desired minimum score.


# This plot shows how the cutoff choices made above will filter hits.

pdf(paste(pdf_out_path, '3.Filtering_by_ModQV.pdf', sep=""), width=11, height=8.5, onefile=T)
equations <- data.frame(CognateBase=c('A','C','G','T'), equation=formatEquationText(cutoffs)) # print formatting, for the plot.
p <- qplot(data=subset(hits, CognateBase %in% c('A','T','G','C')), x=coverage, y=score, colour=CognateBase, size=I(0.5)) + 
  facet_wrap(~CognateBase) + scale_y_continuous(breaks=seq(0,200,50)) +
  labs(title="Subsetting the Potential Hits to Limit False Positives", x="Per Strand Coverage", y="Modification QV", colour="Cognate Base") +
  theme(legend.position = c(0.85,0.3)) +
  geom_abline(data=cutoffs, aes(slope=slope, intercept=intercept)) +
  geom_vline(data=cutoffs, aes(xintercept=x)) +
  geom_text(data=equations, aes(label=equation, x=30, y=200), vjust=0, hjust=0)
show(p)
dev.off()

### Select and order a subset of hits using the above cutoffs.
workHits <- aboveThreshholdHits(hits, cutoffs)
workHits <- workHits[order(workHits$score, decreasing=T),]

###### INPUT MOTIF INFORMATION FROM SMRTPORTAL OR MOTIFFINDER

# SMRTPortal performs automated modification identification and motif finding on
# genomic positions with significantly altered sequencing kinetics.
# You can use the functions below to further characterize motifs of interest - modified or unmodified.
# Input the motifs you wish to examine further below to use R for additional analysis, or to generate
# circos plots.

motifs <- c('GATC','GCACNNNNNNGTT','AACNNNNNNGTGC')
positions <- c(2,3,2) 
motifLabels <- labelContexts(workHits$context, motifs, positions)
table(motifLabels)  # reports how many motifs are found within the 41 base contexts in the gff file.

#  We see that there are a number of hits in the gff file that aren't in the SMRT Portal 6mA motifs.
#  Some of them may be signal occurring 5 positions 5' of 6mA.  SMRT Portal removes these from further analysis
#  when a hit is recognized as a 6mA, but here we have to do it ourselves.

motifs <- c('GATC', 'NNNNGATC', 'GCACNNNNNNGTT', 'NNNGCACNNNNNNGTT', 'AACNNNNNNGTGC', 'NNNNAACNNNNNNGTGC')
positions <- c(2,1,3,1,2,1) 
motifLabels <- labelContexts(workHits$context, motifs, positions)
table(motifLabels)  

#  What about the other motifs SMRT Portal reports?

motifs <- c('GATC', 'GCACNNNNNNGTT', 'AACNNNNNNGTGC','CCTGGYA','CCTGGGRR','CCWGGAAYR','CCAGGYAD','DHCCTGGBB')
positions <- c(2,3,2,1,1,2,1,3)
motifLabels <- labelContexts(workHits$context, motifs, positions)
table(motifLabels) 

#  Let's look more closely at the remaining un-bucketed hits, as well as the rather skethcy looking
#  modified 'C' hits from SMRT Portal

genomeAnnotations <- genomeAnnotation(refPath, motifs, positions) #find the positions of all motifs in the genome.
mm <- merge(workHits, genomeAnnotations, all = TRUE)  
mm$motif[is.na(mm$motif)] <- 'no_motif'
mm$type[is.na(mm$type)] <- 'not_detected'
mm$CognateBase[is.na(mm$CognateBase)] <- 'native'
table(mm$type, mm$motif)

#  We can see that the motifs with modified 'C' are for the most part not recognized as 4mC, and there are
#  many motif instances in the genome which are not detected as modified.

ddply(subset(mm, type!="not_detected"), c('motif','CognateBase'), summarize, meanScore=mean(score))

# The 'C' motifs are also low scoring.

ddply(subset(mm, type!="not_detected"), c('motif','CognateBase'), summarize, meanIPDRatio=mean(IPDRatio))

# ...which reflects their low IPDRatios

ddply(mm, c('motif','CognateBase'), nrow)

# There are a lot of C cognate bases among the 'no_motif' hits.  All of this leads me to believe that the 'C'
# motifs reported in SMRT Portal are hinting at a less restricted, but only weakly detected 5mC motif.  In fact,
# all of the motifs contain 'CCWGG' in their core, and CCWGG is the recognition sequence for Dcm methylase.
# Let's delete SMRT Portal's suggested motifs and look at that instead.  To robustly detect 5mC, we would have to Tet1 treat 
# our sample.

motifs <- c('GATC', 'GCACNNNNNNGTT', 'AACNNNNNNGTGC','CCWGG')
positions <- c(2,3,2,2)
motifLabels <- labelContexts(workHits$context, motifs, positions)
table(motifLabels) 

# As a final indication that this motif is 5mC, the strongest kinetic perturbation for 5mC is not directly over the modified C. Let's look for the
# adjacent footprint signals as well.

motifs <- c('GATC', 'GCACNNNNNNGTT', 'AACNNNNNNGTGC','CCWGG', 'CCWGG','NCCWGG','NNCCWGG')
positions <- c(2,3,2,1,2,1,1)
motifLabels <- labelContexts(workHits$context, motifs, positions)
table(motifLabels) 
genomeAnnotations <- genomeAnnotation(refPath, motifs, positions) #find the positions of all motifs
mm <- merge(workHits, genomeAnnotations, all = TRUE)  
mm$motif[is.na(mm$motif)] <- 'no_motif'
mm$type[is.na(mm$type)] <- 'not_detected'
mm$CognateBase[is.na(mm$CognateBase)] <- 'native'
table(mm$type, mm$motif)
ddply(mm, c('motif','CognateBase'), nrow)

# That's fewer unbinned hits.  

#### GENOME ANNOTATION

# Once you've found the motifs of interest among the hits, calculate how many instances of each motif within
# the entire genome were detected as modified. Merge all-genome motifs (identified using the reference fasta)
# with the results of the gff/contexts search to see what fraction of motifs in the genome are modified. 

###PAGE 4
# This shows the hit filtering plot above, but now faceted by motif and colored by event type.
pdf(paste(pdf_out_path, '4.ModQV_vs_Coverage_by_Motif.pdf', sep=""), width=11, height=8.5, onefile=T)
p <- qplot(data=mm, x=coverage, y=score, colour=type, xlim=c(0,120), ylim=c(0,200), size=I(1.5)) +
  labs(title="ModQV Score vs. Coverage by Motif", x="Per Strand Coverage", y="Modification QV", colour="Type") +
  facet_wrap(~motif) + guides(size=FALSE)
show(p)
# 
p <- qplot(data=subset(mm, CognateBase!="native"), x=coverage, y=score, colour=CognateBase, xlim=c(0,120), ylim=c(0,200), size=I(1.5)) +
  labs(title="ModQV Score vs. Coverage by Motif", x="Per Strand Coverage", y="Modification QV", colour="Cognate Base") +
  facet_wrap(~motif) + guides(size=FALSE)
show(p)
dev.off()

### LOOK AT WHAT *ISN'T* IN THE GFF FILE
# Plot the score vs coverage information for all motif positions found in the reference, including those that aren't hits.
# To do this, we have to get kinetic information on unmethylated motifs from the modifications.csv file.

allKinData <- csv.get(csv)
mga <- getUnmodifiedMotifKin(mm, allKinData)


###PAGE 5
# Compare the score distributions for motif cognate-base positions to scores in the 'no_motif' set
pdf(paste(pdf_out_path, '5.ModificationQV_distribution_by_Motif.pdf', sep=""), width=11, height=8.5, onefile=T)
d=qplot(data=subset(mga, type!="not_detected"), x=score, colour=motif, geom='density') +
  theme(legend.position = c(0.85,0.75)) +
  labs(title="Score Distributions of Each Detected Motif", x="Motif Score", y="Density", colour="Motif")
show(d)
#
# now with facets (usefull when there's lots of motifs)
d=qplot(data=mga, x=score, geom='freqpoly', colour=type) +
  facet_wrap(~ motif, ncol=2, scales="free_y") +
  labs(title="Score Distributions for Motifs by Type", x="Motif Score", y="Density", colour="Type")
show(d)
#
# now zooming into events with scores > 30, since there are so many low scoring / undetected CCWGG sites.
d=qplot(data=subset(mga, score>30), x=score, geom='freqpoly', colour=type) +
  facet_wrap(~ motif, ncol=2, scales="free_y") +
  labs(title="Score Distributions for Motifs by Type", x="Motif Score", y="Density", colour="Type")
show(d)
dev.off()

###PAGE 6
# Compare the IPD ratio distributions of the found motifs.
pdf(paste(pdf_out_path, '6.IPDratio_Distribution_by_Motif.pdf', sep=""), width=11, height=8.5, onefile=T)
d=qplot(data=subset(mga, type!='Not_Detected'), IPDRatio, colour=motif, geom='density') +
  theme(legend.position = c(0.85,0.75)) +
  labs(title="IPD Ratio Distributions by Motif", x="IPD Ratio", y="Density", colour="Motif")
show(d)
#
# now with facets (usefull when there's lots of motifs)
d=qplot(data=mga, x=IPDRatio, colour=type, geom='freqpoly') +
  facet_wrap(~ motif, ncol=2, scales="free_y") +
  labs(x="IPD Ratio", y="Density", colour="Type")
show(d)
#
# now plotting only events with scores > 30, so that we can see the tail of the CCWGG motif data
d=qplot(data=subset(mga, score > 30), x=IPDRatio, colour=type, geom='freqpoly') +
  facet_wrap(~ motif, ncol=2, scales="free_y") +
  labs(x="IPD Ratio", y="Density", colour="Type")
show(d)
dev.off()


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
write.table(summary, file = paste(pdf_out_path, 'MotifSummary.csv', sep=''), quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)


# Since the modifications.csv is so large, you may want to subset the csv file to only the lines you need
# for making plots + circos plotting.  Circos plotting needs baseline datapoints, so you can temporarily add 'GA'
# to your motif list, then regenerate mm for this purpose. This new csv will load faster if you need to revisit the analysis.
allMotifKin <- subset(allKinData, tpl %in% mm$start)
write.table(allMotifKin, file = 'modificationsSubset.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

# R is a very useful tools for making nearly any type of plot.  The 'names' command reveals
# all the information held in the dataframe 'mga' that you may be interested in further examining.
# Check out the ggplot online manual for documentation and additional plotting options:
# docs.ggplot2.org/current/

names(mga)

#### MAKING CIRCOS PLOTS TO VISUALIZE YOUR HITS AND MOTIFS

# Create a karyotype track using information from the reference fasta and the gff header.
makeKaryotypeFile(refPath, gff, circosPrefix)

# Make circos highlights tracks showing the genome positions of each motif,
# both modified and unmodified.
makeCircosMotifHighlights(mm, circosPrefix)

# Make circos basemod score tracks (the red spikes) for each set of motif-associated hits.
makeCircosTracks(mm, circosPrefix)

# Take data points from the base modification csv file to draw a baseline in the circos plot. Choose a 'motif'
# that will be common enough to define a baseline but not so common as to take forever to run. For a plasmid,
# one might choose 'A'; for a bacterial genome 'GA' would be appropriate.
baselineMotif <- c("GT")
makeCircosBaseline(circosPrefix, baselineMotif, refPath, gff, allKinData)
mergeMotifAndBaselineData(circosPrefix, motifs, baselineMotif)

# From here, you can generate a circos config file, which will include pointers to all the files you've just written. 
makeCircosConfFiles(motifs, circosPrefix, confModel, baselineMotif)

# To generate the plot, exit R and do the following from the bash prompt in the directory with your
# config files (assuming you have circos-0.60 installed, which is NOT the most current version):
# $ /tools/circos/bin/circos -conf (word).conf

# From here, you have all the input files needed to assemble the circos config file. 
makeCircosConfFiles(motifs, circosPrefix, confModel, baselineMotif)

# To generate the plot, from the bash prompt (not in R!):
# $ circos -conf (word).conf

