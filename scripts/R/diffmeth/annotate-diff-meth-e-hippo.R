### Annotating differentially methylated bases ###

# First, load package genomation


library(genomation)

# we will annotate the differentially methylated regions according to position in the genome
# we will use mm10 bed file as reference (from refseq) (i got it from:UCSC)

# read the gene BED file:

gene.obj = readTranscriptFeatures("mm10.refseqCurated.bed")


# read the shores and flanking regions and name the flanks as shores (here set as 2000bp of flanking regions) 
# and CpG islands as CpGi

cpg.obj = readFeatureFlank("cgi.mm10.bed",
                         feature.flank.name=c("CpGi","shores"))


# convert methylDiff object to GRanges and annotate

diffCpGann.e = annotateWithFeatureFlank(as(myDiff25p.e,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")

# update 12.15.2017 --> save plots separated hyper and hypomethylation
# I didn~t have to change the dbpath here because I created the myDiff25p variables without saving in tabix files.

diffCpGann.hyper.e = annotateWithFeatureFlank(as(myDiff25p.hyper.e,"GRanges"),
                                              cpg.obj$CpGi,cpg.obj$shores,
                                              feature.name="CpGi",flank.name="shores")

diffCpGann.hypo.e = annotateWithFeatureFlank(as(myDiff25p.hypo.e,"GRanges"),
                                             cpg.obj$CpGi,cpg.obj$shores,
                                             feature.name="CpGi",flank.name="shores")

# 50%

# diffCpGann.50.e = annotateWithFeatureFlank(as(myDiff50p.e,"GRanges"),
#                                         cpg.obj$CpGi,cpg.obj$shores,
#                                         feature.name="CpGi",flank.name="shores")


#### Regional analysis ####

## Promoter/Exon/Intron ##

# summarize methylation information over a set of defined regions such as promoters or CpG islands
# the function below does that and outputs a methylRaw or methylRawList object depending on the input

promoters.e = regionCounts(myobj.e,gene.obj$promoters)

head(promoters.e[[1]])

exons.e = regionCounts(myobj.e,gene.obj$exons)

head(exons.e[[1]])

introns.e = regionCounts(myobj.e,gene.obj$introns)

head(introns.e[[1]])

### Here I decided to store information about the promoters, exons and introns of each sample, so I can 
# have these subsets and futurely do a diff meth analysis of these data separately


# Annotate differentially methylated CpGs with promoter/exon/intron using annotation data

diffAnn.e = annotateWithGeneParts(as(myDiff25p.e,"GRanges"),gene.obj)

# To get the association with distance to TSS and nearest gene name, we use the getAssociationWithTSS
# function from package genomation

# get association differentially methylated all

assocTSS.e = getAssociationWithTSS(diffAnn.e)

# get association hypermethylated

diffAnn.hyper.e = annotateWithGeneParts(as(myDiff25p.hyper.e,"GRanges"),gene.obj)

# capture output

assocTSS.hyper.e = getAssociationWithTSS(diffAnn.hyper.e)

# get association hypomethylated

diffAnn.hypo.e = annotateWithGeneParts(as(myDiff25p.hypo.e,"GRanges"),gene.obj)

# capture output

assocTSS.hypo.e = getAssociationWithTSS(diffAnn.hypo.e)

# get stats about percentage of DMR that overlap with intron/exon/promoters

getTargetAnnotationStats(diffAnn.e, percentage=TRUE, precedence=TRUE)


# # get and save the plot
# 
 jpeg("target-annotation-perc-ovd.e.jpg")
 plotTargetAnnotation(diffAnn.e, precedence=TRUE, col = c("orangered3", "steelblue4", "gold", "antiquewhite" ),
                      main="Annotation of DMR") # perhaps change the main here
 dev.off()
 
 # update 12.15.2017 hyper and hypo --> I have to tell R that this function is from genomation because it is deprecated in methylKit  
 # hyper
 
 png("target-annotation-perc-ovd-hyper.e.png")
 genomation::plotTargetAnnotation(diffAnn.hyper.e, precedence=TRUE, col = c("orangered3", "steelblue4", "gold", "antiquewhite" ),
                                  main="Annotation of DMR - hypermethylated") # perhaps change the main here
 dev.off()
 
 # hypo
 
 png("target-annotation-perc-ovd-hypo.e.png")
 genomation::plotTargetAnnotation(diffAnn.hypo.e, precedence=TRUE, col = c("orangered3", "steelblue4", "gold", "antiquewhite" ),
                                  main="Annotation of DMR - hypomethylated") # perhaps change the main here
 dev.off()
 

# # get and save plot of CpG islands and shores annotation
# 
 jpeg("cpgi-shores-other-ovd.e.jpg")
 plotTargetAnnotation(diffCpGann.e,col=c("blue","gray","white"),
                      main="CpG enrichment of DMR") #perhaps change the main here too
 dev.off()
 
 # update 12.15.2017 hyper and hypo --> I have to tell R that this function is from genomation because it is deprecated in methylKit 
 # hyper
 
 png("cpgi-shores-other-ovd-hyper.e.png")
 genomation::plotTargetAnnotation(diffCpGann.hyper.e,col=c("blue","gray","white"),
                                  main="CpG enrichment of DMR - hypermethylated") #perhaps change the main here too
 dev.off()
 
 # hypo
 
 png("cpgi-shores-other-ovd-hypo.e.png")
 genomation::plotTargetAnnotation(diffCpGann.hypo.e,col=c("blue","gray","white"),
                                  main="CpG enrichment of DMR - hypomethylated") #perhaps change the main here too
 dev.off()

 write.table(assocTSS.hyper.e, "assocTSS_hyper_e_ovd.txt", quote = F,  sep = "\t")
 
 write.table(assocTSS.hypo.e, "assocTSS_hypo_e_ovd.txt", quote = F,  sep = "\t")