### Annotating differentially methylated bases ###

# First, load package genomation


library(genomation)

# we will annotate the differentially methylated regions according to position in the genome
# we will use mm10 bed file as reference (from refseq) (i got it from:
# /zenodotus/dat01/mason_lab_scratch_reference/cmlab/GENOMES/ANNOTATION/mm10/refseq/2015-02-09/mm10.refgene.bed)

# read the gene BED file:

gene.obj = readTranscriptFeatures("mm10.refseqCurated.bed")


# read the shores and flanking regions and name the flanks as shores (here set as 2000bp of flanking regions) 
# and CpG islands as CpGi

cpg.obj = readFeatureFlank("cgi.mm10.bed",
                         feature.flank.name=c("CpGi","shores"))


# convert methylDiff object to GRanges and annotate

diffCpGann.b = annotateWithFeatureFlank(as(myDiff25p.b,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")

#hyper e hypo - update 12.16.2017

diffCpGann.hyper.b = annotateWithFeatureFlank(as(myDiff25p.hyper.b,"GRanges"),
                                              cpg.obj$CpGi,cpg.obj$shores,
                                              feature.name="CpGi",flank.name="shores")

diffCpGann.hypo.b = annotateWithFeatureFlank(as(myDiff25p.hypo.b,"GRanges"),
                                             cpg.obj$CpGi,cpg.obj$shores,
                                             feature.name="CpGi",flank.name="shores")

# 50%

# diffCpGann.50.b = annotateWithFeatureFlank(as(myDiff50p.b,"GRanges"),
#                                         cpg.obj$CpGi,cpg.obj$shores,
#                                         feature.name="CpGi",flank.name="shores")


#### Regional analysis ####

## Promoter/Exon/Intron ##

# summarize methylation information over a set of defined regions such as promoters or CpG islands
# the function below does that and outputs a methylRaw or methylRawList object depending on the input

promoters.b = regionCounts(myobj.b,gene.obj$promoters)

head(promoters.b[[1]])

exons.b = regionCounts(myobj.b,gene.obj$exons)

head(exons.b[[1]])

introns.b = regionCounts(myobj.b,gene.obj$introns)

head(introns.b[[1]])

### Here I decided to store information about the promoters, exons and introns of each sample, so I can 
# have these subsets and futurely do a diff meth analysis of these data separately


# Annotate differentially methylated CpGs with promoter/exon/intron using annotation data

diffAnn.b = annotateWithGeneParts(as(myDiff25p.b,"GRanges"),gene.obj)

# To get the association with distance to TSS and nearest gene name, we use the getAssociationWithTSS
# function from package genomation

# get association differentially methylated all

assocTSS.b = getAssociationWithTSS(diffAnn.b)

# get association hypermethylated

diffAnn.hyper.b = annotateWithGeneParts(as(myDiff25p.hyper.b,"GRanges"),gene.obj)

# capture output

assocTSS.hyper.b = getAssociationWithTSS(diffAnn.hyper.b)

# get association hypomethylated

diffAnn.hypo.b = annotateWithGeneParts(as(myDiff25p.hypo.b,"GRanges"),gene.obj)

# capture output

assocTSS.hypo.b = getAssociationWithTSS(diffAnn.hypo.b)

# get stats about percentage of DMR that overlap with intron/exon/promoters

getTargetAnnotationStats(diffAnn.b, percentage=TRUE, precedence=TRUE)


# # get and save the plot
# 
 jpeg("target-annotation-perc-ovd.b.jpg")
 plotTargetAnnotation(diffAnn.b, precedence=TRUE, col = c("orangered3", "steelblue4", "gold", "antiquewhite" ),
                      main="Annotation of DMR") # perhaps change the main here
 dev.off()
 
 # update 12.15.2017 hyper and hypo --> I have to tell R that this function is from genomation because it is deprecated in methylKit  
 # hyper
 
 png("target-annotation-perc-ovd-hyper.b.png")
 genomation::plotTargetAnnotation(diffAnn.hyper.b, precedence=TRUE, col = c("orangered3", "steelblue4", "gold", "antiquewhite" ),
                                  main="Annotation of DMR - hypermethylated") # perhaps change the main here
 dev.off()
 
 # hypo
 
 png("target-annotation-perc-ovd-hypo.b.png")
 genomation::plotTargetAnnotation(diffAnn.hypo.b, precedence=TRUE, col = c("orangered3", "steelblue4", "gold", "antiquewhite" ),
                                  main="Annotation of DMR - hypomethylated") # perhaps change the main here
 dev.off()

# # get and save plot of CpG islands and shores annotation
# 
 jpeg("cpgi-shores-other-ovd.b.jpg")
 plotTargetAnnotation(diffCpGann.b,col=c("blue","gray","white"),
                      main="CpG enrichment of DMR") #perhaps change the main here too
 dev.off()

 # update 12.15.2017 hyper and hypo --> I have to tell R that this function is from genomation because it is deprecated in methylKit 
 # hyper
 
 png("cpgi-shores-other-ovd-hyper.b.png")
 genomation::plotTargetAnnotation(diffCpGann.hyper.b,col=c("blue","gray","white"),
                                  main="CpG enrichment of DMR - hypermethylated") #perhaps change the main here too
 dev.off()
 
 # hypo
 
 png("cpgi-shores-other-ovd-hypo.b.png")
 genomation::plotTargetAnnotation(diffCpGann.hypo.b,col=c("blue","gray","white"),
                                  main="CpG enrichment of DMR - hypomethylated") #perhaps change the main here too
 dev.off()

 write.table(assocTSS.hyper.b, "assocTSS_hyper_b_ovd.txt", quote = F,  sep = "\t")
 
 write.table(assocTSS.hypo.b, "assocTSS_hypo_b_ovd.txt", quote = F,  sep = "\t")