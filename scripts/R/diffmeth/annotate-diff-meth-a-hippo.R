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

diffCpGann.a = annotateWithFeatureFlank(as(myDiff25p.a,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")

# update 12.15.2017 --> save plots separated hyper and hypomethylation
# here I had to change the attribute dbpath to attr(myDiff25p.hyper.a, "dbpath") <- "/Users/flaviaerius/Documents/methylation_data/cpg_ovd/filtered/methylDiff_1978669e66b8 _type .txt.bgz"
# because I changed the name of the directories (same to hypo :
# attr(myDiff25p.hypo.a, "dbpath") <- "/Users/flaviaerius/Documents/methylation_data/cpg_ovd/filtered/methylDiff_1978669e66b8 _type .txt.bgz")


diffCpGann.hyper.a = annotateWithFeatureFlank(as(myDiff25p.hyper.a,"GRanges"),
                                        cpg.obj$CpGi,cpg.obj$shores,
                                        feature.name="CpGi",flank.name="shores")

diffCpGann.hypo.a = annotateWithFeatureFlank(as(myDiff25p.hypo.a,"GRanges"),
                                              cpg.obj$CpGi,cpg.obj$shores,
                                              feature.name="CpGi",flank.name="shores")
# 50%

# diffCpGann.50.a = annotateWithFeatureFlank(as(myDiff50p.a,"GRanges"),
#                                         cpg.obj$CpGi,cpg.obj$shores,
#                                         feature.name="CpGi",flank.name="shores")


#### Regional analysis ####

## Promoter/Exon/Intron ##

# summarize methylation information over a set of defined regions such as promoters or CpG islands
# the function below does that and outputs a methylRaw or methylRawList object depending on the input

promoters.a = regionCounts(myobj.a,gene.obj$promoters)

head(promoters.a[[1]])

exons.a = regionCounts(myobj.a,gene.obj$exons)

head(exons.a[[1]])

introns.a = regionCounts(myobj.a,gene.obj$introns)

head(introns.a[[1]])

### Here I decided to store information about the promoters, exons and introns of each sample, so I can 
# have these subsets and futurely do a diff meth analysis of these data separately


# Annotate differentially methylated CpGs with promoter/exon/intron using annotation data

diffAnn.a = annotateWithGeneParts(as(myDiff25p.a,"GRanges"),gene.obj)

# To get the association with distance to TSS and nearest gene name, we use the getAssociationWithTSS
# function from package genomation

# get association differentially methylated all

assocTSS.a = getAssociationWithTSS(diffAnn.a)

# get association hypermethylated

diffAnn.hyper.a = annotateWithGeneParts(as(myDiff25p.hyper.a,"GRanges"),gene.obj)

# capture output

assocTSS.hyper.a = getAssociationWithTSS(diffAnn.hyper.a)

# get association hypomethylated

diffAnn.hypo.a = annotateWithGeneParts(as(myDiff25p.hypo.a,"GRanges"),gene.obj)

# capture output

assocTSS.hypo.a = getAssociationWithTSS(diffAnn.hypo.a)

# get stats about percentage of DMR that overlap with intron/exon/promoters

getTargetAnnotationStats(diffAnn.a, percentage=TRUE, precedence=TRUE)


# # get and save the plot
# 
 jpeg("target-annotation-perc-ovd.a.jpg")
 genomation::plotTargetAnnotation(diffAnn.a, precedence=TRUE, col = c("orangered3", "steelblue4", "gold", "antiquewhite" ),
                      main="Annotation of DMR") # perhaps change the main here
 dev.off()
 # update 12.15.2017 hyper and hypo --> I have to tell R that this function is from genomation because it is deprecated in methylKit  
# hyper
 
 png("target-annotation-perc-ovd-hyper.a.png")
 genomation::plotTargetAnnotation(diffAnn.hyper.a, precedence=TRUE, col = c("orangered3", "steelblue4", "gold", "antiquewhite" ),
                      main="Annotation of DMR - hypermethylated") # perhaps change the main here
 dev.off()
 
# hypo
 
 png("target-annotation-perc-ovd-hypo.a.png")
 genomation::plotTargetAnnotation(diffAnn.hypo.a, precedence=TRUE, col = c("orangered3", "steelblue4", "gold", "antiquewhite" ),
                      main="Annotation of DMR - hypomethylated") # perhaps change the main here
 dev.off()

# # get and save plot of CpG islands and shores annotation
# 
 jpeg("cpgi-shores-other-ovd.a.jpg")
 genomation::plotTargetAnnotation(diffCpGann.a,col=c("blue","gray","white"),
                      main="CpG enrichment of DMR") #perhaps change the main here too
 dev.off()
 
# update 12.15.2017 hyper and hypo --> I have to tell R that this function is from genomation because it is deprecated in methylKit 
 # hyper
 
 png("cpgi-shores-other-ovd-hyper.a.png")
 genomation::plotTargetAnnotation(diffCpGann.hyper.a,col=c("blue","gray","white"),
                                  main="CpG enrichment of DMR - hypermethylated") #perhaps change the main here too
 dev.off()
 
 # hypo

 png("cpgi-shores-other-ovd-hypo.a.png")
 genomation::plotTargetAnnotation(diffCpGann.hypo.a,col=c("blue","gray","white"),
                                  main="CpG enrichment of DMR - hypomethylated") #perhaps change the main here too
 dev.off()

 # save assocTSS
 
 write.table(assocTSS.hyper.a, "assocTSS_hyper_a_ovd.txt", quote = F,  sep = "\t")
 
 write.table(assocTSS.hypo.a, "assocTSS_hypo_a_ovd.txt", quote = F,  sep = "\t")
 