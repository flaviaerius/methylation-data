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

diffCpGann.d = annotateWithFeatureFlank(as(myDiff25p.d,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")
# 50%

# diffCpGann.50.d = annotateWithFeatureFlank(as(myDiff50p.d,"GRanges"),
#                                         cpg.obj$CpGi,cpg.obj$shores,
#                                         feature.name="CpGi",flank.name="shores")


#### Regional analysis ####

## Promoter/Exon/Intron ##

# summarize methylation information over a set of defined regions such as promoters or CpG islands
# the function below does that and outputs a methylRaw or methylRawList object depending on the input

promoters.d = regionCounts(myobj.d,gene.obj$promoters)

head(promoters.d[[1]])

exons.d = regionCounts(myobj.d,gene.obj$exons)

head(exons.d[[1]])

introns.d = regionCounts(myobj.d,gene.obj$introns)

head(introns.d[[1]])

### Here I decided to store information about the promoters, exons and introns of each sample, so I can 
# have these subsets and futurely do a diff meth analysis of these data separately


# Annotate differentially methylated CpGs with promoter/exon/intron using annotation data

diffAnn.d = annotateWithGeneParts(as(myDiff25p.d,"GRanges"),gene.obj)

# update 12.15.2017 --> save plots separated hyper and hypomethylation
# here I had to change the attribute dbpath to eg attr(myDiff25p.hyper.a, "dbpath") <- "/Users/flaviaerius/Documents/methylation_data/cpg_ovd/filtered/methylDiff_1978669e66b8 _type .txt.bgz"
# because I changed the name of the directories (same to hypo :
# attr(myDiff25p.hypo.a, "dbpath") <- "/Users/flaviaerius/Documents/methylation_data/cpg_ovd/filtered/methylDiff_1978669e66b8 _type .txt.bgz")


diffCpGann.hyper.d = annotateWithFeatureFlank(as(myDiff25p.hyper.d,"GRanges"),
                                              cpg.obj$CpGi,cpg.obj$shores,
                                              feature.name="CpGi",flank.name="shores")

diffCpGann.hypo.d = annotateWithFeatureFlank(as(myDiff25p.hypo.d,"GRanges"),
                                             cpg.obj$CpGi,cpg.obj$shores,
                                             feature.name="CpGi",flank.name="shores")

# To get the association with distance to TSS and nearest gene name, we use the getAssociationWithTSS
# function from package genomation

# get association differentially methylated all

assocTSS.d = getAssociationWithTSS(diffAnn.d)

# get association hypermethylated

diffAnn.hyper.d = annotateWithGeneParts(as(myDiff25p.hyper.d,"GRanges"),gene.obj)

# capture output

assocTSS.hyper.d = getAssociationWithTSS(diffAnn.hyper.d)

# get association hypomethylated

diffAnn.hypo.d = annotateWithGeneParts(as(myDiff25p.hypo.d,"GRanges"),gene.obj)

# capture output

assocTSS.hypo.d = getAssociationWithTSS(diffAnn.hypo.d)

# get stats about percentage of DMR that overlap with intron/exon/promoters

getTargetAnnotationStats(diffAnn.d, percentage=TRUE, precedence=TRUE)


# # get and save the plot
# 
 jpeg("target-annotation-perc-ovd.d.jpg")
 plotTargetAnnotation(diffAnn.d, precedence=TRUE, col = c("orangered3", "steelblue4", "gold", "antiquewhite" ),
                      main="Annotation of DMR") # perhaps change the main here
 dev.off()
 
 # update 12.15.2017 hyper and hypo --> I have to tell R that this function is from genomation because it is deprecated in methylKit  
 # hyper
 
 png("target-annotation-perc-ovd-hyper.d.png")
 genomation::plotTargetAnnotation(diffAnn.hyper.d, precedence=TRUE, col = c("orangered3", "steelblue4", "gold", "antiquewhite" ),
                                  main="Annotation of DMR - hypermethylated") # perhaps change the main here
 dev.off()
 
 # hypo
 
 png("target-annotation-perc-ovd-hypo.d.png")
 genomation::plotTargetAnnotation(diffAnn.hypo.d, precedence=TRUE, col = c("orangered3", "steelblue4", "gold", "antiquewhite" ),
                                  main="Annotation of DMR - hypomethylated") # perhaps change the main here
 dev.off()

# # get and save plot of CpG islands and shores annotation
# 
 jpeg("cpgi-shores-other-ovd.d.jpg")
 plotTargetAnnotation(diffCpGann.d,col=c("blue","gray","white"),
                      main="CpG enrichment of DMR") #perhaps change the main here too
 dev.off()
 # update 12.15.2017 hyper and hypo --> I have to tell R that this function is from genomation because it is deprecated in methylKit 
 # hyper
 
 png("cpgi-shores-other-ovd-hyper.d.png")
 genomation::plotTargetAnnotation(diffCpGann.hyper.d,col=c("blue","gray","white"),
                                  main="CpG enrichment of DMR - hypermethylated") #perhaps change the main here too
 dev.off()
 
 # hypo
 
 png("cpgi-shores-other-ovd-hypo.d.png")
 genomation::plotTargetAnnotation(diffCpGann.hypo.d,col=c("blue","gray","white"),
                                  main="CpG enrichment of DMR - hypomethylated") #perhaps change the main here too
 dev.off()

 write.table(assocTSS.hyper.d, "assocTSS_hyper_d_ovd.txt", quote = F,  sep = "\t")
 
 write.table(assocTSS.hypo.d, "assocTSS_hypo_d_ovd.txt", quote = F,  sep = "\t")