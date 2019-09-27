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

diffCpGann.c = annotateWithFeatureFlank(as(myDiff25p.c,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")

# update 12.15.2017 --> save plots separated hyper and hypomethylation


diffCpGann.hyper.c = annotateWithFeatureFlank(as(myDiff25p.hyper.c,"GRanges"),
                                              cpg.obj$CpGi,cpg.obj$shores,
                                              feature.name="CpGi",flank.name="shores")

diffCpGann.hypo.c = annotateWithFeatureFlank(as(myDiff25p.hypo.c,"GRanges"),
                                             cpg.obj$CpGi,cpg.obj$shores,
                                             feature.name="CpGi",flank.name="shores")

# 50%

# diffCpGann.50.c = annotateWithFeatureFlank(as(myDiff50p.c,"GRanges"),
#                                         cpg.obj$CpGi,cpg.obj$shores,
#                                         feature.name="CpGi",flank.name="shores")


#### Regional analysis ####

## Promoter/Exon/Intron ##

# summarize methylation information over a set of defined regions such as promoters or CpG islands
# the function below does that and outputs a methylRaw or methylRawList object depending on the input

promoters.c = regionCounts(myobj.c,gene.obj$promoters)

head(promoters.c[[1]])

exons.c = regionCounts(myobj.c,gene.obj$exons)

head(exons.c[[1]])

introns.c = regionCounts(myobj.c,gene.obj$introns)

head(introns.c[[1]])

### Here I decided to store information about the promoters, exons and introns of each sample, so I can 
# have these subsets and futurely do a diff meth analysis of these data separately


# Annotate differentially methylated CpGs with promoter/exon/intron using annotation data

diffAnn.c = annotateWithGeneParts(as(myDiff25p.c,"GRanges"),gene.obj)

# To get the association with distance to TSS and nearest gene name, we use the getAssociationWithTSS
# function from package genomation

# get association differentially methylated all

assocTSS.c = getAssociationWithTSS(diffAnn.c)

# get association hypermethylated

diffAnn.hyper.c = annotateWithGeneParts(as(myDiff25p.hyper.c,"GRanges"),gene.obj)

# capture output

assocTSS.hyper.c = getAssociationWithTSS(diffAnn.hyper.c)

# get association hypomethylated

diffAnn.hypo.c = annotateWithGeneParts(as(myDiff25p.hypo.c,"GRanges"),gene.obj)

# capture output

assocTSS.hypo.c = getAssociationWithTSS(diffAnn.hypo.c)

# get stats about percentage of DMR that overlap with intron/exon/promoters

getTargetAnnotationStats(diffAnn.c, percentage=TRUE, precedence=TRUE)


# # get and save the plot
# 
 jpeg("target-annotation-perc-ovd.c.jpg")
 plotTargetAnnotation(diffAnn.c, precedence=TRUE, col = c("orangered3", "steelblue4", "gold", "antiquewhite" ),
                      main="Annotation of DMR") # perhaps change the main here
 dev.off()
 
 # update 12.15.2017 hyper and hypo --> I have to tell R that this function is from genomation because it is deprecated in methylKit  
 # hyper
 
 png("target-annotation-perc-ovd-hyper.c.png")
 genomation::plotTargetAnnotation(diffAnn.hyper.c, precedence=TRUE, col = c("orangered3", "steelblue4", "gold", "antiquewhite" ),
                                  main="Annotation of DMR - hypermethylated") # perhaps change the main here
 dev.off()
 
 # hypo
 
 png("target-annotation-perc-ovd-hypo.c.png")
 genomation::plotTargetAnnotation(diffAnn.hypo.c, precedence=TRUE, col = c("orangered3", "steelblue4", "gold", "antiquewhite" ),
                                  main="Annotation of DMR - hypomethylated") # perhaps change the main here
 dev.off()

# # get and save plot of CpG islands and shores annotation
 
 jpeg("cpgi-shores-other-ovd.c.jpg")
 plotTargetAnnotation(diffCpGann.c,col=c("blue","gray","white"),
                      main="CpG enrichment of DMR") #perhaps change the main here too
 dev.off()

 # update 12.15.2017 hyper and hypo --> I have to tell R that this function is from genomation because it is deprecated in methylKit 
 # hyper
 
 png("cpgi-shores-other-ovd-hyper.c.png")
 genomation::plotTargetAnnotation(diffCpGann.hyper.c,col=c("blue","gray","white"),
                                  main="CpG enrichment of DMR - hypermethylated") #perhaps change the main here too
 dev.off()
 
 # hypo
 
 png("cpgi-shores-other-ovd-hypo.c.png")
 genomation::plotTargetAnnotation(diffCpGann.hypo.c,col=c("blue","gray","white"),
                                  main="CpG enrichment of DMR - hypomethylated") #perhaps change the main here too
 dev.off()
 
 write.table(assocTSS.hyper.c, "assocTSS_hyper_c_ovd.txt", quote = F,  sep = "\t")
 
 write.table(assocTSS.hypo.c, "assocTSS_hypo_c_ovd.txt", quote = F,  sep = "\t")