#### Comparative analysis ####

# First we have to merge samples by base covered

# The unite() function will return a methylBase object
# that contains methylation information for regions/bases
# covered in all samples.

 meth.c = unite(filtered.myobj.c, destrand = FALSE)

# take a look at the meth object

    head(meth.c)

# create object meth.c as dataframe to try to bind to others meth.n to do a proper correlation graph    
    
# Sample Correlation:
# To do the correlation of samples, use getCorrelation()
# if you choose plot = TRUE it returns both scatterplots and
# correlation coefficients.
# the default is plot = FALSE, so it returns a correlation matrix 
# 
#     jpeg("correlation.c.jpg")
#     getCorrelation(meth.c, plot = TRUE)
#     dev.off()

## Clustering Samples
    
# get a cluster of the samples based in the methylation similarities
# method of this cluster is ward.D(see hclust() function), using pearson as correlation distance

#### WARNING!!!!!!####
# here i will have to investigate how to do it using all 4 samples

    # jpeg("cluster.c.jpg")
    # clusterSamples(meth.c, dist = "correlation", method = "ward", plot = TRUE) 
    # dev.off()

# When I put # # it means it was stracted straight from the user guide of methylKit:

# # Setting the plot=FALSE will return a dendrogram object which can be manipulated
# # by users or fed in to other user functions that can work with dendrograms
# --> maybe the described above can work for me to bind the data of 4 lineages,
# so I'm going to apply this variable construction:

    # hc.c = clusterSamples(meth.c, dist="correlation", method="ward", plot=FALSE)

# this way I am storing this information in the variable hc.c

# PCA Analysis

# PCA analysis will find the pricipal components that explain variance of the samples.
# In this case, with biological replicates, it shows the strength of 
# Scree plot: A scree plot graphs the amount of variation explained by each component

    # jpeg("screeplot.c.jpg")
    # PCASamples(meth.c, screeplot=TRUE)
    # dev.off()

# Scatter plot of PCA analysis

    # jpeg("scatterplot.pca.c.jpg")
    # PCASamples(meth.c)
    # dev.off()

# batch effect (discuss if there is any in our samples)

# Finding differentially methylated bases

# Depending on the sample size per each set, it is taken logistic regression (small),
# or Fisher's exact test (big) to calculate P-values. If there are replicates, it will
# automatically use logistic regression

# P-values will be adjusted to Q-values using SLIM method (read paper)
# the test for differential methylation will test if a model containing the variable
# "treatment" fits better the data than a model that don't account for the existence
# or not of a treatment difference between the samples.

  myDiff.c = calculateDiffMeth(meth.c, overdispersion="MN")

# the myDiff.c variable will contain all the q-values for the selected hypothesis testing

# I would take a look at myDiff.c data with
# head(myDiff.c)

# Save myDiff.c as a dataframe in the folder
   
# out <- capture.output(as.data.frame(myDiff.c))
#  
# cat(out, file="diff-meth-all.c.csv", sep=",", append=TRUE, fill = 5)
    
# It was saving only the first 142 rows and in the end of the file it was written:
# "reached getOption("max.print"). 1247342 rows ommitted"
# So I searched this error in the web and found out that I needed to change the 
# options(max.print) parameter:
    # options(max.print = .Machine$integer.max)
# with this, I set the max.print to the maximum that the machine can do. 

# I will not save this file because it is very big. I will save the subsets

# next step is select the range of q-values that we want to consider as signifficant
# the protocol selects min qvalue of 0.01 and difference of 25 percent, but I would take a 
# look at the data to see maximum, minimum, mean, median of qvalue and of difference
# between samples. 

####### See if there is any way to do this, like a function "stats" maybe! ########

# get hypermethylated bases
    
# difference of 25%

myDiff25p.hyper.c = getMethylDiff(myDiff.c, difference = 25, qvalue = 0.01, type = "hyper", save.db = F)

# out <- capture.output(as.data.frame(myDiff25p.hyper.c))
# options(max.print = .Machine$integer.max) 
# cat(out, file="diff-meth-hyper-25p-q001.c.txt", sep=" ", append=TRUE, fill = 5)
    
# get hypomethylated bases

# difference of 25%

myDiff25p.hypo.c = getMethylDiff(myDiff.c, difference = 25, qvalue = 0.01, type = "hypo", save.db = F)

# out <- capture.output(as.data.frame(myDiff25p.hypo.c))
# options(max.print = .Machine$integer.max)
# cat(out, file="diff-meth-hypo-25p-q001.c.txt", sep=" ", append=TRUE, fill = 5)

# difference of 50%

# myDiff50p.hypo.c = getMethylDiff(myDiff.c, difference = 50, qvalue = 0.01, type = "hypo", save.db = T)

# out <- capture.output(as.data.frame(myDiff50p.hypo.c))
# options(max.print = .Machine$integer.max)
# cat(out, file="diff-meth-hypo-50p-q001.c.txt", sep=" ", append=TRUE, fill = 5)


# get all differentially methylated bases

# 25%

myDiff25p.c = getMethylDiff(myDiff.c, difference = 25, qvalue = 0.01, save.db = F)

# out <- capture.output(as.data.frame(myDiff25p.c))
# options(max.print = .Machine$integer.max)
# cat(out, file="diff-meth-25p-q001.c.txt", sep=" ", append=TRUE, fill = 5)

# 50%

# myDiff50p.c = getMethylDiff(myDiff.c, difference = 50, qvalue = 0.01)

# out <- capture.output(as.data.frame(myDiff50p.c))
# options(max.print = .Machine$integer.max)
# cat(out, file="diff-meth-50p-q001.c.txt", sep=" ", append=TRUE, fill = 5)


# we can check the distribution of hypo/hypermethylated bases per chromosome:

#    outPerChr.c <- capture.output(as.data.frame(
#    diffMethPerChr(myDiff.c ,plot=FALSE,qvalue.cutoff = 0.01, meth.cutoff = 25)
#    ))
#    options(max.print = .Machine$integer.max)
#    cat(outPerChr.c, file="output-diff-meth-per-chr.c.txt", sep=" ", fill = 5)

# get the plot too

     jpeg("diff-meth-per-chr-ovd.c.jpg")
     diffMethPerChr(myDiff.c, plot = TRUE, qvalue.cutoff = 0.01, meth.cutoff = 25)
     dev.off()



