#### Comparative analysis ####

# First we have to merge samples by base covered

# The unite() function will return a methylBase object
# that contains methylation information for regions/bases
# covered in all samples.

 meth.b = unite(filtered.myobj.b, destrand = FALSE)

# take a look at the meth object

    head(meth.b)

# create object meth.b as dataframe to try to bind to others meth.n to do a proper correlation graph    
    
# Sample Correlation:
# To do the correlation of samples, use getCorrelation()
# if you choose plot = TRUE it returns both scatterplots and
# correlation coefficients.
# the default is plot = FALSE, so it returns a correlation matrix 
# 
#     jpeg("correlation.b.jpg")
#     getCorrelation(meth.b, plot = TRUE)
#     dev.off()

## Clustering Samples
    
# get a cluster of the samples based in the methylation similarities
# method of this cluster is ward.D(see hclust() function), using pearson as correlation distance

#### WARNING!!!!!!####
# here i will have to investigate how to do it using all 4 samples

    # jpeg("cluster.b.jpg")
    # clusterSamples(meth.b, dist = "correlation", method = "ward", plot = TRUE) 
    # dev.off()

# When I put # # it means it was stracted straight from the user guide of methylKit:

# # Setting the plot=FALSE will return a dendrogram object which can be manipulated
# # by users or fed in to other user functions that can work with dendrograms
# --> maybe the described above can work for me to bind the data of 4 lineages,
# so I'm going to apply this variable construction:

    # hc.b = clusterSamples(meth.b, dist="correlation", method="ward", plot=FALSE)

# this way I am storing this information in the variable hc.b

# PCA Analysis

# PCA analysis will find the pricipal components that explain variance of the samples.
# In this case, with biological replicates, it shows the strength of 
# Scree plot: A scree plot graphs the amount of variation explained by each component

    # jpeg("screeplot.b.jpg")
    # PCASamples(meth.b, screeplot=TRUE)
    # dev.off()

# Scatter plot of PCA analysis

    # jpeg("scatterplot.pca.b.jpg")
    # PCASamples(meth.b)
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

  myDiff.b = calculateDiffMeth(meth.b, overdispersion="MN")

# the myDiff.b variable will contain all the q-values for the selected hypothesis testing

# I would take a look at myDiff.b data with
# head(myDiff.b)

# Save myDiff.b as a dataframe in the folder
   
# out <- capture.output(as.data.frame(myDiff.b))
#  
# cat(out, file="diff-meth-all.b.bsv", sep=",", append=TRUE, fill = 5)
    
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

myDiff25p.hyper.b = getMethylDiff(myDiff.b, difference = 25, qvalue = 0.01, type = "hyper", save.db = F)

# out <- capture.output(as.data.frame(myDiff25p.hyper.b))
# options(max.print = .Machine$integer.max) 
# cat(out, file="diff-meth-hyper-25p-q001.b.txt", sep=" ", append=TRUE, fill = 5)
    
# get hypomethylated bases

# difference of 25%

myDiff25p.hypo.b = getMethylDiff(myDiff.b, difference = 25, qvalue = 0.01, type = "hypo", save.db = F)

# out <- capture.output(as.data.frame(myDiff25p.hypo.b))
# options(max.print = .Machine$integer.max)
# cat(out, file="diff-meth-hypo-25p-q001.b.txt", sep=" ", append=TRUE, fill = 5)

# difference of 50%

# myDiff50p.hypo.b = getMethylDiff(myDiff.b, difference = 50, qvalue = 0.01, type = "hypo", save.db = T)

# out <- capture.output(as.data.frame(myDiff50p.hypo.b))
# options(max.print = .Machine$integer.max)
# cat(out, file="diff-meth-hypo-50p-q001.b.txt", sep=" ", append=TRUE, fill = 5)


# get all differentially methylated bases

# 25%

myDiff25p.b = getMethylDiff(myDiff.b, difference = 25, qvalue = 0.01, save.db = F)
# out <- capture.output(as.data.frame(myDiff25p.b))
# options(max.print = .Machine$integer.max)
# cat(out, file="diff-meth-25p-q001.b.txt", sep=" ", append=TRUE, fill = 5)

# 50%

# myDiff50p.b = getMethylDiff(myDiff.b, difference = 50, qvalue = 0.01)

# out <- capture.output(as.data.frame(myDiff50p.b))
# options(max.print = .Machine$integer.max)
# cat(out, file="diff-meth-50p-q001.b.txt", sep=" ", append=TRUE, fill = 5)


# we can check the distribution of hypo/hypermethylated bases per chromosome:

#    outPerChr.b <- capture.output(as.data.frame(
#    diffMethPerChr(myDiff.b ,plot=FALSE,qvalue.cutoff = 0.01, meth.cutoff = 25)
#    ))
#    options(max.print = .Machine$integer.max)
#    cat(outPerChr.b, file="output-diff-meth-per-chr.b.txt", sep=" ", fill = 5)

# get the plot too

    jpeg("diff-meth-per-chr-ovd.b.jpg")
    diffMethPerChr(myDiff.b, plot = TRUE, qvalue.cutoff = 0.01, meth.cutoff = 25)
    dev.off()



