## Analysis of DNA methylation data obtained from ERRBS ##
## R version 3.4.0 ##
## New York, NY, December 2nd of 2017 ##
## MasonLab - Flavia E. Rius ##
# Flavia's macbook pro

#### Ma vs 4C - A ####

# install methylkit

source("https://bioconductor.org/biocLite.R")
biocLite("methylKit")

# Warnings:
# 1: package ‘methylKit’ was built under R version 3.4.2 
# 2: package ‘GenomicRanges’ was built under R version 3.4.2 
# 3: package ‘BiocGenerics’ was built under R version 3.4.2 
# 4: package ‘S4Vectors’ was built under R version 3.4.2 
# 5: package ‘IRanges’ was built under R version 3.4.2 
# 6: package ‘GenomeInfoDb’ was built under R version 3.4.2 

# load methylKit

library("methylKit")

# set your working directory as the one where you have the files that
# you want to work with

# setwd("errbs/methylcall10x/")

# 1st step: reading the methylation call files

# make a list with the files
# since I have 4 conditions, I will have to analyze it in pairs:

#   1 - Ma vs 4C - A
#   2 - 4C vs 4C11minus - B
#   3 - 4C11minus vs 4C11plus - C
#   4 - Ma vs 4C11plus - D
#   5 - 4C vs 4C11plus - E
#   6 - Ma vs 4C11minus - F

# so, as in the RNA-seq analysis, I will do a script for each comparison
# this one is for Ma and 4C

# I was trying to do in a way selecting the patterns:
# file.list.chg = list(c(list.files(pattern="^chg.M"),list.files(pattern="^chg.4C_")))
# But then an error has returned: 
# "Error in methRead(file.list.cpg, sample.id = list("Ma_1", "Ma_2", "Ma_3",  : 
# length of 'location'  and 'name' should be same
# Then I checked the object length:
# >length("file.list.cpg")
# [1] 1
# So, my object had length = 1 and I actually have 6 files. So I had to change the way of 
# providing the files' location, and that is why I did the conventional way:

file.list.cpg.a = list("cpg.Ma_1.mincov10.txt", "cpg.Ma_2.mincov10.txt", "cpg.Ma_3.mincov10.txt",
                     "cpg.4C_1.mincov10.txt", "cpg.4C_2.mincov10.txt", "cpg.4C_3.mincov10.txt")

# file.list.chg.a = list("chg.Ma_1.mincov10.txt", "chg.Ma_2.mincov10.txt", "chg.Ma_3.mincov10.txt", 
#                      "chg.4C_1.mincov10.txt", "chg.4C_2.mincov10.txt", "chg.4C_3.mincov10.txt")

# file.list.chh.a = list("chh.Ma_1.mincov10.txt", "chh.Ma_2.mincov10.txt", "chh.Ma_3.mincov10.txt", 
#                      "chh.4C_1.mincov10.txt", "chh.4C_2.mincov10.txt", "chh.4C_3.mincov10.txt")

# read the files to a methylRawList object: myobj

# specify which file you are reading (cpg, chg or chh) and put it in the name of the 
# methylRawList object

#1 - cpg


myobj.a = methRead(file.list.cpg.a,
               sample.id=list("Ma_1","Ma_2", "Ma_3", "4C_1","4C_2", "4C_3"),
               assembly="mm10",
               treatment=c(0,0,0,1,1,1),
               context="CpG"
)

#2 - chg
# myobj.chg.a = methRead(file.list.chg.a,
#                      sample.id=list("Ma_1","Ma_2", "Ma_3", "4C_1","4C_2", "4C_3"),
#                      assembly="mm10",
#                      treatment=c(0,0,0,1,1,1),
#                      context="CpH"
# )

#3 - chh
# myobj.chh.a = methRead(file.list.chh.a,
#                      sample.id=list("Ma_1","Ma_2", "Ma_3", "4C_1","4C_2", "4C_3"),
#                      assembly="mm10",
#                      treatment=c(0,0,0,1,1,1),
#                      context="CHH"
# )

# the objects generated are methylRawList, and appear like this:

# [[1]]
# chr     start       end strand coverage  numCs numTs
# 1        chr1   3020795   3020795      -       20      0    20
# 2        chr1   3020891   3020891      +       42     16    26
# 3        chr1   3020843   3020843      -       20      1    19
# 4        chr1   3020972   3020972      -       26      0    26
# 5        chr1   3020988   3020988      -       26     24     2
# 6        chr1   3020815   3020815      -       20      2    18



# get descriptive statistics of object of samples comparison:
# the statistics is about each sample, and that is why we need to select the column that represent 
# the sample we want to get the stats from.
# the column number is in the order we determined at first place.

# percent of CpG methylation vs frequency 

# summary statistics

# for (n in c(1:6) )
# { 
#   out <- capture.output(getMethylationStats(myobj.cpg.a[[n]],plot=FALSE,both.strands=FALSE)) 
#     
#     cat(n, out, file="output-meth-stats-ma-4c.txt", sep=" ", append=TRUE, fill = 5)
# }

# plots 

# for (n in c(1:6))
# {
#   mypath <- file.path("~/Documents/errbs/methylcall10X", 
#                       paste( "methy-stats-", n, ".jpeg", sep=""))
#   jpeg(file=mypath)
#   getMethylationStats(myobj.cpg.[[n]],plot=TRUE,both.strands=FALSE)
#   dev.off()
# }


# coverage 

# summary statistics

# for (n in c(1:6) )
# { 
#   out <- capture.output(getCoverageStats(myobj.cpg.a[[n]],plot=FALSE,both.strands=FALSE)) 
#   
#   cat(n, out, file="output-coverage-stats-ma-4c.txt", sep=" ", append=TRUE, fill = 5)
# }

# plots

# for (n in c(1:6))
# {
#   mypath.a <- file.path("~/Documents/errbs/methylcall10X/graphs", 
#                       paste( "coverage-stats-ma-4c-", n, ".jpeg", sep=""))
#   jpeg(file=mypath.a)
#   getCoverageStats(myobj.cpg.a[[n]],plot=TRUE,both.strands=FALSE)
#   dev.off()
# }

####### for CHG

# percent of CHG methylation vs frequency 

# summary statistics

# for (n in c(1:6) )
# { 
#   out <- capture.output(getMethylationStats(myobj.chg.a[[n]],plot=FALSE,both.strands=FALSE)) 
#   
#   cat(n, out, file="output-meth-stats-ma-4c-chg.txt", sep=" ", append=TRUE, fill = 5)
# }

# plots 

# for (n in c(1:6))
# {
#   mypath.a <- file.path("~/Documents/errbs/methylcall10X/graphs", 
#                       paste( "meth-stats-ma-4c-chg-", n, ".jpeg", sep=""))
#   jpeg(file=mypath.a)
#   getMethylationStats(myobj.chg.a[[n]],plot=TRUE,both.strands=FALSE)
#   dev.off()
# }


# coverage 

# summary statistics

# for (n in c(1:6) )
# { 
#   out <- capture.output(getCoverageStats(myobj.chg.a[[n]],plot=FALSE,both.strands=FALSE)) 
#   
#   cat(n, out, file="output-coverage-stats-ma-4c-chg.txt", sep=" ", append=TRUE, fill = 5)
# }

# plots

## for (n in c(1:6))
# {
#   mypath.a <- file.path("~/Documents/errbs/methylcall10X/graphs", 
#                       paste( "coverage-stats-chg-", n, ".jpeg", sep=""))
#   jpeg(file=mypath.a)
#   getCoverageStats(myobj.chg.a[[n]],plot=TRUE,both.strands=FALSE)
#   dev.off()
# }

## for CHH

# percent of CpG methylation vs frequency 

# summary statistics

# for (n in c(1:6) )
# { 
#   out <- capture.output(getMethylationStats(myobj.chh.a[[n]],plot=FALSE,both.strands=FALSE)) 
#   
#   cat(n, out, file="output-meth-stats-ma-4c-chh.txt", sep=" ", append=TRUE, fill = 5)
# }

# plots 

# for (n in c(1:6))
# {
#   mypath <- file.path("~/Documents/errbs/methylcall10X/graphs", 
#                       paste( "methy-stats-chh-", n, ".jpeg", sep=""))
#   jpeg(file=mypath)
#   getMethylationStats(myobj.chh.a[[n]],plot=TRUE,both.strands=FALSE)
#   dev.off()
# }


# coverage 

# summary statistics

#for (n in c(1:6) )
#{
#  out <- capture.output(getCoverageStats(myobj.chh.a[[n]],plot=FALSE,both.strands=FALSE))

#  cat(n, out, file="output-coverage-stats-ma-4c-chh.txt", sep=" ", append=TRUE, fill = 5)
#}

# plots

# for (n in c(1:6))
# {
# mypath <- file.path("~/Documents/errbs/methylcall10X/graphs",
#                     paste( "coverage-stats-chh-", n, ".jpeg", sep=""))
# jpeg(file=mypath)
# getCoverageStats(myobj.chh.a[[n]],plot=TRUE,both.strands=FALSE)
# dev.off()
# }

# create directory to send the filtered files

dir.create("filtered")

# filter coverage by excluding the CpG sites with higher coverage than the 99,9 percentile because
# this can be a PCR bias.
# we don't need to filter the lower than 10X coverage because the files were already filtered by
# the EpiCore.

filtered.myobj.a = filterByCoverage(myobj.a, hi.perc = 99.9
                                     , save.db = TRUE, dbdir = "filtered" 
                                    )
