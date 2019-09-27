## Analysis of DNA methylation data obtained from ERRBS ##
## R version 3.4.0 ##
## New York, NY, June 21st - 23rd of 2017 ##
## MasonLab - Flavia E. Rius ##
# Sys.Date()

#### 4C vs 4C11plus - E ####

# install methylkit

# source("https://bioconductor.org/biocLite.R")
# biocLite("methylKit")

# load methylKit

library("methylKit")

# set your working directory as the one where you have the files that
# you want to work with

# setwd("Documents/errbs/methylcall10x/")

# 1st step: reading the methylation call files

# make a list with the files
# since I have 4 conditions, I will have to analyze it in pairs:

#   1 - Ma vs 4C - A
#   2 - 4C vs 4C11minus - B
#   3 - 4C11minus vs 4C11minus - C
#   4 - Ma vs 4C11plus - D
#   5 - 4C vs 4C11plus - E
#   6 - Ma vs 4C11minus - F

# so, as in the RNA-seq analysis, I will do a script for each comparison, identifying each pair with a capital letter from A:F

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

file.list.cpg.e = list("cpg.4C_1.mincov10.txt", "cpg.4C_2.mincov10.txt", "cpg.4C_3.mincov10.txt",
                       "cpg.4C11_plus_1.mincov10.txt", "cpg.4C11_plus_2.mincov10.txt", "cpg.4C11_plus_3.mincov10.txt")

# file.list.chg.e = list("chg.4C_1.mincov10.txt", "chg.4C_2.mincov10.txt", "chg.4C_3.mincov10.txt",
#                        "chg.4C11_plus_1.mincov10.txt", "chg.4C11_plus_2.mincov10.txt", "chg.4C11_plus_3.mincov10.txt")

# file.list.chh.e = list("chh.4C_1.mincov10.txt", "chh.4C_2.mincov10.txt", "chh.4C_3.mincov10.txt",
#                        "chh.4C11_plus_1.mincov10.txt", "chh.4C11_plus_2.mincov10.txt", "chh.4C11_plus_3.mincov10.txt")

# read the files to a methylRawList object: myobj
# specify which file you are reading (cpg, chg or chh) and put it in the name of the 
# methylRawList object

#1 - cpg
myobj.e = methRead(file.list.cpg.e,
               sample.id=list("4C_1","4C_2", "4C_3","4C11_plus_1","4C11_plus_2", "4C11_plus_3"),
               assembly="mm10",
               treatment=c(0,0,0,1,1,1),
               context="CpG"
)

#2 - chg
# myobj.chg.e = methRead(file.list.chg,
#                      sample.id=list("4C_1","4C_2", "4C_3", "4C11_plus_1","4C11_plus_2", 
#                                     "4C11_plus_3"),
#                      assembly="mm10",
#                      treatment=c(0,0,0,1,1,1),
#                      context="CpH"
# )

#3 - chh
# myobj.chh.e = methRead(file.list.chh,
#                      sample.id=list("4C_1","4C_2", "4C_3", "4C11_plus_1","4C11_plus_2", 
#                               "4C11_plus_3"),
#                      assembly="mm10",
#                      treatment=c(0,0,0,1,1,1),
#                      context="CHH")

# the objects generated are methylRawList, and appear like this:

# [[1]]
# chr     start       end strand coverage  numCs numTs
# 1        chr1   3020795   3020795      -       20      0    20
# 2        chr1   3020891   3020891      +       42     16    26
# 3        chr1   3020843   3020843      -       20      1    19
# 4        chr1   3020972   3020972      -       26      0    26
# 5        chr1   3020988   3020988      -       26     24     2
# 6        chr1   3020815   3020815      -       20      2    18


# create directory to send the filtered files

# dir.create("filtered")

# filter coverage by excluding the CpG sites with higher coverage than the 99,9 percentile because
# this can be a PCR bias.
# we don't need to filter the lower than 10X coverage because the files were already filtered by
# the EpiCore.

filtered.myobj.e = filterByCoverage(myobj.e, hi.perc = 99.9, save.db = TRUE, dbdir = "filtered")
