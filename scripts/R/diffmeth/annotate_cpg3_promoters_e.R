library(org.Mm.eg.db)

# from directory cpg_ovd
# read counts with ids

prom1000.counts.hyper <- read.table("promoters/cpg3/cpg3_count_prom1000_assocTSS_hyper_e_ovd.txt", header = F)
prom1000.counts.hypo <- read.table("promoters/cpg3/cpg3_count_prom1000_assocTSS_hypo_e_ovd.txt", header = F)
prom500.counts.hyper <- read.table("promoters/cpg3/cpg3_count_prom500_assocTSS_hyper_e_ovd.txt", header = F)
prom500.counts.hypo <- read.table("promoters/cpg3/cpg3_count_prom500_assocTSS_hypo_e_ovd.txt", header = F)

#select id

prom500.hyper.inputj <- as.character(prom500.counts.hyper$V2)
prom500.hypo.inputj <- as.character(prom500.counts.hypo$V2)

prom1000.hyper.inputj <- as.character(prom1000.counts.hyper$V2)
prom1000.hypo.inputj <- as.character(prom1000.counts.hypo$V2)

# 

prom500.hyper.refseq.nodot <- gsub(pattern = "\\..", "", prom500.hyper.inputj, ignore.case = T) 
prom500.hypo.refseq.nodot <- gsub(pattern = "\\..", "", prom500.hypo.inputj, ignore.case = T) 

prom1000.hyper.refseq.nodot <- gsub(pattern = "\\..", "", prom1000.hyper.inputj, ignore.case = T) 
prom1000.hypo.refseq.nodot <- gsub(pattern = "\\..", "", prom1000.hypo.inputj, ignore.case = T) 

#

prom500.hyper.input <- grep("X", prom500.hyper.refseq.nodot, invert = T, value = T)
prom500.hypo.input <- grep("X", prom500.hypo.refseq.nodot, invert = T, value = T)

prom1000.hyper.input <- grep("X", prom1000.hyper.refseq.nodot, invert = T, value = T)
prom1000.hypo.input <- grep("X", prom1000.hypo.refseq.nodot, invert = T, value = T)


# construct gene annotations

anno3.p500.hyper <- as.data.frame(mapIds(org.Mm.eg.db, keys = prom500.hyper.input, keytype = "REFSEQ",
                                        column = "SYMBOL", multiVals = "first" ))
anno3.p500.hypo <- as.data.frame(mapIds(org.Mm.eg.db, keys = prom500.hypo.input, keytype = "REFSEQ",
                                       column = "SYMBOL", multiVals = "first" ))

anno3.p1000.hyper <- as.data.frame(mapIds(org.Mm.eg.db, keys = prom1000.hyper.input, keytype = "REFSEQ",
                                         column = "SYMBOL", multiVals = "first" ))
anno3.p1000.hypo <- as.data.frame(mapIds(org.Mm.eg.db, keys = prom1000.hypo.input, keytype = "REFSEQ",
                                        column = "SYMBOL", multiVals = "first" ))

# save databases to build gene sets in gsea

write.table(anno3.p500.hyper, "anno_p500_hyper_e.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(anno3.p500.hypo, "anno_p500_hypo_e.txt", sep = "\t", quote = F, row.names = F, col.names = F)

write.table(anno3.p1000.hyper, "anno_p1000_hyper_e.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(anno3.p1000.hypo, "anno_p1000_hypo_e.txt", sep = "\t", quote = F, row.names = F, col.names = F)

save.image("e.RData")
