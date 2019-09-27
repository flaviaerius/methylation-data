# 1. Load the saved workspace with the R variables

# load("f.RData")

# we will use the variables prom 1000. and prom 500. assocTSS.hyper and .hypo.f

# 2. Delete the columns named "V1"(come from row number of the primary prom500 .assocTSS file) and "V2" which is the number of the target feature in the file which contains the TSSs

prom500.assocTSS.hyper.f$V1 <- NULL
prom500.assocTSS.hyper.f$V2 <- NULL
colnames(prom500.assocTSS.hyper.f) <- c("dist.to.feature", "feature.name", "feature.strand")

prom500.assocTSS.hypo.f$V1 <- NULL
prom500.assocTSS.hypo.f$V2 <- NULL
colnames(prom500.assocTSS.hypo.f) <- c("dist.to.feature", "feature.name", "feature.strand")


# 3. Aggregate the data according to the feature.name, collapsing the columns with dist.to.feature and target.row ids

dist.hyper <- aggregate(dist.to.feature ~ feature.name, data = prom500.assocTSS.hyper.f, paste, collapse = ",")

feature.strand.hyper <- aggregate(feature.strand ~ feature.name, data = prom500.assocTSS.hyper.f, paste, collapse = ",")


dist.hypo <- aggregate(dist.to.feature ~ feature.name, data = prom500.assocTSS.hypo.f, paste, collapse = ",")

feature.strand.hypo <- aggregate(feature.strand ~ feature.name, data = prom500.assocTSS.hypo.f, paste, collapse = ",")

# 5. Cbind both datasets and delete the column with the feature.names, that is duplicated

prom500.id.hyper <- cbind(dist.hyper, feature.strand.hyper)

prom500.id.hyper[,3] <- NULL

prom500.id.hypo <- cbind(dist.hypo, feature.strand.hypo)

prom500.id.hypo[,3] <- NULL

# 6. Save the results and the workspace

write.table(prom500.id.hyper, "prom500_assocTSS_hyper_f_aggreg.txt", quote = F, sep = "\t", row.names = F)

write.table(prom500.id.hypo, "prom500_assocTSS_hypo_f_aggreg.txt", quote = F, sep = "\t", row.names = F)

save.image("f.RData")
