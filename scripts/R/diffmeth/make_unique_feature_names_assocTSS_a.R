# 1. Load the saved workspace with the R variables

# load("a.RData")

# 2. Convert row names as a column

#library(data.table)
#setDT(assocTSS.hyper.a, keep.rownames = T)

#setDT(assocTSS.hypo.a, keep.rownames = T)

# This will return the follow:

# > head(assocTSS.hyper.a)
# rn target.row dist.to.feature feature.name feature.strand
# 1:   108          1             817 NM_001011874              -
# 2: 108.1          2             814 NM_001011874              -
# 3: 108.2          3             375 NM_001011874              -
# 4: 108.3          4             337 NM_001011874              -
# 5: 108.4          5             335 NM_001011874              -
# 6: 108.5          6              12 NM_001011874              -

# 3. Delete the row named "target.row" and rename the row "rn" as "target.row"

# assocTSS.hyper.a$target.row <- NULL
# colnames(assocTSS.hyper.a) <- c( "target.row", "dist.to.feature", "feature.name", "feature.strand")

# assocTSS.hypo.a$target.row <- NULL
# colnames(assocTSS.hypo.a) <- c( "target.row", "dist.to.feature", "feature.name", "feature.strand")

# 4. Aggregate the data according to the feature.name, collapsing the columns with dist.to.feature and target.row ids

dist.hyper <- aggregate(dist.to.feature ~ feature.name, data = assocTSS.hyper.a, paste, collapse = ",")

target.row.hyper <- aggregate(target.row ~ feature.name, data = assocTSS.hyper.a, paste, collapse = ",")


dist.hypo <- aggregate(dist.to.feature ~ feature.name, data = assocTSS.hypo.a, paste, collapse = ",")

target.row.hypo <- aggregate(target.row ~ feature.name, data = assocTSS.hypo.a, paste, collapse = ",")

# 5. Cbind both datasets and delete the column with the feature.names, that is duplicated

id.hyper <- cbind(dist.hyper, target.row.hyper)

id.hyper[,3] <- NULL

id.hypo <- cbind(dist.hypo, target.row.hypo)

id.hypo[,3] <- NULL

# 6. Save the results and the workspace

write.table(id.hyper, "assocTSS_hyper_a_aggreg.txt", quote = F, sep = "\t")

write.table(id.hypo, "assocTSS_hypo_a_aggreg.txt", quote = F,  sep = "\t")

save.image("a.RData")



