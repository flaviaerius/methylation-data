# 1. Load the saved workspace with the R variables

# load("c.RData")

# 2. Convert row names as a column

# library(data.table)
# setDT(assocTSS.hyper.c, keep.rownames = T)
# 
# setDT(assocTSS.hypo.c, keep.rownames = T)

# This will return the follow:

# > head(assocTSS.hyper.c)
# rn target.row dist.to.feature feature.name feature.strand
# 1:   108          1             817 NM_001011874              -
# 2: 108.1          2             814 NM_001011874              -
# 3: 108.2          3             375 NM_001011874              -
# 4: 108.3          4             337 NM_001011874              -
# 5: 108.4          5             335 NM_001011874              -
# 6: 108.5          6              12 NM_001011874              -

# 3. Delete the row named "target.row" and rename the row "rn" as "target.row"

# assocTSS.hyper.c$target.row <- NULL
# colnames(assocTSS.hyper.c) <- c( "target.row", "dist.to.feature", "feature.name", "feature.strand")
# 
# assocTSS.hypo.c$target.row <- NULL
# colnames(assocTSS.hypo.c) <- c( "target.row", "dist.to.feature", "feature.name", "feature.strand")

# 4. Aggregate the data according to the feature.name, collapsing the columns with dist.to.feature and target.row ids

dist.hyper <- aggregate(dist.to.feature ~ feature.name, data = assocTSS.hyper.c, paste, collapse = ",")

target.row.hyper <- aggregate(target.row ~ feature.name, data = assocTSS.hyper.c, paste, collapse = ",")


dist.hypo <- aggregate(dist.to.feature ~ feature.name, data = assocTSS.hypo.c, paste, collapse = ",")

target.row.hypo <- aggregate(target.row ~ feature.name, data = assocTSS.hypo.c, paste, collapse = ",")

# 5. Cbind both datasets and delete the column with the feature.names, that is duplicated

id.hyper <- cbind(dist.hyper, target.row.hyper)

id.hyper[,3] <- NULL

id.hypo <- cbind(dist.hypo, target.row.hypo)

id.hypo[,3] <- NULL

# 6. Save the results and the workspace

write.table(id.hyper, "assocTSS_hyper_c_aggreg.txt", quote = F, sep = "\t")

write.table(id.hypo, "assocTSS_hypo_c_aggreg.txt", quote = F, sep = "\t")

save.image("c.RData")



