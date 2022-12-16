library (edgeR)
library (limma)
library (RColorBrewer)
library (mixOmics)
library (HTSFilter)

# Referenced (http://www.nathalievialaneix.eu/doc/html/solution_edgeR-tomato-withcode.html)

# Import data (below example)
# GeneID	ST5	ST6	ST7	STx
# Wdr13		370	431	453	504
# Npnt		599	582	655	731
# Twnk		901	969	1000	1165
# Gene X	.	.	.	.

rawCountTable <- read.table ("S0-5678-exp-ID.txt", header=TRUE, sep="\t", row.names = 1)


# (Optional) If you want to subset your samples 
rawCountTable1234 <- rawCountTable[-c(1,2,3,4)]


# Import study design (below example)
# File  condition
# ST5   Mut
# ST6   WT
# ST7   Mut
# STx   .

sampleInfo <- read.table("design2.csv", header = TRUE, sep=",", row.names = 1)


# Create a DGEList data object
dgeFull <- DGEList (rawCountTable, group=sampleInfo$condition)

# == Data exploration ==

# Create pesudo-counts
pseudoCounts <- log2(dgeFull$counts+1)

# histogram & box plots
hist(pseudoCounts[,"ST5"])
boxplot (pseudoCounts, col="gray", las=3)

# MA plot
par(mfrow = c(1,2))
avalues <- (pseudoCounts[,1] + pseudoCounts[,2])/2
mvalues <- (pseudoCounts[,1] - pseudoCounts[,2])
plot (avalues, mvalues, xlab = "A", ylab ="M", pch = 19, main = "Sample1 vs Sample2")

# MDS plot
plotMDS (pseudoCounts)

# HeatMap
sampleDists <- as.matrix(dist(t(pseudoCounts)))
sampleDists
cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(16)
cim (sampleDists, color = cimColor, symkey=FALSE)


# == DEG test ==

# Remove genes with zero counts for all samples
dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum)!=0, ], group=dgeFull$samples$group)
# (e.g. 24984 elements to 21404 elements)

head(dgeFull$counts)

# Estimate normalization factor 
# (Trimmed Mean of M-values (TMM): 
# assuming the majority of 'housekeeping' genes have the same expression levels across conditions)
dgeFull <- calcNormFactors(dgeFull, method="TMM")
eff.lib.size <- dgeFull$samples$lib.size*dgeFull$samples$norm.factors

# I think, 
# normCounts = {count/(lib.size/1000000)}/norm.factor
normCounts <- cpm (dgeFull)
pseudoNormCounts <- log2 (normCounts +1)
boxplot (pseudoNormCounts, col="gray", las =3)

plotMDS (pseudoNormCounts)

# Estimate common and tagwise dispersion
dgeFull <- estimateCommonDisp(dgeFull)
dgeFull <- estimateTagwiseDisp(dgeFull)


# Perform exact test for differential expression
dgeTest <- exactTest(dgeFull)


resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table))

write.csv (resNoFilt, file="testResult.csv")
write.csv (pseudoNormCounts, file="NormCPM.csv")
