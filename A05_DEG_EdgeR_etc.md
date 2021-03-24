# Followed guide below

Practical statistical analysis of RNA-Seq data - edgeR - tomato data
Annick Moisan, Ignacio Gonzales, Nathalie Villa-Vialaneix
March, 10th 2015; updated September 6th

http://www.nathalievialaneix.eu/doc/html/solution_edgeR-tomato-withcode.html


# Copy and paste of R commands used

<pre>
> rawCountTable <- read.table("CountTotal.txt", header=TRUE, sep="\t", row.names=1)
> head(rawCountTable)
> sampleInfo <- read.table("design.csv", header=TRUE, sep=",", row.names=1)
> dgeFull <- DGEList()
> dgeFull <- DGEList(rawCountTable, group = sampleInfo$Condition)
> dgeFull
> pseudoCounts <- log2(dgeFull$counts+1)
> head(pseudoCounts)
> hist(pseudoCounts[,"S01"])
> boxplot (pseudoCounts, col="gray", las=3)
> plotMDS(pseudoCounts)
> dgeFull2 <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum)!= 0,], group=dgeFull$samples$group)
> head (dgeFull2$counts)
> View(dgeFull)
> View(dgeFull2)
> dgeFull2 <- calcNormFactors(dgeFull2, method ="TMM")
> eff.lib.size <- dgeFull2$samples$lib.size*dgeFull$samples$norm.factors
> normConts <- cpm (dgeFull2)
> pseudoNormCounts <- log2(normConts+1)
> boxplot(pseudoNormCounts, col='gray', las=3)
> dgeFull2 <- estimateCommonDisp(dgeFull2)
> dgeFull2 <- estimateTagwiseDisp(dgeFull2)
> dgeTest <- exactTest(dgeFull2, pair = c("WT", "Homo"))
> hist(dgeTest$Table[,"PValue"], breaks = 50)
> hist(dgeTest$table[,"PValue"], breaks = 50)
> resNoFilt <- topTags (dgeTest, n=nrow(dgeTest$table))
> sigDownReg <-resNoFilt$table[resNoFilt$table$PValue<0.05,]
> sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
> sigUpReg <- sigDownReg[order(sigDownReg$logFC, decreasing=TRUE),]
> write.csv(sigDownReg, file="sigDownReg.csv")
> write.csv(sigUpReg, file="sigUpReg.csv")
</pre>

# Additional way

1. Raw transcripts counts 
2. Quantile normalized in matlab 
   - NormData = quantilenorm(Data) 
3. Explore data using Matlab and Excel. 
