setwd("~/phd/celine/biotypes/")
library(DESeq2)

mirna <- read.csv("mirbase_subread_matrix.txt", sep = "\t", row.names = 1)
pirna <- read.csv("piRNA_subread_matrix.txt", sep = "\t",row.names = 1)

meta <- read.table("../celine_samples.txt", row.names = 1, header = T)
meta$Age <- factor(meta$Age, levels = c("12","20"))
meta$Design.Age <- factor(meta$Design.Age, 
                          levels = c("WT-SH.12","HD-SH.12","WT-EE.12","HD-EE.12",
                                     "WT-SH.20","HD-SH.20","WT-EE.20","HD-EE.20"))

#----------------------------------------------------------------------------------------------------
# Investigating EE effect in WT mice at 12 weeks of age
colnames(mirna) <- rownames(meta)
index <- which(meta$Age == "12" & meta$Genotype == "WT")
mirna.12.wt <- mirna[,index]
meta.12.wt <- meta[index,]
dds.12.wt <- DESeqDataSetFromMatrix(countData = mirna.12.wt,
                                 colData = meta.12.wt,
                                 design= ~ Housing)

keep <- rowSums(counts(dds.12.wt)) >= 250
dds.12.wt <- dds.12.wt[keep,]
dim(dds.12.wt)
dds.12.wt$Design <- relevel(dds.12.wt$Design, ref = "WT-SH")

dds.12.wt <- DESeq(dds.12.wt)
res.12.WT.SHvsEE <- results(dds.12.wt, contrast = c("Housing","SH", "EE"))
res.12.WT.SHvsEE <- results(dds.12, contrast = c("Design","HD-EE", "WT-EE"))

res.12.WT.SHvsEE <- res.12.WT.SHvsEE[order(res.12.WT.SHvsEE$pvalue),]
sum(res.12.WT.SHvsEE$padj < 0.05, na.rm=TRUE)

# Investigating EE effect in WT mice at 20 weeks of age
index <- which(meta$Age == "20" & meta$Genotype == "WT")
mirna.20.wt <- mirna[,index]
meta.20.wt <- meta[index,]
dds.20.wt <- DESeqDataSetFromMatrix(countData = mirna.20.wt,
                                    colData = meta.20.wt,
                                    design= ~ Housing)

keep <- rowSums(counts(dds.20.wt)) >= 250
dds.20.wt <- dds.20.wt[keep,]
dim(dds.20.wt)
dds.20.wt$Design <- relevel(dds.20.wt$Design, ref = "WT-SH")

dds.20.wt <- DESeq(dds.20.wt)
res.20.WT.SHvsEE <- results(dds.20.wt, contrast = c("Housing","SH", "EE"))
res.12.WT.SHvsEE <- results(dds.12, contrast = c("Design","HD-EE", "WT-EE"))

res.20.WT.SHvsEE <- res.20.WT.SHvsEE[order(res.20.WT.SHvsEE$pvalue),]
sum(res.20.WT.SHvsEE$padj < 0.05, na.rm=TRUE)
resSig.20.WT.SHvsEE <- subset(res.20.WT.SHvsEE, padj < 0.05)


# Investigate early vs late WT age effect on miRNA
index <- which(meta$Design == "WT-SH")
mirna.wt.sh <- mirna[,index]
meta.wt.sh <- meta[index,]

dds.wt.sh <- DESeqDataSetFromMatrix(countData = mirna.wt.sh,
                              colData = meta.wt.sh,
                              design= ~ Age)
keep <- rowSums(counts(dds.wt.sh)) >= 250
dds.wt.sh <- dds.wt.sh[keep,]
dim(dds.wt.sh)
dds.wt.sh$Design <- relevel(dds.wt.sh$Design.Age, ref = "12")

dds.wt.sh <- DESeq(dds.wt.sh)
res.12vs20.WT.SH <- results(dds.wt.sh, contrast = c("Age","20", "12"))

res.12vs20.WT.SH <- res.12vs20.WT.SH[order(res.12vs20.WT.SH$pvalue),]
sum(res.12vs20.WT.SH$padj < 0.05, na.rm=TRUE)
sum(res.12vs20.WT.SH$padj < 0.05)
sum(res.12vs20.WT.SH$log2FoldChange < -0.7 & res.12vs20.WT.SH$padj < 0.05)

# Investigate early vs late WT EE effect on miRNA
index <- which(meta$Design == "WT-EE")
mirna.wt.ee <- mirna[,index]
meta.wt.ee <- meta[index,]

dds.wt.ee <- DESeqDataSetFromMatrix(countData = mirna.wt.ee,
                                    colData = meta.wt.ee,
                                    design= ~ Age)
keep <- rowSums(counts(dds.wt.ee)) >= 250
dds.wt.ee <- dds.wt.ee[keep,]
dim(dds.wt.ee)
dds.wt.ee$Design <- relevel(dds.wt.ee$Design.Age, ref = "12")

dds.wt.ee <- DESeq(dds.wt.ee)
res.12vs20.WT.EE <- results(dds.wt.ee, contrast = c("Age","20", "12"))

res.12vs20.WT.EE <- res.12vs20.WT.EE[order(res.12vs20.WT.EE$pvalue),]
sum(res.12vs20.WT.EE$padj < 0.05, na.rm=TRUE)
sum(res.12vs20.WT.SH$log2FoldChange < -0.7 & res.12vs20.WT.SH$padj < 0.05)

#---------------------------------------------------------------------------------
# Investigating early vs late HD effect on miRNA
colnames(mirna) <- rownames(meta)
index <- which(meta$Design == "HD-SH")
mirna.hd.sh <- mirna[,index]
meta.hd.sh <- meta[index,]

dds <- DESeqDataSetFromMatrix(countData = mirna.hd.sh,
                              colData = meta.hd.sh,
                              design= ~ Age)
keep <- rowSums(counts(dds)) >= 250
dds <- dds[keep,]
dim(dds)
dds$Design <- relevel(dds$Design.Age, ref = "12")

dds <- DESeq(dds)
res.12vs20.HD.SH <- results(dds, contrast = c("Age","20", "12"))

res.12vs20.HD.SH <- res.12vs20.HD.SH[order(res.12vs20.HD.SH$pvalue),]
sum(res.12vs20.HD.SH$padj < 0.05, na.rm=TRUE)
sum(res.12vs20.HD.SH$padj < 0.05)
sum(res.12vs20.HD.SH$log2FoldChange < -0.7 & res.12vs20.HD.SH$padj < 0.05)

resultsNames(dds)
plotMA(res.12vs20.HD.SH, ylim=c(-2,2))

plotCounts(dds, gene=which.min(res.12vs20.HD.SH$padj), intgroup="Age")
res.12vs20.HD.SH <- subset(res.12vs20.HD.SH, padj < 0.05 & log2FoldChange > 1)
resSig.12vs20 <- subset(res.12vs20.HD.SH, padj < 0.05 & log2FoldChange < -1)
resSig.12vs20 <- subset(res.12vs20.HD.SH, padj < 0.05)
res.12vs20.HD.SH <- rbind(res.12vs20.HD.SH, resSig.12vs20)

write.table(as.data.frame(resSig.12), file="mirna_HD_SH_12vs20_padj0.05.txt", sep = "\t", quote=FALSE)

which(rownames(resSig.12) %in% rownames(resSig.20))
mir_names <- c("MIMAT0000612", "	MIMAT0017044")
res.12vs20.HD.SH["MIMAT0000127",]

#----------------------------------------------------------------------------------------------------
# Investigating HD effect at 12 weeks of age
index <- which(meta$Age == "12")
mirna.12 <- mirna[,index]
meta.12 <- meta[index,]
dds.12 <- DESeqDataSetFromMatrix(countData = mirna.12,
                              colData = meta.12,
                              design= ~ Design)

keep <- rowSums(counts(dds.12)) >= 250
dds.12 <- dds.12[keep,]
dim(dds.12)
dds.12 <- DESeq(dds.12)
res.12.WTvsHD.SH <- results(dds.12, contrast = c("Design","HD-SH", "WT-SH"))
sum(res.12.WTvsHD.SH$padj < 0.05, na.rm=TRUE)
res.12.WTvsHD.SH <- res.12.WTvsHD.SH[order(res.12.WTvsHD.SH$pvalue),]
sum(res.12.WTvsHD.SH$log2FoldChange < -0.5 & res.12.WTvsHD.SH$padj < 0.05)
sum(res.12.WTvsHD.SH$log2FoldChange > 0.5 & res.12.WTvsHD.SH$padj < 0.05)

res.12.SHvsEE.HD <- results(dds.12, contrast = c("Design","HD-EE", "HD-SH"))
sum(res.12.SHvsEE.HD$padj < 0.05, na.rm=TRUE)

resultsNames(dds.12)
plotMA(res.12.WTvsHD.SH, ylim=c(-2,2))

plotCounts(dds.12, gene=which.min(res.12.WTvsHD.SH$padj), intgroup="Design")
resSig.12.WTvsHD.SH <- subset(res.12.WTvsHD.SH, padj < 0.05 & log2FoldChange > 1)
resSig.12 <- subset(res.12.WTvsHD.SH, padj < 0.05 & log2FoldChange < -1)
resSig.12 <- subset(res.12.WTvsHD.SH, padj < 0.05)
resSig.12.WTvsHD.SH <- rbind(resSig.12.WTvsHD.SH, resSig.12)

resSig.12.WTvsHD.SH["MIMAT0000612",]

write.table(as.data.frame(resSig.12), file="mirna_HDvsWT_SH_12_padj0.05.txt", sep = "\t", quote=FALSE)

#----------------------------------------------------------------------------------------------------
# Investigating HD effect at 20 weeks of age
index <- which(meta$Age == "20")
mirna.20 <- mirna[,index]
meta.20 <- meta[index,]
dds.20 <- DESeqDataSetFromMatrix(countData = mirna.20,
                                 colData = meta.20,
                                 design= ~ Design)

keep <- rowSums(counts(dds.20)) >= 250
dds.20 <- dds.20[keep,]
dim(dds.20)
dds.20$Design <- relevel(dds.20$Design, ref = "WT-SH")

dds.20 <- DESeq(dds.20)
res.20.WTvsHD.SH <- results(dds.20, contrast = c("Design","HD-SH", "WT-SH"))
sum(res.20.WTvsHD.SH$padj < 0.05, na.rm=TRUE)
res.20.WTvsHD.SH <- res.20.WTvsHD.SH[order(res.20.WTvsHD.SH$pvalue),]
sum(res.20.WTvsHD.SH$log2FoldChange > 0.5 & res.20.WTvsHD.SH$padj < 0.05)
sum(res.20.WTvsHD.SH$log2FoldChange < -0.5 & res.20.WTvsHD.SH$padj < 0.05)

res.20.SHvsEE.HD <- results(dds.20, contrast = c("Design","HD-EE", "HD-SH"))
sum(res.20.SHvsEE.HD$padj < 0.05, na.rm=TRUE)

res.20.SHvsEE.HD["MIMAT0016984",] #mmu-miR-132-5p
res.20.SHvsEE.HD["MIMAT0004522",] #mmu-miR-27b-5p
res.20.SHvsEE.HD["MIMAT0025138",] #mmu-miR-378c

resultsNames(dds.20)
plotMA(res.20.WTvsHD.SH, ylim=c(-2,2))

plotCounts(dds.20, gene=which.min(res.20.WTvsHD.SH$padj), intgroup="Design")
resSig.20.WTvsHD.SH <- subset(res.20.WTvsHD.SH, padj < 0.05 & log2FoldChange > 1)
resSig.20 <- subset(res.20.WTvsHD.SH, padj < 0.05 & log2FoldChange < -1)
resSig.20 <- subset(res.20.WTvsHD.SH, padj < 0.05)
resSig.20.WTvsHD.SH <- rbind(resSig.20.WTvsHD.SH, resSig.20)

write.table(as.data.frame(resSig.20), file="mirna_HDvsWT_SH_20_padj0.05.txt", sep = "\t", quote=FALSE)

#----------------------------------------------------------------------------------------------------
# Investigating age effect in piRNAs
colnames(pirna) <- rownames(meta)

index <- which(meta$Genotype == "WT")
pirna.subset <- pirna[,index]
meta.subset <- meta[index,]
dds.pi <- DESeqDataSetFromMatrix(countData = pirna.subset,
                                 colData = meta.subset,
                                 design= ~ Age)

keep <- rowSums(counts(dds.pi)) >= 250
dds.pi<- dds.pi[keep,]
dim(dds.pi)

dds.pi <- DESeq(dds.pi)
res.pi <- results(dds.pi, contrast = c("Age","20", "12"))
sum(res.pi$padj < 0.05, na.rm=TRUE)

res.pi <- res.pi[order(res.pi$pvalue),]
res.pi <- subset(res.pi, padj < 0.05)

res.pi.WT.12vs20 <- res.pi
write.table(res.pi.WT.12vs20, file = "piRNA_DE_WT_12vs20__Deseq2.txt", sep = "\t", quote = F)

res.pi.WT.SH.12vs20 <- res.pi
res.pi.WT.EE.12vs20 <- res.pi
res.pi.HD.SH.12vs20 <- res.pi
res.pi.HD.EE.12vs20 <- res.pi

# Investigating HD effect in pirnas at 12 weeks of age
index <- which(meta$Age == "12")
pirna.12 <- pirna[,index]
meta.12 <- meta[index,]
dds.12.pi <- DESeqDataSetFromMatrix(countData = pirna.12,
                                    colData = meta.12,
                                    design= ~ Design)

keep <- rowSums(counts(dds.12.pi)) >= 250
dds.12.pi <- dds.12.pi[keep,]
dim(dds.12.pi)

dds.12.pi <- DESeq(dds.12.pi)
res.12.WTvsHD.SH.pi <- results(dds.12.pi, contrast = c("Design","HD-SH", "WT-SH"))
sum(res.12.WTvsHD.SH.pi$padj < 0.05, na.rm=TRUE)
res.12.WTvsHD.SH.pi <- res.12.WTvsHD.SH.pi[order(res.12.WTvsHD.SH.pi$pvalue),]
subset(res.12.WTvsHD.SH.pi, padj < 0.05)

res.12.SHvsEE.HD.pi <- results(dds.12.pi, contrast = c("Design","HD-EE", "HD-SH"))
sum(res.12.SHvsEE.HD.pi$padj < 0.05, na.rm=TRUE)

res.12.SHvsEE.WT.pi <- results(dds.12.pi, contrast = c("Design","WT-EE", "WT-SH"))
sum(res.12.SHvsEE.WT.pi$padj < 0.05, na.rm=TRUE)

sum(res.12.WTvsHD.SH.pi$log2FoldChange > 0.5 & res.12.WTvsHD.SH.pi$padj < 0.05)
sum(res.12.WTvsHD.SH.pi$log2FoldChange < -0.5 & res.12.WTvsHD.SH.pi$padj < 0.05)

# Investigating HD effect in pirnas at 20 weeks of age
index <- which(meta$Age == "20")
pirna.20 <- pirna[,index]
meta.20 <- meta[index,]
dds.20.pi <- DESeqDataSetFromMatrix(countData = pirna.20,
                                 colData = meta.20,
                                 design= ~ Design)

keep <- rowSums(counts(dds.20.pi)) >= 250
dds.20.pi <- dds.20.pi[keep,]
dim(dds.20.pi)
dds.20.pi <- DESeq(dds.20.pi)

res.20.WTvsHD.SH.pi <- results(dds.20.pi, contrast = c("Design","HD-SH", "WT-SH"))
sum(res.20.WTvsHD.SH.pi$padj < 0.05, na.rm=TRUE)
sum(res.20.WTvsHD.SH.pi$log2FoldChange > 0.5 && res.20.WTvsHD.SH.pi$padj < 0.05)
subset(res.20.WTvsHD.SH.pi, padj < 0.05)

res.20.WTvsHD.SH.pi <- res.20.WTvsHD.SH.pi[order(res.20.WTvsHD.SH.pi$pvalue),]
resSig.20.pi.WTvsHD.SH <- subset(res.20.WTvsHD.SH.pi, padj < 0.05)
write.table(resSig.20.pi.WTvsHD.SH, file = "piRNA_DE_20wks_HDvsWT_SH_Deseq2.txt", sep = "\t", quote = F)

res.20.SHvsEE.HD.pi <- results(dds.20.pi, contrast = c("Design","HD-EE", "HD-SH"))
sum(res.20.SHvsEE.HD.pi$padj < 0.05, na.rm=TRUE)

res.20.SHvsEE.WT.pi <- results(dds.20.pi, contrast = c("Design","WT-EE", "WT-SH"))
sum(res.20.SHvsEE.WT.pi$padj < 0.05, na.rm=TRUE)


resultsNames(dds.20)
plotMA(res.20.WTvsHD.SH, ylim=c(-2,2))

plotCounts(dds.20, gene=which.min(res.20.WTvsHD.SH$padj), intgroup="Design")
resSig.20.WTvsHD.SH <- subset(res.20.WTvsHD.SH, padj < 0.05 & log2FoldChange > 1)
resSig.20 <- subset(res.20.WTvsHD.SH, padj < 0.05 & log2FoldChange < -1)
resSig.20 <- subset(res.20.WTvsHD.SH, padj < 0.05)
resSig.20.WTvsHD.SH <- rbind(resSig.20.WTvsHD.SH, resSig.20)

# Determining composition difference between WT_SH and WT-EE at 12 and 20 wks
colnames(pirna) = rownames(meta)
index <- which(meta$Age == "20" & meta$Genotype == "WT")
pirna.20.wt <- pirna[,index]
meta.20.wt <- meta[index,]


library(phyloseq)
phy.obj = phyloseq(otu_table(pirna.20.wt, taxa_are_rows = TRUE), sample_data(meta.20.wt))

## Rarefy for alpha-div
set.seed(223)
phy.obj.rare = rarefy_even_depth(phy.obj, sample.size = min(sample_sums(phy.obj)))

phy.obj.ra = transform_sample_counts(phy.obj, function(x) x/sum(x))
### Bray Curtis distance
distBC.20 = distance(phy.obj.ra, method = "bray")
ordBC.20 = ordinate(phy.obj.ra, method = "PCoA", distance = distBC)
vegan::adonis(distBC.12 ~ Housing, data = meta.12.wt)
vegan::adonis(distBC.20 ~ Housing, data = meta.20.wt)

jpeg('Bray_curtis.jpg', units = "in", width = 9, height = 4.5, res = 300)
plot_ordination(phy.obj.ra, ordBC.20, color = "Housing", shape = "Housing") +
  ggtitle("Bray Curtis distance") +
  geom_jitter(aes(color = Housing)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(name = "Component 1") +
  scale_y_continuous(name = "Component 2") 
#scale_color_manual(values = c('#388ECC','#F68B33'))#+ geom_point(size = 5)
plot_ordination(phy.obj.ra, ordBC, color = "Status", shape = "Gene") +
  ggtitle("Bray Curtis distance") +
  geom_jitter(aes(color = Status)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_grid(~Gender) +
  scale_x_continuous(name = "Component 1") +
  scale_y_continuous(name = "Component 2") 
dev.off()

# Determining composition difference between HD_SH and HD_EE at 12 and 20 wks
index <- which(meta$Age == "12" & meta$Genotype == "HD")
pirna.12.hd <- pirna[,index]
meta.12.hd <- meta[index,]


library(phyloseq)
phy.obj = phyloseq(otu_table(pirna.12.hd, taxa_are_rows = TRUE), sample_data(meta.12.hd))

## Rarefy for alpha-div
set.seed(223)
phy.obj.rare = rarefy_even_depth(phy.obj, sample.size = min(sample_sums(phy.obj)))

phy.obj.ra = transform_sample_counts(phy.obj, function(x) x/sum(x))
### Bray Curtis distance
distBC.12.hd = distance(phy.obj.ra, method = "bray")
ordBC.12.hd = ordinate(phy.obj.ra, method = "PCoA", distance = distBC)
vegan::adonis(distBC.12.hd ~ Housing, data = meta.12.hd)
vegan::adonis(distBC.20.hd ~ Housing, data = meta.20.hd)


