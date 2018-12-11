################################################################################
# Load libraries 
################################################################################
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tweeDEseq))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(DESeq))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(GOstats))
suppressPackageStartupMessages(library(org.At.tair.db))

################################################################################
# Read data 
################################################################################
# Get samples from folder names
sample.dir <- "/mnt/518D6BCF3ECC578E/RNAseq/results/STAR_mapping"
sample.names <- list.dirs(sample.dir, recursive = F, full.names = F)
# Read first sample
###################
cat("Reading sample", sample.names[1], "\n")
counts.matrix <- read.table(paste0(sample.dir, "/", sample.names[1], "/", "gene_assigned_P"), header = T)
colnames(counts.matrix)[7] <- sample.names[1]
# Read other samples
####################
for (sample in sample.names[-1]){
  cat("Reading sample", sample, "\n")
  temp.matrix <- read.table(paste0(sample.dir, "/", sample, "/", "gene_assigned_P"), header = T)
  colnames(temp.matrix)[7] <-sample
  counts.matrix <- merge(counts.matrix, temp.matrix, by = colnames(counts.matrix)[1:6], all = T)
}
rm(temp.matrix)

################################################################################
# Add sample data 
################################################################################
coldata <- data.frame("condition" = c("Ws_90", "Ws_90", "Ws_90", "mkp1_90", 
                                      "mkp1_90", "mkp1_90"), 
                      row.names = sample.names)
Ws_90.samples <- which(coldata$condition == "Ws_90")
mkp1_90.samples <- which(coldata$condition == "mkp1_90")
cts <- counts.matrix[-c(1:6)]
row.names(cts) <- counts.matrix$Geneid
# check sample names match
all(rownames(coldata) == colnames(cts))

################################################################################
# Normalization 
################################################################################
# TMM
cts.TMM <- normalizeCounts(cts, method="TMM")

# FPKM
width <- counts.matrix$Length
cts.FPKM <- t(t(cts / width * 1000)/colSums(cts) * 1e6)

# Check normalization
#####################
# Raw counts
maPlot(rowMeans(cts[,mkp1_90.samples]),
       rowMeans(cts[,Ws_90.samples]),
       pch=19, cex=.5, ylim=c(-8,8), allCol="darkgray", lowess=TRUE,
       main="Raw counts",
       xlab=expression(A == log[2](sqrt(mkp1_90 %.% Ws_90))),
       ylab=expression(M == log[2](mkp1_90)-log[2](Ws_90)))
grid(col="black")

# TMM
maPlot(rowMeans(cts.TMM[,mkp1_90.samples]),
       rowMeans(cts.TMM[,Ws_90.samples]),
       pch=19, cex=.5, ylim=c(-8,8), allCol="darkgray", lowess=TRUE,
       main="TMM normalization",
       xlab=expression(A == log[2](sqrt(mkp1_90 %.% Ws_90))),
       ylab=expression(M == log[2](mkp1_90)-log[2](Ws_90)))
grid(col="black")

# FPKM
maPlot(rowMeans(cts.FPKM[,mkp1_90.samples]),
       rowMeans(cts.FPKM[,Ws_90.samples]),
       pch=19, cex=.5, ylim=c(-8,8), allCol="darkgray", lowess=TRUE,
       main="FPKM normalization",
       xlab=expression(A == log[2](sqrt(mkp1_90 %.% Ws_90))),
       ylab=expression(M == log[2](mkp1_90)-log[2](Ws_90)))
grid(col="black")


################################################################################
# Differential gene expression analysis 
################################################################################
# create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = cts.TMM,
                              colData = coldata,
                              design= ~ condition)

# set control sample
dds$condition <- relevel(dds$condition, ref = "Ws_90")

# Analysis
dds <- DESeq(dds)
res <- results(dds)
res

# Summary
summary(res)

# Select significantly regulated genes (FDR < 0.05)
res.sig <- res[ res$padj < 0.05 & !is.na(res$padj), ]
summary(res.sig)

################################################################################
# Heatmap
################################################################################
# Counts of the identified DEG
DEG.names <- as.character(rownames(res.sig))
cts.TMM.DEG <- cts.TMM[DEG.names,]

# Generate the plot
heatmap.2(as.matrix(cts.TMM.DEG), dendrogram = "column", scale = "row", 
          labRow = FALSE, trace = "none")

################################################################################
# Ring Plot 
################################################################################
# Define UP and DOWN genes
UP.DEG <- DEG.names[res.sig$log2FoldChange > 0]
DOWN.DEG <- DEG.names[res.sig$log2FoldChange < 0]

# Prepare data for the plot
plot.df <- data.frame("genes"=c(length(UP.DEG), length(DOWN.DEG)), 
                      "condition"=c("UP","DOWN"))
plot.df$fraction <- plot.df$genes / sum(plot.df$genes)
plot.df <- plot.df[order(plot.df$fraction), ]
plot.df$ymax <- cumsum(plot.df$fraction)
plot.df$ymin = c(0, head(plot.df$ymax, n=-1))
non.DEG <- length(counts.matrix$Geneid) - length(DEG.names)

# PLot
ggplot(plot.df, aes(fill=condition, ymax=ymax, ymin=ymin, xmax=4, xmin=2.5)) +
  geom_rect(colour="grey30") +
  coord_polar(theta="y") +
  scale_fill_manual("",values=c('DOWN'='green','UP'='red')) +
  theme_bw() +
  theme(panel.grid=element_blank(), axis.text=element_blank()) +
  theme(axis.ticks=element_blank()) +
  labs(title="Gene Expression Ring Plot", x="", y="") +
  annotate("text", x = 0, y = 0, label = paste("non-DEGs:", non.DEG)) +
  geom_text(aes(label=plot.df$genes,x=4.75,y=(ymin+ymax)/2)) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("DEG.ring.plot.png")

################################################################################
# Gene Onthology Enrichment Analysis 
################################################################################
# Define universe of genes
geneUniverse <- as.character(counts.matrix$Geneid)
length(geneUniverse)

# Hypergeometrical test of Biological Process (BP)
##################################################
# set analysis parameters
params.GO <- new("GOHyperGParams", geneIds=DEG.names,
                 universeGeneIds=geneUniverse,
                 annotation="org.At.tair.db", ontology="BP",
                 pvalueCutoff=0.05, conditional=FALSE,
                 testDirection="over")
# Functinal enrichment
DEG.GO.BP <- hyperGTest(params.GO)
DEG.GO.BP

# conditional test
conditional(params.GO) <- TRUE
DEG.GO.BP.Cond <- hyperGTest(params.GO)
DEG.GO.BP.Cond

# top ten most significant GO terms
head(summary(DEG.GO.BP.Cond))

# Save results as html report
#############################
htmlReport(DEG.GO.BP.Cond, file="DEG_GO_BP_Cond.html")

# REVIGO summary and visualization of GOs
#########################################
# Prepare input data
write.table(summary(DEG.GO.BP.Cond)[,c(1,2)], "REVIGO.input.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
# Load into http://revigo.irb.hr/ website


