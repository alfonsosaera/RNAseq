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
# Volcano plot
################################################################################
#
# Prepare data
#
volcano.plot.df <- data.frame(gene=rownames(res)[!is.na(res$pvalue)],
                              log2FC=res$log2FoldChange[!is.na(res$pvalue)],
                              logpval=-log10(res$pvalue[!is.na(res$pvalue)]),
                              fdr=res$padj[!is.na(res$pvalue)],
                              stringsAsFactors = FALSE)

# replace infinite -log10(pvalue) with the maximum value that is not infinite + 15%
max.pvalue <- max(volcano.plot.df$logpval[is.finite(volcano.plot.df$logpval)]) * 1.15
infinite.pvalue <- which(is.infinite(volcano.plot.df$logpval))
volcano.plot.df$logpval[infinite.pvalue] <- max.pvalue

# replace NA values in fdr with 1
volcano.plot.df$fdr[which(is.na(volcano.plot.df$fdr))] <- 1

# reorder dataframe
volcano.plot.df <- volcano.plot.df[order(-volcano.plot.df$fdr),]

#
# Generate volcano plot
#
# Use color to show UP and Down regulation
##########################################
# Define color
volcano.plot.df$color <- ifelse(volcano.plot.df$fdr < 0.0001 & volcano.plot.df$log2FC < -1, "Down-regulated",
                                ifelse(volcano.plot.df$fdr < 0.0001 & volcano.plot.df$log2FC > 1, "UP-regulated",
                                       "Non-DEG"))
volcano.plot.df$color <- factor(volcano.plot.df$color, levels = c("UP-regulated", "Down-regulated", "Non-DEG"))
# Define fill
volcano.plot.df$fill <- ifelse(volcano.plot.df$fdr < 0.0001 & volcano.plot.df$log2FC < -1, "green",
                               ifelse(volcano.plot.df$fdr < 0.0001 & volcano.plot.df$log2FC > 1, "red",
                                      "grey70"))

# Label genes with high significance
####################################
# Select labels with -log10(p-value) > 200
volcano.plot.df$label <- ifelse(volcano.plot.df$logpval > 200, volcano.plot.df$gene, "")

# Use shape to show some Gene Ontologies
########################################
# Simulate Gene Ontology 1
Gene.Ontology.1 <- sample(volcano.plot.df$gene[which(volcano.plot.df$color == "Down-regulated" &
                                                       volcano.plot.df$logpval > 100 &
                                                       volcano.plot.df$logpval <= 200)],
                          10, replace = FALSE)
# Simulate Gene Ontology 2
Gene.Ontology.2 <- sample(volcano.plot.df$gene[which(volcano.plot.df$color == "UP-regulated" &
                                                       volcano.plot.df$logpval > 50 &
                                                       volcano.plot.df$logpval <= 300)],
                          10, replace = FALSE)
# Assign shape to gene ontolgies
volcano.plot.df$shape <- ifelse(volcano.plot.df$gene %in% Gene.Ontology.1, "GO:00001",
                                ifelse(volcano.plot.df$gene %in% Gene.Ontology.2, "GO:00002",
                                       "NO_GO"))
volcano.plot.df$shape <- factor(volcano.plot.df$shape,
                                levels = c("GO:00001", "GO:00002", "NO_GO"))

# Volcano plot
##############
axis.size <- 1
font.size <- 10

ggplot(volcano.plot.df, aes(x=log2FC, y = logpval, label = label)) +
  geom_vline(xintercept = 0, color = "black",
             linetype = "solid", size = 1) +
  geom_point(aes(color = color, shape = shape),  #aes(shape = shape, color = color),
             fill = alpha(volcano.plot.df$fill, 0.2),
             size = 2, stroke = 1) +
  geom_text_repel(min.segment.length = unit(0, 'lines'),
                  color = volcano.plot.df$fill, size = 4) +
  annotate("text", x = -0.65, y = 280,
           label = 'bold("Ws 90")',
           color = "black", parse = T) +
  annotate("text", x = 0.85, y = 280,
           label = 'bold("mkp1 90")',
           color = "black", parse = T) +
  scale_shape_manual(breaks = c("GO:00001", "GO:00002"),
                     values = c("GO:00001" = 22,
                                "GO:00002" = 24,
                                "NO_GO" = 21)) +
  scale_color_manual(breaks = c("UP-regulated", "Down-regulated"),
                     values = c("UP-regulated" = "red",
                                "Down-regulated" = "green",
                                "Non-DEG" = "grey70")) +
  labs(color = "", shape = "",
       x = expression(log[2]("Fold Change")),
       y = expression(-log[10](italic(p)-value))) +
  theme_classic() +
  theme(axis.line = element_line(size = axis.size, color = "black"),
        axis.ticks = element_line(size = axis.size, color = "black"),
        axis.ticks.length = unit(axis.size * 5, "points"),
        plot.title = element_text(hjust = (0.5), size = font.size + 8),
        axis.title.y = element_text(size = font.size + 5),
        axis.title.x = element_text(size = font.size + 5),
        axis.text = element_text(size = font.size + 2))

ggsave("volcano.plot.png")

################################################################################
# Heatmap
################################################################################
# Counts of the identified DEG
DEG.names <- as.character(rownames(res.sig))
cts.TMM.DEG <- cts.TMM[DEG.names,]

# Generate the plot
png("DEG.heatmap.png")
heatmap.2(as.matrix(cts.TMM.DEG), dendrogram = "column", scale = "row", 
          labRow = FALSE, trace = "none")
dev.off()

################################################################################
# PCA
################################################################################
# calculate PCAs
pca <- prcomp(t(cts.TMM.DEG))

# Generate biplot
ggbiplot(pca, obs.scale = 1, var.scale = 1, var.axes = F,
         labels = colnames(cts.TMM.DEG),
         groups = c(rep("Ws_90",3), rep("mkp1_90",3)), ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme_bw()

ggsave("DEG.PCA.pmg")
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
# Gene Ontology Enrichment Analysis 
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


