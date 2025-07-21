library(DESeq2)
library(stringr)
library(ggplot2)
library(dplyr)

packageVersion("DESeq2")


cts <- read.csv("counts_gene.txt", sep="\t", header=TRUE, row.names=1)

meta <- read.csv("meta3.csv", header=TRUE)

to_int <- function(x){
  as.integer(x*1000)
}
counts <- select(cts, meta$ID)
cts_count <- apply(counts,c(1,2),as.integer)
cts_count[is.na(cts_count)] <- 0

dds <- DESeqDataSetFromMatrix(countData = cts_count,
                              colData = meta,
                              design=~group)

dds <- DESeq(dds)
res <- results(dds)

summary(res)
res
write.csv(res, "mockvs24_DE_result_gene_v2.csv")
df = as.data.frame(res)

df %>% filter(padj < 0.01)

df$GENE = row.names(df)
head(df)
plotMA(res, ylim=c(-2,2))

df %>% filter(grepl('ISG15|MX1|OAS2|TRIM22|IFITMx', GENE))
ISG_trans <- rownames(df %>% filter(grepl('ISG15|MX1|OAS2|TRIM22|IFITMx', GENE)))
length(ISG_trans)

df$diffexpressed <- "NO"
# if log2Foldchange > 2 and pvalue < 0.01, set as "UP" 
df$diffexpressed[df$log2FoldChange > 2 & df$padj < 0.01] <- "UP"
# if log2Foldchange < -2 and pvalue < 0.01, set as "DOWN"
df$diffexpressed[df$log2FoldChange < -2 & df$padj < 0.01] <- "DOWN"
df$dflabel <- NA
df$dflabel[df$GENE %in% ISG_trans] <- df$GENE[df$GENE %in% ISG_trans]
p <- ggplot(data=df, aes(x=-log2FoldChange, y=-log10(padj), col=diffexpressed, label=dflabel)) + geom_point(alpha=0.5) + 
  theme_minimal() +
  geom_text()
  #geom_text_repel() +
  #scale_color_manual(values=c("blue", "black", "red"))
p
p2 <- p + geom_vline(xintercept=c(-2, 2), col="yellow") +
  geom_hline(yintercept=-log10(0.01), col="yellow")
p2
