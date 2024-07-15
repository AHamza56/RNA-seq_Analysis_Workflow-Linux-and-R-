# Installing BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

# Installing DESeq2
BiocManager::install("DESeq2")

# Uploading Feature Count Matrix and Meta-Data 
cnt <- read.csv("counts.csv")
met <- read.csv("metadata2.csv")

# Checking structure of the data
str(cnt)
str(met)

# Making sure the row names in ColData matches to column names in counts data
all(colnames(cnt) %in% rownames(met))

# Checking order of row names and column names
all(colnames(cnt)== rownames(met))

# Loading DESeq2
library(DESeq2)

# Building DESeq data set 
dds <- DESeqDataSetFromMatrix(countData = cnt, colData = met,
                              design = ~dexamethasone)
dds

# Removal of low count reads (optional step)
keep <- rowSums(counts(dds)) >= 10

# Generating subset
dds <- dds [keep, ]
dds

# Setting reference for DEG analysis
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

# Running the DESeq function to performs the differential expression analysis
deg <- DESeq(dds) 

# Extracting the results
res <- results(deg)

# Saving results in CSV format
write.csv(res, "test_udemy.csv")

# Summary statistics of results
summary(res) # It will give results based on p-value < 0.1

# Using p-value 0.05

res0.05 <- results(deg, alpha = 0.05)
summary(res0.05) # now the number up regulated and down regulated genes is different

# Saving latest results in CSV format
write.csv(res0.05, "test_udemy_0.05.csv")

# Converting gene IDs to gene names
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library("AnnotationDbi")
# Transforming results in data frame
res0.05.df <- as.data.frame(res0.05)
str(res0.05.df)
# Creating a new column with gene symbols
res0.05.df$symbol <- mapIds(org.Hs.eg.db, rownames(res0.05.df), keytype = "ENSEMBL" , column = "SYMBOL")
res0.05.df
# Saving in CSV format
write.csv(res0.05.df, "final_test_udemy.csv")

### Quality check parameters 
# Building PCA plot
library(ggplot2)
vsd <- vst(deg, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "dexamethasone", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Using ggplot2 to create PCA plot
p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = dexamethasone )) +
  geom_point(size = 3) + # Changing size to make points larger
  xlab(paste0("PC1: ", percentVar [1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  theme_bw() + # Using a clean theme
  theme(
  axis.title = element_text(size = 14), # Increasing axis title size
  axis.text = element_text(size = 12), # Increasing axis text size
  legend.title = element_text(size = 12), # Increasing legend title size
  legend.text = element_text(size = 10) # # Increasing legend text size
  ) +
  ggtitle("PCA of Variance Stabilized Data") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) # Customizing title

# Printing the plot
print(p)

# Saving the plot as high-resolution images
ggsave("PCA_plot.png", plot = p, width = 8, height = 6, dpi = 300) # Saving as PNG
ggsave("PCA_plot.pdf", plot = p, width = 8, height = 6) # Saving as PDF

# Size factor estimation
size_factors <- sizeFactors(deg)
print(size_factors) # Sample 4 has very low size factor estimation 

# Estimating the dispersion
plot_dispersion <- plotDispEsts(deg)

# Building the mean_average plot
library("ggplot2")
# res0.05 is our DESeq2 result object so converting it to data frame for ggplot
res_df <- as.data.frame(res0.05)
res_df$gene <- rownames(res_df)  # Add row names as a column

# Adding a column for significant genes
res_df$significant <- ifelse(res_df$padj < 0.05, "Significant", "Not Significant")

# Create the MA plot
ma_plot <- ggplot(res_df, aes(x = log2(baseMean), y = log2FoldChange, color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  theme_minimal() +
  labs(
    title = "MA Plot of Differential Expression",
    x = "Log2 Mean Expression",
    y = "Log2 Fold Change",
    color = "Gene Significance"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8), 
    panel.background = element_rect(fill = "white"), 
    plot.background = element_rect(fill = "white")
  )

print(ma_plot)

# Save the plot to a file
ggsave("MA_Plot.png", plot = ma_plot, width = 10, height = 8, dpi = 300)
ggsave("MA_Plot.pdf", plot = ma_plot, width = 10, height = 8, dpi = 300)
ggsave("MA_Plot.tiff", plot = ma_plot, width = 10, height = 8, dpi = 600)

# Getting idea about best genes
install.packages("dplyr")
library(dplyr)
best_genes <- res0.05.df %>%
  arrange(padj) %>% 
  head(10)
write.csv(best_genes, "best_genes.csv")
   
# Generating volcano plot
vol <- res0.05.df %>%
  filter(! is.na(padj))
library(ggplot2)
vol <- res0.05.df %>%
  filter(!is.na(padj)) %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )
# Define a custom theme
custom_theme <- theme_minimal() +
  theme(
    text = element_text(family = "Arial", size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )
# Creating the volcano plot
volcano_plot <- ggplot(vol, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black")) +
  labs(
    title = "Volcano Plot",
    x = expression(Log[2]*" Fold Change"),
    y = expression(-Log[10]*" Adjusted P-value"),
    color = "Significance"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  custom_theme + 
  geom_text(data = best_genes, aes(label = symbol), hjust = -0.2, vjust = 0.5)

# Saving the plot
ggsave("Volcano_plot.tiff", plot = volcano_plot, width = 10, height = 8, dpi = 600)

# Building heat map
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

# Getting top 30 genes
top_genes <- res0.05.df %>%
  arrange(padj)%>%
  head(30)

# Getting the normalized counts from the DESeq2 object (deg) of top 30 genes
mat <- counts(deg, normalized = T) [row.names(top_genes),]
head(mat, 5)

# Computing Z value
mat.z <- t(apply(mat, 1, scale))
head(mat.z, 5)
# Fixing column names
colnames(mat.z) <- rownames(met)
head(mat.z, 5)

# Creating heat map
heat_map <- Heatmap(mat.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat.z),row_labels = top_genes$symbol)

# Saving the heat map to a file
pdf("heatmap1.pdf", width = 10, height = 8)
draw(heat_map, heatmap_legend_side = "bottom")
dev.off()

# Installing necessary packages 
BiocManager::install("clusterProfiler")

# Loading the necessary libraries 
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Extracting gene symbols with log2FoldChange > 0.5
genes_to_test <- rownames(res0.05.df[res0.05.df$log2FoldChange > 0.5,])

# Performing GO enrichment analysis
go_results <- enrichGO(gene = genes_to_test, 
                       OrgDb = org.Hs.eg.db, 
                       keyType = "ENSEMBL", 
                       ont = "BP", 
                       pAdjustMethod = "BH",  # Adjust p-values using Benjamini-Hochberg method
                       pvalueCutoff = 0.05,   # Cutoff for p-values
                       qvalueCutoff = 0.10,   # Cutoff for q-values
                       readable = TRUE)       # Convert IDs to readable format

# Converting enrichment results to data frame for viewing
go_results_df <- as.data.frame(go_results)

# Ploting the bar plot for GO results, showing top 10 categories
fit <- plot(barplot(go_results, showCategory = 10))

# Saving the plot as a PNG file
png("GO_Barplot.png", res = 250, width = 1200, height = 1000)
print(fit)
dev.off()

# Performing KEGG pathway analysis

# Loading the necessary libraries for KEGG pathway analysis and visualization
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)

# Read the input data from a CSV file
final_test_udemy <- read.csv("final_test_udemy.csv")

# Extracting gene symbols from the input data
gene_symbols <- final_test_udemy$symbol

# Converting gene symbols to Entrez IDs using the org.Hs.eg.db database
entrez_ids <- bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# Merging the Entrez IDs with the original data based on the gene symbol
final_data <- merge(final_test_udemy, entrez_ids, by.x="symbol", by.y="SYMBOL", all.x=TRUE)

# Removing rows with NA Entrez IDs
final_data <- final_data[!is.na(final_data$ENTREZID), ]

# Removing duplicated gene symbols to ensure each gene is unique
final_data <- final_data[!duplicated(final_data$symbol), ]

# Creating a named vector of log2 fold changes with Entrez IDs as names
gene_list <- final_data$log2FoldChange
names(gene_list) <- final_data$ENTREZID

# Performing KEGG pathway enrichment analysis using the Entrez IDs
kegg_enrich <- enrichKEGG(gene = names(gene_list),
                          organism = 'hsa',
                          pvalueCutoff = 0.05)

# Displaying the top results of the KEGG pathway enrichment analysis
head(kegg_enrich)

# Generating a dot plot to visualize the top 10 enriched KEGG pathways
dotplot(kegg_enrich, showCategory=10)

# Visualizing a specific KEGG pathway (e.g., "hsa04110") using the pathview package
pathview(gene.data = gene_list, pathway.id = "hsa04110", species = "hsa")