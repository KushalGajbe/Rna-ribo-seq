# Load the DESeq2 library
library(DESeq2)

# Read count data and meta data
countData <- read.csv("final_counts_rna.csv", row.names = "Geneid")
colData <- read.csv("Metadata_rna.csv", header = TRUE, sep = ",", row.names = 1)

# Read TPM (Transcripts Per Million) data
tpm_data <- read.csv("TPM.csv", header = TRUE, row.names = 1)

# Check if all row names in meta data are present in count data
print(all(rownames(colData) %in% colnames(countData)))

# Subset count data to keep only columns with matching row names in meta data
countData <- countData[, rownames(colData)]

# Check if row names in meta data and count data match exactly
print(all(rownames(colData) == colnames(countData)))

# Filter count data based on TPM condition for both conditions
filtered_countData <- countData[rowSums(tpm_data[, c("Fluconazole_rna1", "Fluconazole_rna2", "Fluconazole_rna3", "Fluconazole_rna4")] > 1) >= 3 &
                                  rowSums(tpm_data[, c("Nodrug_rna1", "Nodrug_rna2", "Nodrug_rna3", "Nodrug_rna4")] > 1) >= 3,]

# Remove rows with NA values
filtered_countData<- na.omit(filtered_countData)

# Replace zero counts with a small positive value (1 in this case)
filtered_countData[filtered_countData == 0] <- 1

# Convert the 'Condition' variable in meta data to a factor
colData$Condition <- as.factor(colData$Condition)


# Create a DESeqDataSet object
design_formula <- ~ factor(Condition)
dds <- DESeqDataSetFromMatrix(countData = filtered_countData, colData = colData, design = design_formula)

# Run DESeq normalization and differential expression analysis
dds <- DESeq(dds)
res <- results(dds, tidy = TRUE)

write.csv(res, file="raw_results.csv", row.names = FALSE)

# Filter for differentially expressed genes
significant_results <- subset(res, abs(log2FoldChange) >= 1 & padj <= 0.05)

# Remove rows with NA values
significant_results <- na.omit(significant_results)

# Write the results to a CSV file
write.csv(significant_results, file = "results1.csv", row.names = FALSE)




##############################################################################

library(pheatmap)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(dplyr)

#Dispersion Plot
plotDispEsts(dds)

#PCA plot
#variance stabalizing transformation
vsd<- vst(dds, blind = FALSE)

#use transformed values to generate pca plot
plotPCA(vsd, intgroup= c("Condition"))


#heatmaps

#Heatmap of sample to sample distance matrix (with clustering) based on the normalised counts
#generate distance matrix 
sampleDist<-dist(t(assay(vsd)))
sampleDistMatrix<- as.matrix(sampleDist)

#set a color scheme

colors<- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)

#generate the heatmap

pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDist, 
         clustering_distance_cols = sampleDist, col=colors)

#heatmap of log transformed normalised counts. We will use the top 10 genes. 



top_hits<-res_df[order(res_df$padj),][1:10,]
top_hits<-row.names(top_hits)
top_hits


#plot MA
plotMA(dds,ylim=c(-2,2))

#remove noise 

resLFC<- lfcShrink(dds, coef ="factor.Condition.NoDrug"
                   , type="apeglm")

# use resultsNames(dds), to fins which coeff
# need to install apeglm package to run the above code. 
#BiocManager::install("apeglm")

plotMA(resLFC, ylim=c(-2,2))




#Volcano Plot

#cahnge resLFC to a dataframe
resLFC<- as.data.frame(resLFC)

#LABEL THE GENES
resLFC$diffexpressed<-"NO"

resLFC$diffexpressed[resLFC$log2FoldChange>0.1 & resLFC$padj<0.05]<-"UP"
resLFC$diffexpressed[resLFC$log2FoldChange<0.1 & resLFC$padj<0.05]<-"DOWN"
resLFC$delabel<-NA

dev.off()
ggplot(data=resLFC, aes(x=log2FoldChange,y=-log10(pvalue),col=diffexpressed,label=delabel))+
  geom_point()+
  theme_minimal()+
  geom_text_repel()+
  scale_color_manual(values = c('blue','black','red'))+
  theme(text = element_text(size=20))



