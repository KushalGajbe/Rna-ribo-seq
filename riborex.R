
#BiocManager::install("edgeR")
#install.packages("fdrtool")


#install.packages("riborex-2.4.0.tar.gz", repos=NULL, type="source")

library(riborex)

RNACntTable<- read.csv("counts_rna.csv", row.names = "Geneid")
RiboCntTable<- read.csv("counts_ribo.csv", row.names = "Geneid")
tpm_data<- read.csv("TPM.csv", row.names = "Geneid")

head(RNACntTable, 5)
head(RiboCntTable, 5)

filtered_rnacnt <- RNACntTable[rowSums(tpm_data[, c("Fluconazole_rna1", "Fluconazole_rna2", "Fluconazole_rna3", "Fluconazole_rna4")] > 1) >= 3 &
                                  rowSums(tpm_data[, c("Nodrug_rna1", "Nodrug_rna2", "Nodrug_rna3", "Nodrug_rna4")] > 1) >= 3,]

filtered_ribocnt<- RiboCntTable[rowSums(tpm_data[, c("Fluconazole_ribo1", "Fluconazole_ribo2", "Fluconazole_ribo3", "Fluconazole_ribo4")] > 1) >= 3 &
                                                    rowSums(tpm_data[, c("Nodrug_ribo1", "Nodrug_ribo2", "Nodrug_ribo3", "Nodrug_ribo4")] > 1) >= 3,]
rnaCond<- c("Fluconazole", "Fluconazole", "Fluconazole", "Fluconazole",
            "Nodrug", "Nodrug", "Nodrug", "Nodrug")
riboCond<- c("Fluconazole", "Fluconazole", "Fluconazole", "Fluconazole",
            "Nodrug", "Nodrug", "Nodrug", "Nodrug")

# Get the common genes between filtered_rnacnt and filtered_ribocnt
common_genes <- intersect(row.names(filtered_rnacnt), row.names(filtered_ribocnt))

# Subset filtered_rnacnt and filtered_ribocnt to keep only the common genes
filtered_rnacnt_common <- filtered_rnacnt[common_genes, ]
filtered_ribocnt_common <- filtered_ribocnt[common_genes, ]


res.deseq2<-riborex(filtered_rnacnt_common, filtered_ribocnt_common, rnaCond, riboCond)


# Filter for differentially expressed genes
significant_results <- subset(res.deseq2, abs(log2FoldChange) >= 1 & padj <= 0.05)

# Remove rows with NA values
significant_results <- na.omit(significant_results)

hist(res.deseq2$pvalue, main = 'DESeq2 unadjusted p-values',
     xlab='Unadjusted p-values')

write.table(res.deseq2, "riborex_res_deseq2.txt", quote=FALSE)
write.csv(significant_results, "riborex_final.csv", quote=FALSE)
