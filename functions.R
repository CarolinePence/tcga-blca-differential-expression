#### nicePCA function ####
# makes a PCA plot using pca output and the group2 logical vector 
nicePCA <- function(pcaData, groupLogical){
    
    # pull out PC1 and PC2 scores
    pc1 <- pcaData$x[,1]
    pc2 <- pcaData$x[,2]
    
    # make a group column from the logical vector
    # FALSE = group1, TRUE = group2 
    groupLabs <- factor(groupLogical, levels = c(FALSE, TRUE), labels = c("1", "2"))
    
    # put this info all in a data frame 
    plotDf <- data.frame(PC1 = pc1, PC2 = pc2, Group = groupLabs)
    
    # get percent variance explained
    pve1 <- round(pcaData$pve[1] * 100, 2)
    pve2 <- round(pcaData$pve[2] * 100, 2)
    
    # make the plot
    p <- ggplot(plotDf, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
        geom_point(size = 2) +
        coord_fixed(ratio = 1) +
        labs(x = paste0("PC1 (", pve1, "% variance explained)"),
             y = paste0("PC2 (", pve2, "% variance explained)")) +
        theme_minimal()
    
    return(p)
}





#### niceVolcano function ####
# makes a volcano plot from resOutput 
niceVolcano <- function(resData, pCutoff = 0.05, fcCutoff = 1.5){
    df <- as.data.frame(resData)
    
    # calculate -log10 of adjusted p values (y-axis)
    df$negLog10P <- -log10(df$padj)
    
    # make column of genes that are significant (less than p)
    df$isSig <- !is.na(df$padj) & df$padj < pCutoff
    
    # figure out which genes pass both the p value and fold change cutoffs
    df$passesBoth <- df$isSig & (df$log2FoldChange > fcCutoff | df$log2FoldChange < -fcCutoff)
    
    # label genes that pass both cutoffs and keep top 5 to display on the graph
    labelDf <- df[which(df$passesBoth), ]
    labelDf <- labelDf[order(labelDf$padj), ]
    labelDf <- head(labelDf, 5)
    
    # make the plot
    p <- ggplot(df, aes(x = log2FoldChange, y = negLog10P)) +
        geom_point(aes(color = passesBoth, shape = isSig), size = 1, na.rm = TRUE) +
        geom_hline(yintercept = -log10(pCutoff), linetype = "dashed") +
        geom_vline(xintercept = c(-fcCutoff, fcCutoff), linetype = "dashed", color = "red") +
        geom_text_repel(data = labelDf,
                        aes(label = gene_name),
                        size = 3) +
        labs(x = "log2(Fold Change)",
             y = "-log10(FDR-Adjusted p-values)",
             title = paste0(nrow(df), " loci analyzed"),
             color = "Exceeds both cutoffs",
             shape = "Sig. DE") +
        theme_minimal()
    
    return(p)
}





#### resSummary function ####
# prints out a summary of the DESeq2 results
resSummary <- function(resData){
    
    df <- as.data.frame(resData)
    
    message("How many genes were included in the DESeq2 analysis after filtering?")
    print(nrow(df))
    
    message("How many genes included in the DESeq2 analysis had non-missing p-adjusted values?")
    nonNA <- sum(!is.na(df$padj))
    print(nonNA)
    
    message("Note: using 0.05 as the p-adjusted cutoff for significance for the following calculations.")
    
    message("How many genes were significantly downregulated in group2 relative to group1?")
    downreg <- sum(!is.na(df$padj) & df$padj < 0.05 & df$log2FoldChange < 0)
    print(downreg)
    
    message("How many genes were significantly upregulated in group2 relative to group1?")
    upreg <- sum(!is.na(df$padj) & df$padj < 0.05 & df$log2FoldChange > 0)
    print(upreg)
    
    message("Names of the top 10 most significantly differentially expressed genes from smallest to largest p-adjusted:")
    sortedDf <- df[order(df$padj), ]
    top10 <- head(sortedDf$gene_name, 10)
    print(top10)
}
