library(Seurat)

set.seed(666)

# perform differential expression (DE) testing using Seurat on the positive and the negative cells
de_testing <- function(positive_cells, negative_cells) {
    markers <- FindMarkers(pbmc, ident.1=positive_cells, ident.2=negative_cells)
    markers <- markers[markers$p_val_adj< 0.05,] # select genes with adjusted p-value < 0.05
    markers <- markers[order(markers$avg_log2FC,decreasing = TRUE),] # sort by average avg_log2FC
    return(markers)
}

# output the table from (DE) analysis with a LaTex formatting
output_markers <- function(markers, ind_min, ind_max) {
    for(ind_row in ind_min:ind_max) {
        line <- ""
        line <- paste(line, rownames(markers)[ind_row])
        for(ind_col in 1:5) {
            line <- paste(line, "&")
            val <- markers[ind_row,ind_col]
            line <- paste(line, "$", format(val, digits=3), "$")
        }
        line <- paste(line, "\\\\ \\hline\n")
        cat(line)
    }
}


# load the raw counts
raw_counts <- read.csv("../../dataset/Imagine/raw_matrix.csv", header=FALSE)
colnames(raw_counts) <- as.character(unlist(raw_counts[1,]))
rownames(raw_counts) <- as.character(unlist(raw_counts[,1]))
raw_counts <- raw_counts[-c(1),-c(1)]
pbmc <- CreateSeuratObject(counts = raw_counts, project = "pbmc")
pbmc <- NormalizeData(pbmc)


###################################################################################
# (DE) analysis for cluster 20 (dynamics)
df_selection <- read.csv("selected_cells_dynamics_c20.csv")
print(df_selection)
positive_cells <- df_selection[df_selection$selection == 1,'X']
negative_cells <- df_selection[df_selection$selection == 0,'X']
markers <- de_testing(positive_cells, negative_cells) # perform DE testing
print(head(markers)) # print the first rows
write.table(markers, 'c20_dynamics_markers.csv', sep=",", row.names = TRUE, quote=FALSE)
