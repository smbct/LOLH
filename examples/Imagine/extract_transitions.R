library(Seurat)

set.seed(666)

# load the raw counts
raw_counts <- read.csv("../../dataset/Imagine/raw_matrix.csv", header=FALSE)
colnames(raw_counts) <- as.character(unlist(raw_counts[1,]))
rownames(raw_counts) <- as.character(unlist(raw_counts[,1]))
raw_counts <- raw_counts[-c(1),-c(1)]
pbmc <- CreateSeuratObject(counts = raw_counts, project = "pbmc")
pbmc <- NormalizeData(pbmc)

# scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# find variable features
pbmc <- FindVariableFeatures(object = pbmc)

# principal component analysis
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# if there is an error with Mnumeric
# install.packages("SeuratObject")
# -> https://github.com/satijalab/seurat/issues/4436

# computation of the neighborhood graph
pbmc <- FindNeighbors(
  pbmc,
  query = NULL,
  distance.matrix = FALSE,
  k.param = 3,
  prune.SNN = 1/15,
  nn.method = "annoy",
  n.trees = 50,
  annoy.metric = "euclidean",
  nn.eps = 0,
  verbose = TRUE,
  force.recalc = FALSE,
  l2.norm = FALSE,
  cache.index = FALSE,
  index = NULL
)


# extraction of the graph
names(pbmc)
graph <- pbmc@graphs$RNA_nn

dp <- diff(graph@p)
indexes <- cbind(graph@i+1,rep(seq_along(dp),dp))

barcodes <- Cells(pbmc)
left <- sapply(indexes[,1], function(i) barcodes[i])
right <- sapply(indexes[,2], function(i) barcodes[i])

# add the reverse transitions
left_comb <- c(left, right)
right_comb <- c(right, left)

df_transitions <- data.frame(left_comb, right_comb)
df_transitions <- df_transitions[-c(1),]

names(df_transitions) <- c("T-1", "T")

# row names are now T + row index
rownames_transitions <- sprintf("T%s",seq(0,(dim(df_transitions)[1]-1)))
rownames(df_transitions) <- rownames_transitions

write.csv(df_transitions, "../../dataset/Imagine/transitions.csv", row.names=TRUE, quote=FALSE)
