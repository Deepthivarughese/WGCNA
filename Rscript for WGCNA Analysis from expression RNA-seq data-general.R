# Install the impute package from Bioconductor if you haven't already
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("WGCNA")

library(WGCNA)
library(DESeq2)
library(ggplot2)

Expression_data <- read.csv("Mydataset_Expression_Data.csv", row.names = 1)
#head(Expression_data)
#head(Expression_data$ID)
head(Expression_data[-1,])

############## CHECKING IF THE DATA IS NORMALISED OR NOT ################
# Summary statistics
summary(as.vector(as.matrix(Expression_data)))

# Histogram
hist(as.vector(as.matrix(Expression_data)), breaks = 50, main="Expression Distribution", xlab="Expression Value")

# Check for large values (typical of raw counts)
max(Expression_data)


################# NORMALISATION ############
normalized_data <- log2(Expression_data + 1)
head(normalized_data)
write.csv(normalized_data, "Mydataset_Log2_normalized_expressiondata.csv")

#normalized_data <- Expression_data
########### VERIFY NORMALIZATION WITH BOX PLOT #######
boxplot(as.matrix(normalized_data), main="Boxplot of Normalized Expression Data", las=2)
hist(as.matrix(normalized_data), breaks=50, main="Histogram of Normalized Expression Data")


######### LOAD NORMALISED DATA AND CHECK FOR MISSING VALUES ##########
dataExpr <- normalized_data
datExpr <- read.csv("Mydataset_Log2_normalized_expressiondata", row.names = 1)
datExpr <- as.data.frame(t(datExpr))  # Transpose to match WGCNA format
dim(datExpr)

head(datExpr)

########## Check for missing values & REMOVAL OF BAD SAMPLES ###########

gsg <- goodSamplesGenes(datExpr, verbose = 3)
print(gsg$allOK)

# Remove bad genes and samples
datExpr_clean <- datExpr[gsg$goodSamples, gsg$goodGenes]
dim(datExpr_clean)


########### CHECK FOR OUTLIERS ##########
sampleTree <- hclust(dist(datExpr_clean), method = "average")
plot(sampleTree, main = "Sample Clustering to Detect Outliers", sub="", xlab="")


########## REMOVE OUTLIERS BY SETTING CUTOFF HEIGHT ACCRODINGLY ##############

clust <- cutreeStatic(sampleTree, cutHeight = 100, minSize = 10)
keepSamples <- (clust == 1)
datExpr_clean1 <- datExpr_clean[keepSamples, ]
dim(datExpr_clean1)


############ CHOOSE SOFT THRESHOLDING POWER ############
# A crusial step in WGCNA is choosing soft thresholding power
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr_clean1, powerVector = powers, verbose = 5)

# Plot scale-free topology model fit
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit",
     type="b", main="Soft Threshold Selection")

# Print the chosen power
print(paste("Optimal power:", sft$powerEstimate))





############ CONSTRUCT THE NETWORK AND IDENTIFY MODULES ##########

softPower <- 20# Adjust based on the previous step
adjacency <- adjacency(datExpr_clean, power = softPower)

# Turn adjacency matrix into Topological Overlap Matrix (TOM)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# Hierarchical clustering
geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, main = "Clustering of Genes Based on TOM", sub="", xlab="")


########## CONSTRUCT THE NETWORK ##########

net <- blockwiseModules(
  datExpr_clean,
  power = 17,  # Replace with the best power from pickSoftThreshold()
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  corType = "bicor",  # More robust than Pearson
  networkType = "signed",  # More biologically relevant
  saveTOMs = TRUE,  # Save TOM for later use
  verbose = 3
)

moduleColors <- net$colors  # Extract module assignments from WGCNA result

table(moduleColors)

hubGenes <- chooseTopHubInEachModule(datExpr_clean, net$colors)
hubGenes


########### TRAIT - MODULE RELATIONSHIP ##########

#Create a trait data in csv format using script
trait_data <- data.frame(
  row.names = c("1","2","3","5","6","13","14","15","16","17"),
  Disease_Status =c(0,0,0,0,0,1,1,1,1,1)
)
trait_data
rownames(datExpr_clean)
trait_match_data <- trait_data[rownames(datExpr_clean), , drop = FALSE]

MEs = moduleEigengenes(datExpr_clean, colors = net$colors)$eigengenes

moduleTraitCor <- cor(MEs, trait_match_data, use = "p")
moduleTraitP <- corPvalueStudent(moduleTraitCor, nrow(datExpr_clean))
MEs


plotDendroAndColors(net$dendrograms[[1]], net$colors[net$blockGenes[[1]]])



###############  Co-EXPRESSED GENE NETWORK ############

########## To get co-expressed gene set we need expression data that is not transposed. So use the original Exp_data

### Display hub gene in each module
table(moduleColors)
hubGenes  # This will display the modules or the hub genes


# Assuming 'moduleColors' contains module assignments for each gene
genes_in_M5 <- names(moduleColors[moduleColors == 5])

# Print the genes
print(genes_in_M5)

# Taking the hub gene of M5
hub_gene <- "Syt17"

selected_genes <- unique(c(genes_in_M5, hub_gene))  # Combine ME5 genes with hub gene

# Subset expression data
expr_subset <- normalized_data[selected_genes, ]
expr_subset

############ Step COMPUTE the Weighted Co-Expression Matrix ############
# Choose an appropriate soft-thresholding power (determined by analysis)
power <- softPower

# Compute adjacency matrix (weighted co-expression network)
adjacency_matrix <- adjacency(t(expr_subset), power = power)

# Convert adjacency matrix to a data frame
adjacency_df <- as.data.frame(as.table(adjacency_matrix))
colnames(adjacency_df) <- c("Gene1", "Gene2", "Weight")

# Remove self-pairs (Gene1 == Gene2)
adj_df <- adjacency_df[adjacency_df$Gene1 != adjacency_df$Gene2, ]

tail(adj_df)

# Save the weighted co-expression network to a CSV file
write.csv(adj_df, "Mydataset_ME5_hub_gene_Syt17_network.csv", row.names = FALSE)

# Set a threshold to retain strong connections (adjust as needed)
threshold <- 0.7

filtered_network <- adj_df[adj_df$Weight > threshold, ]

# Save filtered network
write.csv(filtered_network, "Filtered_Mydataset_ME5_hub_gene_Syt17_network.csv", row.names = FALSE)


