setwd("~/Desktop/Sheff/NEOF-Zuloagaea/FunctionalAnnotation/BGC_permutation_test")

#### Get all genes that are within a BGC (biosynthetic and other genes)
# Load BGC table
bgc <- read.delim("plantgeneclusters.txt", header = FALSE, stringsAsFactors = FALSE)
# The last column contains all gene IDs in each BGC
gene_lists <- bgc[, ncol(bgc)]
# Split the semicolon-separated lists and remove underscores if needed
bgc_genes <- unlist(strsplit(gene_lists, ";"))
# Unique genes
bgc_genes <- unique(bgc_genes)

#### Get all the horizontally transferred genes
#Load HGT table
hgt <- read.delim("HGT_genes.txt", header=FALSE)
#Convert it into a vector
hgt_genes <- hgt[[1]]

#### Get all the genes in the genome
#Load gene coordinates and IDs
all_genes_coordinates <- read.delim("gene_coordinates.bed", header=FALSE)
#Convert gene IDs into a vector
all_genes <- all_genes_coordinates[[4]]

#### Remove genes in organelle contigs
#Load organelle contigs
organelle_contigs <- read.delim("organelle_contigs.txt", header=FALSE)
#Convert gene IDs into a vector
organelle_contigs_vector <- organelle_contigs[[1]]
# Keep only genes NOT on organelle scaffolds
bgc_genes_filtered <- bgc_genes[!grepl(paste(organelle_contigs_vector, collapse="|"), bgc_genes)]
hgt_genes_filtered <- hgt_genes[!grepl(paste(organelle_contigs_vector, collapse="|"), hgt_genes)]
all_genes_filtered <- all_genes[!grepl(paste(organelle_contigs_vector, collapse="|"), all_genes)]

#### Perform the permutation test
#Set seed to make it reproducible
set.seed(31)
#Get number of HGTs
n_hgt <- length(hgt_genes_filtered)
#Set the number of permutations to be done, and do it
n_perm <- 10000
perm_counts <- replicate(n_perm, {
  random_genes <- sample(all_genes_filtered, n_hgt, replace = FALSE)
  sum(random_genes %in% bgc_genes_filtered)
})
#Number of observed HGT candidates into BGCs
observed <- sum(hgt_genes_filtered %in% bgc_genes_filtered)
#p-value: proportion of permutations >= observed
p_val <- mean(perm_counts >= observed)

#Plot with ggplot
library(ggplot2)
df <- data.frame(count = perm_counts)
xmax <- max(max(perm_counts), observed)
p <- ggplot(df, aes(x = count)) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "left") +
  geom_vline(xintercept = observed, linetype = "dashed", size = 1, color="red") +
  scale_x_continuous(limits = c(0, xmax),
                     breaks = seq(0, xmax, by = 1)) +
  labs(
    x = "Number of genes in BGCs (n = 237)",
    y = "Frequency"
  ) +
  annotate("text", x = 10, y = 1600, label = "Random gene\n sample", color = "gray40", fontface = "bold") +
  annotate("text", x = 29, y = 1000, label = "HGTs", color = "red", fontface = "bold") +
  theme_classic()
ggsave(
  filename = "HGT_BGC_permutation.png",
  plot = p,
  width = 6,
  height = 4
)
