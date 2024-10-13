library(pheatmap)
library(ComplexHeatmap)
library(reshape2)

# SET UP
# Genes of interest
genes_of_interest <- c("ATM", "TRIM28", "TOP1", "TOP2A", "TOP2B", "TP53", "TP73", "TP63")
# Desired order of genes for heatmap
desired_order <- c("TP63", "TP73", "TP53", "TOP1", "TOP2B", "TRIM28", "ATM", "TOP2A")
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
ann_colors = list( Mark = c(H3K9Ac="green", H3K9me3="orange"))

create_matrix <- function(files) {
  data_list <- lapply(files, function(file) {
    df <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    df_subset <- df[, c(4, 5)]
    colnames(df_subset) <- c("gene_name", "value")
    df_filtered <- df_subset[df_subset$gene_name %in% genes_of_interest, ]
    df_filtered$value <- as.numeric(as.character(df_filtered$value))
    df_filtered <- df_filtered[!is.na(df_filtered$value), ]
    df_max <- aggregate(value ~ gene_name, data = df_filtered, max)
    file_info <- sub("\\..*", "", basename(file))
    replicate_info <- sub(".*\\.(H3K9[^\\.]+\\.rep\\d+).*", "\\1", basename(file))
    combined_name <- paste0(file_info, ".", replicate_info)
    df_max$combined_name <- combined_name
    return(df_max)
  })

  combined_df <- do.call(rbind, data_list)
  combined_df_matrix <- dcast(combined_df, gene_name ~ combined_name, value.var = "value", fun.aggregate = max, na.rm = TRUE)

  # Set row names make sure all entries are #
  rownames(combined_df_matrix) <- combined_df_matrix$gene_name
  combined_df_matrix <- combined_df_matrix[, -1]  # Remove gene_name column
  combined_df_matrix <- as.matrix(combined_df_matrix)

  # Print row names and dimensions
  print("Row names of combined_df_matrix:")
  print(rownames(combined_df_matrix))
  print("Dimensions of combined_df_matrix:")
  print(dim(combined_df_matrix))

  return(combined_df_matrix)
}



#RUN
rm(h3k9ac_matrix)
rm(h3k9me_matrix)
h3k9ac_path <- "promoter_5kb"
h3k9me_path <- "promoter_5kb"
h3k9ac_csv <- "promoter_5kb.h3k9ac.heatmap.csv"
h3k9me_csv <- "promoter_5kb.h3k9ac.heatmap.csv"
basename <- "H3K9Ac_H3K9me3.promoter_5kb_max"

h3k9ac_files <- list.files(path = h3k9ac_path, full.names = TRUE, pattern = "H3K9ac")
h3k9me_files <- list.files(path = h3k9me_path, full.names = TRUE, pattern = "H3K9me")

# Create matrices
h3k9ac_matrix <- create_matrix(h3k9ac_files)
h3k9me_matrix <- create_matrix(h3k9me_files)
tmp <- gsub(".H3K9ac", "", colnames(h3k9ac_matrix))
tmp <- gsub("rep", "", tmp)
colnames(h3k9ac_matrix) <- tmp
gsub(".H3K9me3", "", colnames(h3k9me_matrix))
tmp <- gsub("rep", "", tmp)
colnames(h3k9me_matrix) <- tmp

# Combine matrices side by side
combined_matrix <- cbind(h3k9ac_matrix, h3k9me_matrix)

## Custom column names to differentiate the sections
#colnames(combined_matrix) <- c(paste("H3K9Ac", colnames(h3k9ac_matrix), sep = "_"), 
#                               paste("H3K9Me", colnames(h3k9me_matrix), sep = "_"))
#
## Reorder genes in the combined matrix
#combined_matrix <- combined_matrix[match(desired_order, rownames(combined_matrix)), ]

# CSV files with values
write.csv(as.data.frame(h3k9ac_matrix), file = h3k9ac_csv, row.names = TRUE)
write.csv(as.data.frame(h3k9me_matrix), file = h3k9me_csv, row.names = TRUE)


# Combined Heatmap without clustering
#png(filename=paste(basename, ".v3.png", sep=""), width=2000, height=500, pointsize=100)
#tiff(filename=paste(basename, ".v4.tiff", sep=""), width=2000, height=500, pointsize=100)
tiff(filename=paste(basename, ".v4.tiff", sep=""), width=1000, height=300, res=100)
pheatmap(combined_matrix, 
         cluster_cols = FALSE,  # Disable column clustering
         scale = "none",  # Disable scaling
         breaks = seq(0, 1, length.out = 101),  # Define breaks for colors
         show_rownames = TRUE,  # Ensure row names (gene names) are shown
         show_colnames = TRUE,  # Ensure column names are shown
         gaps_col=26,
         annotation_col=data.frame(Mark=c(rep("H3K9Ac", 26), rep("H3K9me3", 26))),
         annotation_colors=ann_colors,
         fontsize = 10,
         color = color_palette)  # Apply custom color palette
dev.off()

