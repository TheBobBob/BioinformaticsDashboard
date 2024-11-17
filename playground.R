
genecounts <- read_excel("/Users/navan/Downloads/RANBP1/RANBP1/240531_combined_count_files_for_r_import_RS4.xlsx", sheet =1, col_names = TRUE)
genecounts <- as.data.frame(genecounts) 
rownames(genecounts) <- genecounts[,1]
genecounts$Gene_Name <- NULL
condition <- data.frame(c('C','C','C','C','C','C','R','R','R','R','R','R'), row.names = c('C1','C2','C3','C4','C5','C6','R1','R2','R3','R4','R5','R6'))
names(condition)[1]<-paste("genotype")

dds <- DESeqDataSetFromMatrix(countData = genecounts, colData = condition, design = ~genotype)

de <- DESeq(dds)
res <- results(de)
#############################################################

res$pvalue_log10 <- -log10(res$pvalue)

pvalue_threshold <- 0.05
fold_change_threshold <- 2
fold_change_min <- -2
fold_change_max <- 2

res$significance <- ifelse(res$pvalue < pvalue_threshold,
                           "Significant", "Not Significant")

res$new_column <- rownames(res)

res$names <- ifelse(
  res$log2FoldChange & res$pvalue & 
    (res$log2FoldChange > fold_change_min & res$log2FoldChange < fold_change_max) &
    (res$pvalue < pvalue_threshold) & res$pvalue_log10 > 100, 
  res$new_column, 
  ""
)

# Set "UP" if log2FoldChange > 0
res$diffexpressed[res$log2FoldChange > 0] <- "UP"

# Set "DOWN" if log2FoldChange < 0
res$diffexpressed[res$log2FoldChange < 0] <- "DOWN"

# Optionally set "NO" if log2FoldChange == 0 or p-value is above threshold (insignificant)
res$diffexpressed[res$pvalue > pvalue_threshold | res$log2FoldChange == 0] <- "NO_CHANGE"
############################################################
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

original_gene_list <- res$log2FoldChange

names(original_gene_list) <- res$new_column

gene_list <- na.omit(original_gene_list)

gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

df_inverted<- gse@result %>%
  separate_rows(core_enrichment, sep = "/")


pvalue_threshold <- 0.05
y_threshold <- -log10(pvalue_threshold)
x_threshold <- 2

saveRDS(list(genecounts = genecounts, dds = dds, res = res, df_inverted = df_inverted), 
        file = "analysis_results.rds")


##################################################
#checking different plots

library(enrichplot)
dotplot(gse, showCategory = 35) #alter size later

# Set the readable GSEA result
gsex <- setReadable(gse, 'org.Hs.eg.db', 'ENTREZID')
geneList <- gse@geneList

# Create the plots with adjusted max.overlaps
p1 <- cnetplot(gsex, foldChange = geneList, max.overlaps = 50)  # Adjust as needed
p2 <- cnetplot(gsex, categorySize = "pvalue", foldChange = geneList, max.overlaps = 50)
p3 <- cnetplot(gsex, foldChange = geneList, circular = TRUE, colorEdge = TRUE, max.overlaps = 50)

# Combine the plots
cowplot::plot_grid(p1, p2, p3, ncol = 3, labels = LETTERS[1:3], rel_widths = c(.8, .8, 1.2))

# Create heatmaps with adjusted max.overlaps
p1 <- heatplot(gsex, showCategory=5)
p2 <- heatplot(gsex, foldChange=geneList, showCategory=5)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

BiocManager::install("clusterProfiler")
upsetplot(gse)

############################
library(rentrez)
library(dplyr) 
library(ggplot2)  

results <- data.frame(Term = character(), Year = integer(), Count = integer(), stringsAsFactors = FALSE)

terms <- gse$Description[1:5]  

for (term in terms) {
  for (year in 2023:2024) {
    query <- paste(term, "[Title/Abstract] AND", year, "[PDAT]")
    
    search_results <- entrez_search(db = "pubmed", term = query, retmax = 0)
    
    results <- rbind(results, data.frame(Term = term, Year = year, Count = search_results$count))
  }
}

total_counts <- results %>%
  group_by(Year) %>%
  summarize(Total_Count = sum(Count), .groups = 'drop')

results <- results %>%
  left_join(total_counts, by = "Year")

results <- results %>%
  mutate(Ratio = Count / Total_Count)

print(results)

ggplot(results, aes(x = Year, y = Ratio, color = Term)) +
  geom_line() +
  geom_point(size = 3, shape = 20, fill = "white", stroke = 1) +  # Bolded dots
  scale_x_continuous(limits = c(2013, 2025), breaks = seq(2013, 2025, by = 2.5)) +  # 2.5-year breaks
  labs(title = "Publication Ratio for Enriched Terms", x = "Year", y = "Publication Ratio") +
  theme_minimal()

ggsave("pubmed_analysis_plot.png", width = 12, height = 8, dpi = 300)

############################################
#prep script

#1. need to take a tsv file, rename the columns depending on whether it was even or odd, with the first half being wild type and the second half being mutant 
#2. rename the rows from ENGS ids to get the actual name of the gene 
#3. return an xlsx file
# Load necessary library
library(readr)

# Function to extract row names and save to a text file
extract_rownames <- function(input_file, rownames_file) {
  # Read the TSV file
  data <- read_tsv(input_file)
  
  # Extract row names (assuming they are the first column)
  row_names <- data[[1]]
  
  # Save row names to a text file
  writeLines(row_names, rownames_file)
}

# Specify input and output file paths
input_file <- "C:/Users/navan/Downloads/RANBP1/RANBP1/cseed_without_preprocessing.tsv"  # Replace with your TSV file path
rownames_file <- "C:/Users/navan/Downloads/RANBP1/RANBP1/names.txt"  # Replace with desired output path

# Call the function
extract_rownames(input_file, rownames_file)

#################

library(biomaRt)
library(org.Hs.eg.db)
#library(EnsDb.Hsapiens.v86)
library(tidyverse)

ensembl.ids <- read.delim('names.txt', header = F)


listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)

ensembl.con <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl')

attr <- listAttributes(ensembl.con)


mappings <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
      filters = "ensembl_gene_id", 
      values = ensembl.ids, 
      mart = ensembl.con)

####################
library(readr)
library(writexl)

# Replace these with your actual file paths
tsv_file_path <- "C:/Users/navan/Downloads/RANBP1/RANBP1/cseed_without_preprocessing.tsv"
row_names_file_path <- "C:/Users/navan/Downloads/RANBP1/RANBP1/names.txt"
output_excel_file_path <- "C:/Users/navan/Downloads/RANBP1/RANBP1/final_cseed_output.xlsx.xlsx"

# Load the TSV file, treating the first row as the header
df <- read_tsv(tsv_file_path)

# Load the row names from the text file
new_row_names <- readLines(row_names_file_path)

# Check if the number of row names matches the number of rows in the DataFrame
if (length(new_row_names) != nrow(df)) {
  stop("The number of row names does not match the number of rows in the DataFrame.")
}

# Replace the row names
rownames(df) <- new_row_names

# Check if the number of columns is even, and trim if necessary
num_cols <- ncol(df)
if (num_cols %% 2 != 0) {
  df <- df[, -ncol(df)]  # Remove the last column if odd
}

# Rename the first half of the columns to 'C' and the second half to 'R'
half_cols <- ncol(df) / 2
colnames(df) <- c(rep("C", half_cols), rep("R", half_cols))

# Write to an Excel file
write_xlsx(df, output_excel_file_path)

cat("Process completed successfully!\n")


########################################################################################################
#need to get rid of duplicates for rnaseq and proteomics dataset
genecounts <- read_excel("C:/Users/navan/Downloads/RANBP1/RANBP1/cseed_data/PROTEOMICS_final_cseed_output.xlsx")
genecounts <- as.data.frame(genecounts)
# Step 1: Convert genecounts to a data frame (if it's not already one)
genecounts <- as.data.frame(genecounts)

# Step 2: Set the first column as rownames for genecounts
rownames(genecounts) <- genecounts[, 1]

# Step 3: Remove the first column from genecounts (it's now the row names)
genecounts <- genecounts[, -1]  # This removes the first column (which has been set as rownames)

# Step 4: Create descriptions from the second column of genecounts
descriptions <- data.frame(Description = genecounts[, 1])  # Create descriptions dataframe from second column

# Step 5: Set the rownames of descriptions to be the same as genecounts
rownames(descriptions) <- rownames(genecounts)

# Step 6: Remove the second column from genecounts (it is now in descriptions)
genecounts <- genecounts[, -1]  # Remove the second column which is now in descriptions

# Assuming df_final_rnaseq is your data frame
df_rnaseq <- df_final_rnaseq[!duplicated(df_final_rnaseq[, 1]) & !is.na(df_final_rnaseq[, 1]), ]

write.xlsx(df_rnaseq, "C:/Users/navan/Downloads/RANBP1/RANBP1/cseed_data/PROTEOMICS_final_cseed_output.xlsx")



##################################
library(UpSetR)

gse_df <- gse@result
top_terms <- gse_df %>% arrange(pvalue) %>% head(10)  # Adjust based on your data
top_gene_sets <- strsplit(top_terms$core_enrichment, "/")
gene_sets_list <- lapply(top_gene_sets, function(x) unique(trimws(x)))

gene_sets_df <- fromList(setNames(gene_sets_list, top_terms$ID))

# Create the UpSet plot
upset(gene_sets_df, 
      sets = names(gene_sets_df), 
      main.bar.color = "steelblue",
      sets.bar.color = "darkred", 
      order.by = "freq", 
      matrix.color = "gray",
      keep.order = TRUE)

genecounts <- genecounts[, -1]
genecounts <- genecounts[, -1]
descriptions <- genecounts[,1]

num_samples <- ncol(genecounts)
