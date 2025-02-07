#### Reformat counts matrix #####

# To integrate transcriptome data with genotyping and FACS indexing data, we need a shared set of cell IDs that can join the data together.
# Transcriptome cell barcodes are derived from libraries with unique i5/i7 plate indexes and oligo-dT cell barcodes.
# This script joins the transcriptome cell barcodes to the descriptive cell IDs from FACS that match those in the other metadata.

# It then takes the outputs from STARsolo and reformats the transcriptome counts matrix with gene name and cell IDs.
# The resulting counts matrix can be used in single-cell genomics analysis packages.


## Import libraries ####
library(tidyverse)

#### READ IN FILES #####

# Read in the gene features file
features_file <- read_tsv("data/rna/Solo.out/GeneFull_Ex50pAS/raw/features.tsv.gz",
                          col_names = c("gene_id", "gene_name", "assay"))


# Read in the tsv file with cell barcodes:
barcodes_file <- read_tsv("data/rna/Solo.out/GeneFull_Ex50pAS/raw/barcodes.tsv.gz", col_names = "barcode")


# Import the table of transcriptome i5/i7 plate indexes
# This needs have the details of plate indexes used during tagmentation PCR in the transcriptome library preparation.
plate_indexes <- read_csv("data/rna/transcriptome_plate_indexes.csv")

# Import the list of oligo-dT cell barcodes with details of which barcode is found in each well
oligodt_barcodes <- read_csv("data/rna/oligodT_384_barcodes.csv")

# Import the genotyping metadata which has the identity of every well
genotyping_metadata <- read_tsv("data/genotyping/genotyping_metadata.txt")


## Combine plate indexes and single-cell barcodes ####

# Join the cell IDs to transcriptome i5/i7 plate indexes and oligo-dT cell barcodes
transcriptome_indexes <- genotyping_metadata %>% 
  left_join(plate_indexes, by = join_by(Plate)) %>% 
  left_join(oligodt_barcodes, by = join_by(Index, OligodT_barcode_plate))


# Create a column with the cell barcode.
# The barcode from the transcriptome preprocessing pipeline is a combination of the Library (plates with shared i5/i7 plate indexes) and the OligodT_barcode_name
transcriptome_indexes <- transcriptome_indexes %>% 
  select(Cell, Library, i7_index_name, i5_index_name, OligodT_barcode_name) %>% 
  mutate(barcode = paste(Library, OligodT_barcode_name, sep = "_"))


## Reformat the barcodes file to use descriptive cell names ####

# This includes "unknown" cell IDs, which correspond to unassigned counts
barcodes_file <- barcodes_file %>% left_join(transcriptome_indexes) %>% mutate(Cell = coalesce(Cell, barcode))

# Write a barcodes file with the descriptive cell IDs
barcodes_file %>% pull(Cell) %>% write_lines("data/rna/Solo.out/GeneFull_Ex50pAS/raw/cell_id_barcodes.tsv.gz")



## Reformat the RNA-seq counts to make a counts matrix ###

# We add gene names as the rownames and descriptive cell IDs as colnames

# Read in the tsv file with RNA-seq counts:
raw_counts_file <- Matrix::readMM("data/rna/Solo.out/GeneFull_Ex50pAS/raw/matrix.mtx.gz")

# Function for adding cell barcodes and gene names
format_counts_matrix <- function(matrix, barcodes, features){
  if (dim(matrix)[2] != length(barcodes$barcode)) {
    stop("Mismatch dimension between barcodes file and the matrix file")
  }else{
    print("Cell numbers match")
    colnames(matrix) <- barcodes$Cell
    }
  
  if (dim(matrix)[1] != length(features$gene_name)) {
    stop("Mismatch dimension between gene file and the matrix file")
  }else{
    print("Gene numbers match")
    rownames(matrix) <- make.unique(features$gene_name)
  }
  return(matrix)
}

# Add cell barcodes and gene names
counts_matrix <- format_counts_matrix(raw_counts_file, barcodes_file, features_file)


# Note if multiple instances of the pre-processing pipeline have been run, counts matrices can be joined with cbind:
# counts_matrix <- matrix1 %>% 
#   cbind(matrix2) %>% 


# Convert to a dataframe
counts_matrix <- counts_matrix %>% as.matrix() %>% as.data.frame()

dim(counts_matrix)

# Remove the unassigned counts
counts_matrix <- counts_matrix %>% 
  select(-contains("unknown"))
dim(counts_matrix)


# Export this to a file
counts_matrix %>% data.table::fwrite("data/rna/transcriptome_counts_matrix.txt.gz", sep = "\t", row.names = T)

