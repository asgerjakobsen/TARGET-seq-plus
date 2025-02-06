# Make genotyping demultiplexing metadata file for use in TARGET-seq pipeline

## Import libraries ####

library(tidyverse)
library(stringr)

## Combine plate and single-cell barcodes ####

# Make a table of the plate barcodes used in the Genotyping PCR1 step as a csv file and import this table
plate_barcodes <- read_csv("data/genotyping/sample_NOC131/NOC131_plate_barcodes.csv")

# Create the combined plate_barcode_sequence
plate_barcodes <- plate_barcodes %>% mutate(plate_barcode_sequence = str_c(Forward_barcode_sequence, "\\\\+", Reverse_barcode_sequence))

# Import the table of single-cell Fluidigm barcodes, which are shared across every plate
fluidigm_barcodes <- read_csv("data/genotyping/Fluidigm_384_barcodes.csv")

# Join the plate barcodes to the Fluidigm barcodes to make 384 single-cell barcodes for each plate barcode
genotyping_indexes <- plate_barcodes %>% 
  right_join(fluidigm_barcodes, by = "Fluidigm_barcode_plate", relationship = "many-to-many")


## Join to FACS indexing metadata ####

# We join the barcodes to the FACS indexing data, which provides the sample donor ID for each cell.

# Import FACS indexing data
index_data_all_plates <- read_tsv("data/FACS_indexing/FACS_index_data.tsv")

cell_ids <- index_data_all_plates %>% dplyr::select(Cell, Sample, Sample_type, Sort, Plate, Index, TIME)

# Join cell IDs to the genotyping barcodes
genotyping_indexes <- genotyping_indexes %>% left_join(cell_ids)


## Assign No Template Controls ####

# On each plate, some wells were left empty as no template controls (blanks). 
# These wells have no associated FACS indexing data. We will therefore assign the NA's as blanks.

# In this experiment, wells 23O and 23P were blanks on each plate:
genotyping_indexes %>% filter(is.na(TIME)) %>% pull(Index) %>% unique()

# Assign Sample, Sample_type and Sort for the blanks
genotyping_indexes <- genotyping_indexes %>% 
  mutate(Sample = if_else(is.na(TIME), "Blank", Sample),
         Sample_type = if_else(is.na(TIME), "No_template_control", Sample_type),
         Sort = if_else(is.na(Sort), str_extract(Plate, "S\\d+?"), Sort))

# Create a cell ID for the blank wells
genotyping_indexes <- genotyping_indexes %>% 
  mutate(Cell = if_else(is.na(Cell), paste(Plate, Index, Sample, sep = "-"), Cell))


## Output metadata files ####

# Output a metadata file for the genotyping demultiplexing script
demultiplexing_barcodes <- genotyping_indexes %>% 
  select(Cell, plate_barcode_sequence, PCR2_plate, well_barcode_sequence, well_id) %>% 
  rename(plate_id = "PCR2_plate", cell_id = Cell) %>% 
  mutate(plate_id = paste("Plate", plate_id, sep = "_"))

demultiplexing_barcodes %>% write_tsv("data/genotyping/genotyping_demultiplexing.txt")


# Output an file with all the metadata
genotyping_indexes %>% 
  dplyr::select(Cell, Sample, Sample_type, Sort, Plate, PCR2_plate, Index) %>% 
  write_tsv("data/genotyping/genotyping_metadata.txt")

