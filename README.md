# TARGET-seq+

TARGET-seq+ is a modified version of TARGET-seq, a multi-omic method for simultaneous capture of single-cell genotype, RNA-seq, and surface protein expression. By incorporating elements of the Smart-seq3 chemistry, TARGET-seq+ increases the number of cells passing quality filters and the number of genes detected per cell, thus improving the detection of lowly expressed genes, whilst retaining high-fidelity single-cell genotyping. 

For details of the method see our manuscript:

[Jakobsen, Turkalj, Zeng et al,. Selective advantage of mutant stem cells in human clonal hematopoiesis is associated with attenuated response to inflammation and aging, Cell Stem Cell 2024](https://doi.org/10.1016/j.stem.2024.05.010)

This repository contains a collection of scripts for the analysis of TARGET-seq+ data.

## Analysis overview

Broadly, the analysis workflow follows the following steps:
1.	Analysis of FACS indexing data
2.	Analysis of single-cell genotyping data within each patient sample and integration with flow cytometry indexing data to generate metadata for the full dataset
3.	Pre-processing of transcriptome data
4.	Integration of the full dataset and downstream analysis

We provide an example dataset for testing the analysis workflow in the `data` directory.

## 1. Analysis of FACS indexing data

This step gathers cell surface immunophenotyping data from FACS index sorting files into a single table. 
From these data, we know which cell type was sorted into each well and the sample donor ID. This is important for defining which cells came from control and test samples when performing downstream analysis.

Using the cell surface immunofluorescence values, we can gate for positive and negative cells to define immunophenotypic populations.

The output of this step is a table of cell IDs and FACS indexing values for the full dataset.

## 2. Analysis of single-cell genotyping

These steps take targeted amplicon sequencing FASTQ files and assign a genotype to each cell at each locus of interest.

1. Run `01_Make_genotyping_demultiplexing_file.R` to generate a demultiplexing file for use in the TARGET-seq SCpipeline and a metadata file for genotyping analysis.

2. Run the TARGET-seq SCpipeline (https://github.com/albarmeira/TARGET-seq) to demultiplex and map the single-cell genotyping data and generate tables of allelic counts for each mutation locus per cell. 

3. Perform genotyping calling for each mutation locus: `02_Genotype_calling_NOC131_DNMT3A_I780T.R` and `03_Genotype_calling_NOC131_DNMT3A_Q606X.R`

4. Integrate the genotypes within single cells to assign clonal identities within each sample: `04_Clone_assignment_NOC131.R`

5. Integrate the FACS indexing and single-cell genotyping data to make a metadata file for multiomic analysis: `05_Integrate_metadata.R`

## 3. Pre-processing of transcriptome data

To pre-process TARGET-seq+ transcriptome data, we have written a custom python pipeline: [TARGET-seq-plus-RNA](https://github.com/asgerjakobsen/TARGET-seq-plus-RNA). 

This takes FASTQ files and performs single-cell demultiplexing, mapping, and gene feature counting. The gene counts from STARsolo are used in downstream analysis.

## 4. Integration of transcriptome data with genotyping and FACS indexing

Integration of transcriptome data with genotyping and FACS indexing data relies on a shared cell identifiers. 

The transcriptome cell barcodes in the outputs from the pre-processing pipeline are derived from the combination of libraries with unique i5/i7 plate indexes and oligo-dT cell barcodes.

The `Reformat_counts_matrix.R` script joins the transcriptome cell barcodes to the descriptive cell IDs from FACS that match those in the other metadata. It then takes the outputs from STARsolo and reformats the transcriptome counts matrix with gene name and cell IDs. 

The resulting counts matrix can be used with the integrated metadata file in single-cell genomics analysis packages for downstream analysis.

## Contact

All scripts were written by [Asger Jakobsen](https://www.imm.ox.ac.uk/people/asger-jakobsen). 

In case of problems or doubts about the scripts please raise an issue.

