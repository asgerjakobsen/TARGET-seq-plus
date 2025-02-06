# TARGET-seq-plus


## Overview

Broadly, the analysis workflow follows the following steps:
1.	Analysis of FACS indexing data
2.	Analysis of single-cell genotyping data within each patient sample
3.	Integration of genotyping and flow cytometry indexing data to generate metadata for the full dataset
4.	Pre-processing of RNA-seq data
5.	Integration of the full dataset and downstream analysis


## 1. Analysis of FACS indexing data

This step gathers cell surface immunophenotyping data from FACS index sorting files into a single table. 
From these data, we know which cell type was sorted into each well and the sample donor ID. This is important for defining which cells came from control and test samples when performing downstream analysis.

Using the cell surface immunofluorescence values, we can gate for positive and negative cells to define immunophenotypic populations.

The output of this step is a table of cell IDs and FACS indexing values for the full dataset.

## 2. Analysis of single-cell genotyping

Preprocessing of single-cell genotyping data is performed using the pipeline from the original TARGET-seq manuscript: https://github.com/albarmeira/TARGET-seq

Steps:

1. Run `01_make_genotyping_demultiplexing_file.R` to generate a demultiplexing file for use in the TARGET-seq SCpipeline and a metadata file for genotyping analysis.

2. Run the TARGET-seq SCpipeline to map the single-cell genotyping data and generate tables of allelic counts for each mutation locus. 

3. Perform genotyping calling for each mutation locus: `02_NOC131_DNMT3A_I780T.R` and `03_NOC131_DNMT3A_Q606X.R`

4. Integrate the genotypes within single cells to assign clonal identities within each sample: `04_clone_assignment_NOC131.R`

5. Integrate the FACS indexing and sngle-cell genotyping data to make a metadata file for multiomic analysis: `05_make_metadata.R`

## 3. Analysis of RNA-seq


