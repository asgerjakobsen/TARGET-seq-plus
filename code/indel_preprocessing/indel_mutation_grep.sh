# Indel mutation grep in separated cDNA and gDNA amplicons

# In this example, demultiplexed FASTQ files are in a directory called fastq_NOC171/
# We grep for the DNMT3A p.W581fs mutation

#!/bin/bash

#### Separate cDNA and gDNA ####

# First we use cutadapt to separate cDNA and gDNA amplicons based on the PCR1 primer sequences provided in fasta files.

module load cutadapt

mkdir -p split_fastq/report

for i in $(ls fastq_NOC171/*_R1.fq.gz); \
do read1=`basename $i`; \
sample=`echo $read1 | sed 's/_R[12].fq.gz//'`; \
fastq2=`echo $i | sed 's/_R1/_R2/'`; \
echo $sample; \
cutadapt --cores=2 \
    -g file:Fwd_primers.fasta -G file:Rev_primers.fasta \
    --pair-adapters --discard-untrimmed \
    -o split_fastq/${sample}_{name}_R1.fq.gz \
    -p split_fastq/${sample}_{name}_R2.fq.gz \
    $i  $fastq2 \
    > split_fastq/report/${sample}_trimming_report.txt; \
done


#### Count WT and mutant alleles ####

# Then fastq-grep is used to count the instances of the WT and mutant sequences in each cell.

module load fastq-tools

mkdir -p mutation_grep


## gDNA Rev R2 WT

# grep for the WT sequence

for i in $(ls split_fastq/*gDNA_R2.fq.gz); \
do fastq=`basename $i`; \
sample=`echo $fastq | sed 's/_.*_R[12].fq.gz//'`; \
echo -e "$sample\t`zcat $i | fastq-grep -c AGACCCCTGGA`" >> mutation_grep/DNMT3A_W581fs_gDNA_WT.txt; \
done

## gDNA Rev R2 MUT

# grep for the MUT sequence

for i in $(ls split_fastq/*gDNA_R2.fq.gz); \
do fastq=`basename $i`; \
sample=`echo $fastq | sed 's/_.*_R[12].fq.gz//'`; \
echo -e "$sample\t`zcat $i | fastq-grep -c AGACCCCCTGGA`" >> mutation_grep/DNMT3A_W581fs_gDNA_MUT.txt; \
done



## cDNA Fwd R1 WT

# grep for the WT sequence

for i in $(ls split_fastq/*cDNA_R1.fq.gz); \
do fastq=`basename $i`; \
sample=`echo $fastq | sed 's/_.*_R[12].fq.gz//'`; \
echo -e "$sample\t`zcat $i | fastq-grep -c AGACCCCTGGA`" >> mutation_grep/DNMT3A_W581fs_cDNA_WT.txt; \
done



## mRNA Fwd R1 MUT

# grep for the MUT sequence

for i in $(ls split_fastq/*cDNA_R1.fq.gz); \
do fastq=`basename $i`; \
sample=`echo $fastq | sed 's/_.*_R[12].fq.gz//'`; \
echo -e "$sample\t`zcat $i | fastq-grep -c AGACCCCCTGGA`" >> mutation_grep/DNMT3A_W581fs_cDNA_MUT.txt; \
done


