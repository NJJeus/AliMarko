
# AliMarko (alignment & markov models)
v0.5  
A Pipeline for Virus Sequence Detection  
Requirements:
* snakemake
* mamba
* a kraken database

AliMarko is a bioinformatics pipiline which recognises viral sequences in metagenomic data.

It uses two methods: alignment to a database of viral genomes and performing homology analysis with HMM.

AliMarko creates an easy-to-read HTML report

## Commmon scheme
1. The pipeline gets FASTQ files
2. Filtration of the FASTQ files with Kraken2 and deleting all reads that are recognised as belonging to cellular organisms
3. Depuplication with fastp
4. Mapping the deduplicated fastq to a database of viral genomes.
5. Calculating coverage and snp for all viral genomes
6. Drawing mapping visualizations for significatly covered genomes
7. Assemble reads from the deduplicated FASTQ with SPADes, getting contigs.
8. Analyse the contigs with HMM using our script
9. Drawing visualisations for all contigs
10. Create an HTML report for all the results

![Common scheme](Documentation/common_scheme.png)

## Screenshot of a multisample report
![Pic_2](Documentation/PIC_2.png)

## Screenshot of a sample report
![Pic_3](Documentation/PIC_3.png)






