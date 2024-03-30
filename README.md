
# AliMarko (alignment & markov models)
v0.5  
A Pipeline for Virus Sequence Detection  
Requirements:
* snakemake
* mamba
* a kraken database


AliMarko is designed to automate essential analysis steps for viral sequence data processing. It combines read alignment to viral reference genomes and Hidden Markov Models (HMM) with futher phylogenetic analysis. 

AliMarko recieve FATQ files and creates an easy-to-read HTML reports. 

A sample HTML report contains tabular and graphical representations of alignment to all sequences from a reference database and HMM analysis (see Fig. 2)

A multisample HTML report contains tabular results and heatmaps for alignment to a reference database and HMM analyzis


The pipeline can process several FASTQ files in parallel due using of Snakemake. 

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

## Fig. 2. Screenshot of a sample report
![Pic_3](Documentation/PIC_3.png)
*Fig. 2. Screenshots of the one-sample HTML report of AliMarko. A - a screenshot of the html report with visualization of alignment to a virus reference genome (simulated). Information about the reference and general information about alignment is shown in a table. If the reference genome contains several fragments. B - visualization of probably concatenation-caused alignment to a genome of Murine leukemia virus. While one fragment of the reference had good coverage depth, other parts of the genome aren’t covered at all. C- visualization of matches of HMM against a contig. The matches are colored by their score. Several models matched the contig. D - phylogenetic tree of contig of presumable Caliciviridae origin.(see C). Sequences in the tree are colored with their taxonomic group.*

## Fig. 3. Screenshot of a multisample report
![Pic_2](Documentation/PIC_2.png)








