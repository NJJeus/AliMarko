import os
import re
import glob
import random


base = "DATA/test_not_paired/"

files, = glob_wildcards(base+'raw_fastq/'+"{file}"+'.fastq.gz')

want_all = (expand(base + 'ictv_coverage/' + '{file}' + '.tsv', file=files))

rule all:
    input: base + 'result_coverage_table.tsv'

rule map_raw_fastq:
    input: 
        read1 = base+'raw_fastq/'+"{file}"+'.fastq.gz',
        reference = base + 'all_virus_reference/' + 'all_genomes.fasta'
    output: base+'bam_sorted/'+'{file}'+'.sorted.bam'
    threads: 10
    conda:
        "envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} {input.reference} {input.read1} | samtools sort -@ {threads} -o {output} -
        """
    
    
rule extract_mapped_sam_sorted:  
    input: base+'bam_sorted/'+'{file}'+'.sorted.bam'
    output: base+'mapped_bam_sorted/'+'{file}'+'_mapped.sorted.bam'
    threads: 10
    conda:
        "envs/bwa.yaml"
    shell:
        """
        samtools view -b -F 4 {input} --threads {threads} > {output}
        """

rule extract_mapped_fastq:
    input: base+'mapped_bam_sorted/'+'{file}'+'_mapped.sorted.bam'
    output: 
        read1=base+'mapped_fastq/'+"{file}"+'.fastq.gz'
    threads: 10
    conda:
        "envs/bwa.yaml"
    shell:
        """
        samtools fastq {input} -0 {output.read1} --threads {threads} -c 8;
        """
    
    
rule kraken2:
    input: read1=base+'mapped_fastq/'+"{file}"+'.fastq.gz'
    output: kraken2_report = base+"kraken_results/{file}.report",
            kraken2_out = base+"kraken_results/{file}.out",
            read1=base+'mapped_not_classified_fastq/'+"{file}"+'.fastq.gz'
    params: 
            kraken2_db = "/mnt/disk1/DATABASES/kraken2/pro_and_eu",
            sample = lambda wildcards: wildcards.file
    conda:
        "envs/kraken2.yaml"
    threads: 10
    shell: 
        f"""
        kraken2 --threads {{threads}} --confidence 0.9 --db {{params.kraken2_db}} {{input.read1}} --use-names --report {{output.kraken2_report}} --output {{output.kraken2_out}} --unclassified-out {base}/mapped_not_classified_fastq/{{params.sample}}.fastq.gz.tmp 
        gzip -c {{output.read1}}.tmp > {{output.read1}}; rm {{output.read1}}.tmp
        """

rule map_extracted_fastq:
    input: 
        read1=base+'mapped_not_classified_fastq/'+"{file}"+'.fastq.gz',
        reference = base + 'all_virus_reference/' + 'all_genomes.fasta'
    output: base + 'unclassified_sorted_bam/' + '{file}' + '.sorted.bam'
    threads: 10
    conda:
        "envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} {input.reference} {input.read1} | samtools sort -@ {threads} -o {output} -
        """
        
rule calculate_coverage:
    input: base + 'unclassified_sorted_bam/' + '{file}' + '.sorted.bam'
    output: base + 'calculated_coverage/' + '{file}' + '.txt'
    conda:
        "envs/bwa.yaml"
    shell:
        """
        samtools coverage {input} > {output}
        """

rule convert_coverage:
    input: base + 'calculated_coverage/' + '{file}' + '.txt'
    output: base + 'ictv_coverage/' + '{file}' + '.csv'
    shell:
        """
        python scripts/convert_ictv.py -c {input} -o {output} -t ictv_tables
        """
        
rule generalize_coverage:
    input: 
        expand(base + 'ictv_coverage/' + '{file}' + '.csv', file=files)
    output: base + 'result_coverage_table.tsv'
    shell:
        """
        python scripts/generalize_ictv.py -f {input} -o {output}
        """