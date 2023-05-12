import os
import re
import glob
import random


base = "DATA/test/"

files, = glob_wildcards(base+'raw_fastq/'+"{file}"+'_R1.fastq.gz')

suffix_1 = "_R1.fastq.gz"
suffix_2 = "_R2.fastq.gz"


want_all = (expand(base + 'ictv_coverage/' + '{file}' + '.tsv', file=files))

rule all:
    input: base + 'result_coverage_table.tsv'

rule map_raw_fastq:
    input: 
        read1=base+'raw_fastq/'+"{file}"+'_R1.fastq.gz',
        read2=base+'raw_fastq/'+"{file}"+'_R2.fastq.gz',
        reference = base + 'all_virus_reference/' + 'all_genomes.fasta'
    output: base+'sam_sorted/'+'{file}'+'.sam.sorted'
    threads: 10
    conda:
        "envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} {input.reference} {input.read1} {input.read2} | samtools sort -@ {threads} -o {output} -
        """
    
    
rule extract_mapped_sam_sorted:  
    input: base+'sam_sorted/'+'{file}'+'.sam.sorted'
    output: base+'mapped_sam_sorted/'+'{file}'+'_mapped.sam.sorted'
    threads: 10
    conda:
        "envs/bwa.yaml"
    shell:
        """
        samtools view -b -F 4 -f 1 {input} --threads {threads} > {output}
        """

rule extract_mapped_fastq:
    input: base+'mapped_sam_sorted/'+'{file}'+'_mapped.sam.sorted'
    output: 
        read1=base+'mapped_fastq/'+"{file}"+'_1.fastq.gz',
        read2=base+'mapped_fastq/'+"{file}"+'_2.fastq.gz'
    threads: 10
    conda:
        "envs/bwa.yaml"
    shell:
        """
        samtools fastq {input} -1 {output.read1} -2 {output.read2} --threads {threads} -c 8;
        """
    
    
rule kraken2:
    input: read1=base+'mapped_fastq/'+"{file}"+'_1.fastq.gz',
           read2=base+'mapped_fastq/'+"{file}"+'_2.fastq.gz'
    output: kraken2_report = base+"kraken_results/{file}.report",
            kraken2_out = base+"kraken_results/{file}.out",
            read1=base+'mapped_not_classified_fastq/'+"{file}"+'_1.fastq.gz',
            read2=base+'mapped_not_classified_fastq/'+"{file}"+'_2.fastq.gz'
    params: 
            kraken2_db = "/mnt/disk1/DATABASES/kraken2/pro_and_eu",
            sample = lambda wildcards: wildcards.file
    conda:
        "envs/kraken2.yaml"
    threads: 10
    shell: 
        f"""
        kraken2/KRAKEN_DIR/kraken2 --threads {{threads}} --confidence 0.9 --db {{params.kraken2_db}} {{input.read1}} {{input.read2}} --use-names --report {{output.kraken2_report}} --output {{output.kraken2_out}} --unclassified-out {base}/mapped_not_classified_fastq/{{params.sample}}#.fastq.gz.tmp --paired
        gzip -c {{output.read1}}.tmp > {{output.read1}}; rm {{output.read1}}.tmp
        gzip -c {{output.read2}}.tmp > {{output.read2}}; rm {{output.read2}}.tmp
        """


rule pair_unclassified_fastq:
    input:
        read1=base+'mapped_not_classified_fastq/'+"{file}"+'_1.fastq.gz',
        read2=base+'mapped_not_classified_fastq/'+"{file}"+'_2.fastq.gz'
    output:
        read1=base+'mapped_not_classified_paired_fastq/'+"{file}"+'_1.fastq.gz',
        read2=base+'mapped_not_classified_paired_fastq/'+"{file}"+'_2.fastq.gz'
    threads: 10
    conda:
        "envs/seqkit.yaml"
    shell:
        f"""
        seqkit pair --force -j {{threads}} -1 {{input.read1}} -2 {{input.read2}}  -O {base}/mapped_not_classified_paired_fastq/
        """

rule map_extracted_fastq:
    input: 
        read1=base+'mapped_not_classified_paired_fastq/'+"{file}"+'_1.fastq.gz',
        read2=base+'mapped_not_classified_paired_fastq/'+"{file}"+'_2.fastq.gz',
        reference = base + 'all_virus_reference/' + 'all_genomes.fasta'
    output: base + 'unclassified_sorted_sam/' + '{file}' + '.sorted.sam'
    threads: 10
    conda:
        "envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} {input.reference} {input.read1} {input.read2} | samtools sort -@ {threads} -o {output} -
        """
        
rule calculate_coverage:
    input: base + 'unclassified_sorted_sam/' + '{file}' + '.sorted.sam'
    output: base + 'calculated_coverage/' + '{file}' + '.txt'
    conda:
        "envs/bwa.yaml"
    shell:
        """
        ~/samtools/samtools coverage {input} > {output}
        """

rule convert_coverage:
    input: base + 'calculated_coverage/' + '{file}' + '.txt'
    output: base + 'ictv_coverage/' + '{file}' + '.csv'
    shell:
        """
        python scripts/convert_ictv.py -c {input} -o {output} -t DATA/sheets/
        """
        
rule generalize_coverage:
    input: 
        expand(base + 'ictv_coverage/' + '{file}' + '.csv', file=files)
    output: base + 'result_coverage_table.tsv'
    shell:
        """
        python scripts/generalize_ictv.py -f {input} -o {output}
        """