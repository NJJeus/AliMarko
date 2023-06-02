import os
import re
import glob
import random


base = "DATA/test/"
input_folder = "DATA/test/raw_fastq/"
ictv_db_folder = 'DATA/bats/' + 'all_virus_reference/'



suffix_1 = "_1.fq"
suffix_2 = "_2.fq"

files, = glob_wildcards(input_folder+"{file}"+suffix_1)

want_all = (expand(f'{base}drawings/{{file}}/', file=files))

rule all:
    input: base + 'result_coverage_table.tsv',  want_all


    

rule map_raw_fastq:
    input: 
        read1 = input_folder+"{file}"+suffix_1,
        read2 = input_folder+"{file}"+suffix_2,
        reference = ictv_db_folder + 'all_genomes.fasta'
    output: temp(base+'bam_sorted/'+'{file}'+'.sorted.bam')
    threads: 20
    priority: 1
    conda:
        "envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} {input.reference} {input.read1} {input.read2} | samtools sort -@ {threads} -o {output} -
        
        """
    
    
rule extract_mapped_bam_sorted:  
    input: base+'bam_sorted/'+'{file}'+'.sorted.bam'
    output: temp(base+'mapped_bam_sorted/'+'{file}'+'_mapped.sorted.bam')
    threads: 20
    priority: 2
    conda:
        "envs/bwa.yaml"
    shell:
        """
        samtools view -b -F 4 -f 1 {input} --threads {threads} > {output}
        """

rule extract_mapped_fastq:
    input: base+'mapped_bam_sorted/'+'{file}'+'_mapped.sorted.bam'
    output: 
        read1=base+'mapped_fastq/'+"{file}"+'_1.fastq.gz',
        read2=base+'mapped_fastq/'+"{file}"+'_2.fastq.gz'
    threads: 20
    priority: 3
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
            kraken2_out = temp(base+"kraken_results/{file}.out"),
            read1=base+'mapped_not_classified_fastq/'+"{file}"+'_1.fastq.gz',
            read2=base+'mapped_not_classified_fastq/'+"{file}"+'_2.fastq.gz'
    params: 
            kraken2_db = "/mnt/disk1/DATABASES/kraken2/pro_and_eu",
            sample = lambda wildcards: wildcards.file
    priority: 4
    conda:
        "envs/kraken2.yaml"
    threads: 10
    shell: 
        f"""
        kraken2 --threads {{threads}} --confidence 0.9 --db {{params.kraken2_db}} {{input.read1}} {{input.read2}} --use-names --report {{output.kraken2_report}} --output {{output.kraken2_out}} --unclassified-out {base}/mapped_not_classified_fastq/{{params.sample}}#.fastq.gz.tmp --paired
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
    priority: 5
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
        reference = ictv_db_folder + 'all_genomes.fasta'
    output: base + 'unclassified_sorted_bam/' + '{file}' + '.sorted.bam'
    threads: 10
    priority: 6
    conda:
        "envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} {input.reference} {input.read1} {input.read2} | samtools sort -@ {threads} -o {output} -
        samtools index {output}
        """
        
rule calculate_coverage_and_quality:
    input: base + 'unclassified_sorted_bam/' + '{file}' + '.sorted.bam'
    output: 
        coverage=base + 'calculated_coverage_and_quality/' + '{file}_coverage' + '.txt',
        quality=base + 'calculated_coverage_and_quality/' + '{file}_quality' + '.txt'
    priority:
        6
    conda:
        "envs/bwa.yaml"
    shell:
        """
        samtools coverage {input} > {output.coverage}
        samtools view {input} | awk '{{print $3 "\t" $5}}' > {output.quality}
        sed -i '/^\*/d' {output.quality}
        """

rule convert_coverage:
    input: 
        coverage=base + 'calculated_coverage_and_quality/' + '{file}_coverage' + '.txt',
        quality=base + 'calculated_coverage_and_quality/' + '{file}_quality' + '.txt'
    output: 
        base + 'ictv_coverage/' + '{file}' + '.csv'
    priority:
        7
    shell:
        """
        python scripts/convert_ictv.py -c {input.coverage} -q {input.quality} -o {output} -t ictv_tables
        """
rule generate_tmp_coverage_files:
    input: 
        coverage=base + 'calculated_coverage_and_quality/' + '{file}_coverage' + '.txt'
    output: f'{base}tmp/cov_tmp/{{file}}.txt'
    priority:
        8
    shell:
        """
        python scripts/draw_coverage.py -c {input.coverage} -t ictv_tables -o {output}
        """
        
rule plot_coverage:
    input:
        bam=base + 'unclassified_sorted_bam/' + '{file}' + '.sorted.bam',
        coverage=f'{base}tmp/cov_tmp/{{file}}.txt'
    output: directory(f'{base}drawings/{{file}}/')
    priority:
        9
    params:
        reference=ictv_db_folder + 'all_genomes.fasta'
    shell:
        """
        bash scripts/plot_coverage.sh -f={params.reference} -b={input.bam} -l={input.coverage} -o={output}
        """


rule generalize_coverage:
    input: 
        expand(base + 'ictv_coverage/' + '{file}' + '.csv', file=files)
    output: base + 'result_coverage_table.tsv'
    shell:
        """
        python scripts/generalize_ictv.py -f {input} -o {output}
        """
        
