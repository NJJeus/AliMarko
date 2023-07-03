import os
import re
import glob
import random


base = "DATA_test/"
input_folder = "DATA_test/raw_fastq/"
genome_reference = 'ictv_virus_reference.fa'

suffix = '.fastq.gz'

files, = glob_wildcards(input_folder+"{file}"+ suffix)

want_all = (expand(f'{base}/htmls/{{file}}.html', file=files))

print('Files:', end=' ')
print(files)

rule all:
    input: base + 'result_coverage_table.tsv',  want_all
    
rule index_reference:
    input: genome_reference
    output : genome_reference + '.amb'
    conda:
         "envs/bwa.yaml"
    shell:
        """
        bwa index {input}
        """

rule map_raw_fastq:
    input: 
        read1 = input_folder + "{file}"+'.fastq.gz',
        reference = genome_reference,
        reference_index = genome_reference + '.amb'
    output: temp(base+'bam_sorted/'+'{file}'+'.sorted.bam')
    threads: 10
    conda:
        "envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} {input.reference} {input.read1} | samtools sort -@ {threads} -o {output} -
        """
    
    
rule extract_mapped_sam_sorted:  
    input: base+'bam_sorted/'+'{file}'+'.sorted.bam'
    output: temp(base+'mapped_bam_sorted/'+'{file}'+'_mapped.sorted.bam')
    threads: 10
    priority:
        2
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
    priority:
        3
    conda:
        "envs/bwa.yaml"
    shell:
        """
        samtools fastq {input} -0 {output.read1} --threads {threads} -c 8;
        """
    
    
rule kraken2:
    input: read1=base+'mapped_fastq/'+"{file}"+'.fastq.gz'
    output: kraken2_report = base+"kraken_results/{file}.report.gz",
            kraken2_out = temp(base+"kraken_results/{file}.out"),
            read1=temp(base+'mapped_not_classified_fastq/'+"{file}"+'.fastq.gz')
    params: 
            kraken2_db = "/mnt/disk1/DATABASES/kraken2/pro_and_eu",
            sample = lambda wildcards: wildcards.file
    conda:
        "envs/kraken2.yaml"
    threads: 10
    priority:
        4
    shell: 
        f"""
        kraken2 --threads {{threads}} --confidence 0.9 --db {{params.kraken2_db}} {{input.read1}} --use-names --report {{output.kraken2_report}}.tmp --output {{output.kraken2_out}} --unclassified-out {base}/mapped_not_classified_fastq/{{params.sample}}.fastq.gz.tmp 
        gzip -c {{output.read1}}.tmp > {{output.read1}}; rm {{output.read1}}.tmp
        gzip -c {{output.kraken2_report}}.tmp > {{output.kraken2_report}}; rm {{output.kraken2_report}}.tmp
        """
        
rule deduplicate_fastq:
    input:
        read1=base+'mapped_not_classified_fastq/'+"{file}"+'.fastq.gz'
    output:
        read1=base+'deduplicated_fastq/'+"{file}"+'.fastq.gz'
    threads: 10
    priority: 13
    conda:
        "envs/fastp.yaml"
    shell:
        """
        fastp -w {threads} -i {input.read1}  -o {output.read1}
        """

rule map_extracted_fastq:
    input: 
        read1=base+'deduplicated_fastq/'+"{file}"+'.fastq.gz',
        reference = genome_reference
    output: base + 'unclassified_sorted_bam/' + '{file}' + '.sorted.bam'
    threads: 10
    priority:
        5
    conda:
        "envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} {input.reference} {input.read1} | samtools sort -@ {threads} -o {output} -
        samtools index {output}
        """
        
rule calculate_coverage_and_quality:
    input: base + 'unclassified_sorted_bam/' + '{file}' + '.sorted.bam'
    output: 
        coverage=base + 'calculated_coverage_and_quality/' + '{file}_coverage' + '.txt',
        quality=base + 'calculated_coverage_and_quality/' + '{file}_quality' + '.txt'
    priority:
        24
    conda:
        "envs/bwa.yaml"
    shell:
        """
        samtools coverage {input} > {output.coverage}
        samtools view {input} | awk '{{print $3 "\t" $5}}' > {output.quality}
        sed -i '/^\*/d' {output.quality}
        """

rule count_snp:
    input: 
        coverage = base + 'calculated_coverage_and_quality/' + '{file}_coverage' + '.txt',
        bam= base + 'unclassified_sorted_bam/' + '{file}' + '.sorted.bam'
    output:
        base + '/snps/{file}.csv'
    conda:
        "envs/freebayes.yaml"
    params:
        reference = genome_reference,
        threshold=10
    priority:
        29
    shell:
        """
        bash scripts/calculate_variance.sh -f={params.reference} -b={input.bam} -l={input.coverage} -o={output} -t={params.threshold}
        """
    
        
rule convert_coverage:
    input: 
        coverage=base + 'calculated_coverage_and_quality/' + '{file}_coverage' + '.txt',
        snps= base + '/snps/{file}.csv'
    output: 
        main=base + 'ictv_coverage/' + '{file}' + '.csv',
        tmp=f'{base}tmp/cov_tmp/{{file}}.txt'
    conda: 'envs/scripts.yaml'
    priority:
        32
    shell:
        """
        python scripts/convert_ictv.py -c {input.coverage} -o {output.main} -t ictv_tables -s {input.snps} -m {output.tmp}
        """
        
rule plot_coverage:
    input:
        bam=base + 'unclassified_sorted_bam/' + '{file}' + '.sorted.bam',
        coverage=f'{base}tmp/cov_tmp/{{file}}.txt'
    output: directory(f'{base}drawings/{{file}}/')
    priority:
        37
    params:
        reference=genome_reference
    conda:
        'envs/bamsnap.yaml'
    shell:
        """
        bash scripts/plot_coverage.sh -f={params.reference} -b={input.bam} -l={input.coverage} -o={output}
        """
        
rule make_html:
    input:
        coverage=base + 'ictv_coverage/' + '{file}' + '.csv',
        drawings=base+'drawings/'+ '{file}'
    output: f'{base}/htmls/{{file}}.html'
    priority:
        39
    shell:
        """
        python scripts/make_html.py -c {input.coverage} -d {input.drawings} -o {output}
        """


rule generalize_coverage:
    input: 
        expand(base + 'ictv_coverage/' + '{file}' + '.csv', file=files)
    output: base + 'result_coverage_table.tsv'
    shell:
        """
        python scripts/generalize_ictv.py -f {input} -o {output}
        """
        
