import os
import re
import glob
import random


base = "DATA/bats/"

files, = glob_wildcards(base+'raw_fastq/'+"{file}"+'_R1.fastq.gz')
print(files)

suffix_1 = "_R1.fastq.gz"
suffix_2 = "_R2.fastq.gz"


want_all = (expand(base + 'ictv_coverage/' + '{file}' + '.tsv', file=files))
print(want_all[:1])

rule all:
    input: base + 'result_coverage_table.tsv'

rule map_raw_fastq:
    input: 
        read1=base+'raw_fastq/'+"{file}"+'_R1.fastq.gz',
        read2=base+'raw_fastq/'+"{file}"+'_R2.fastq.gz',
        reference = base + 'all_virus_reference/' + 'all_genomes.fasta'
    output: base+'sam_sorted/'+'{file}'+'.sam.sorted'
    threads: 10
    shell:
        """
        bwa mem -t {threads} {input.reference} {input.read1} {input.read2} | samtools sort -@ {threads} -o {output} -
        """
    
    
rule extract_mapped_sam_sorted:  
    input: base+'sam_sorted/'+'{file}'+'.sam.sorted'
    output: base+'mapped_sam_sorted/'+'{file}'+'_mapped.sam.sorted'
    threads: 10
    shell:
        """
        samtools view -b -F 4 -f 1 {input} --threads {threads} > {output}
        """

rule extract_mapped_fastq:
    input: base+'mapped_sam_sorted/'+'{file}'+'_mapped.sam.sorted'
    output: 
        read1=base+'mapped_fastq/'+"{file}"+'_1.fastq',
        read2=base+'mapped_fastq/'+"{file}"+'_2.fastq'
    threads: 10
    shell:
        """
        samtools fastq {input} -1 {output.read1} -2 {output.read2} --threads {threads};
        """
    
    
rule kraken2:
    input: read1=base+'mapped_fastq/'+"{file}"+'_1.fastq',
           read2=base+'mapped_fastq/'+"{file}"+'_2.fastq'
    output: kraken2_report = base+"kraken_results/{file}.report",
            kraken2_out = base+"kraken_results/{file}.out"
    params: kraken2_db = "/mnt/disk1/DATABASES/kraken2/pro_and_eu"
    threads: 10
    shell: 
        """
        kraken2/KRAKEN_DIR/kraken2 --threads {threads} --confidence 0.9 --db {params.kraken2_db} {input.read1} {input.read2} --use-names --report {output.kraken2_report} --output {output.kraken2_out}
        """
        
rule extract_seq_ids:
    input: base+"kraken_results/{file}.out"
    output: base + 'unclassified_ids/' + '{file}' + '.txt'
    shell:
        """
        grep 'taxid 0' {input} --no-filename | cut -f2 > {output}
        """

rule extract_unclassified_fastq:
    input: 
        ids = base + 'unclassified_ids/' + '{file}' + '.txt',
        read1 = base+'mapped_fastq/'+"{file}"+'_1.fastq',
        read2 = base+'mapped_fastq/'+"{file}"+'_2.fastq'    
    output: 
        read1=base+'mapped_not_classified_fastq/'+"{file}"+'_1.fastq.gz',
        read2=base+'mapped_not_classified_fastq/'+"{file}"+'_2.fastq.gz'
    threads: 10
    shell:
        """
        seqkit grep -j {threads} -f {input.ids} {input.read1} | gzip > {output.read1}
        seqkit grep -j {threads} -f {input.ids} {input.read2} | gzip > {output.read2}
        """
rule pair_unclassified_fastq:
    input:
        read1=base+'mapped_not_classified_fastq/'+"{file}"+'_1.fastq.gz',
        read2=base+'mapped_not_classified_fastq/'+"{file}"+'_2.fastq.gz'
    output:
        read1=base+'mapped_not_classified_paired_fastq/'+"{file}"+'_1.fastq.gz',
        read2=base+'mapped_not_classified_paired_fastq/'+"{file}"+'_2.fastq.gz'
    threads: 10
    shell:
        """
        seqkit pair --force -j {threads} -1 {input.read1} -2 {input.read2}  -O DATA/bats/mapped_not_classified_paired_fastq/
        """

rule map_extracted_fastq:
    input: 
        read1=base+'mapped_not_classified_paired_fastq/'+"{file}"+'_1.fastq.gz',
        read2=base+'mappedpp_not_classified_paired_fastq/'+"{file}"+'_2.fastq.gz',
        reference = base + 'all_virus_git remote add origin https://github.com/NJJeus/surprise_ictv_pipeline.gitreference/' + 'all_genomes.fasta'
    output: base + 'unclassified_sorted_sam/' + '{file}' + '.sorted.sam'
    threads: 10
    shell:
        """
        bwa mem -t {threads} {input.reference} {input.read1} {input.read2} | samtools sort -@ {threads} -o {output} -
        """
        
rule calculate_coverage:
    input: base + 'unclassified_sorted_sam/' + '{file}' + '.sorted.sam'
    output: base + 'calculated_coverage/' + '{file}' + '.txt'
    shell:
        """
        ~/samtools/samtools coverage {input} > {output}
        """

rule convert_coverage:
    input: 
        expand(base + 'ictv_coverage/' + '{file}' + '.csv', file=files)
    output: base + 'result_coverage_table.tsv'
    shell:
        """
        python scripts/generalize_ictv.py -f {input} -o {output}
        """