import os
import re
import glob
import random



base = "DATA_test/" # A folder where output files is written
input_folder = "DATA_test/raw_fastq/" # A folder with input fastq files 
genome_reference = 'ictv_virus_reference.fa' # A fasta file with reference sequnces for alignment
HMM_folder = 'HMM_folder/' # A folder with HMM for analyzis
HMM_info = 'ictv_tables/hmm_info.csv' # A folder with taxonomy information of HMM


suffix_1 = "_1.fq.gz" # An ending and extension of FASTQ files. They may be comressed with gz or not
suffix_2 = "_2.fq.gz"


files, = glob_wildcards(input_folder+"{file}"+suffix_1)

want_all = (expand(f'{base}/htmls/{{file}}.html', file=files))

rule all:
    input: f"{base}/htmls/general_html.html",  want_all, f'{base}/phylo/ictv_report.csv'


rule index_reference:
    input: genome_reference
    output : genome_reference + '.amb'
    conda:
         "envs/bwa.yaml"
    shell:
        """
        bwa index {input}
        """

        
rule kraken2:
    input: read1 = input_folder+"{file}"+suffix_1,
            read2 = input_folder+"{file}"+suffix_2
    output: kraken2_report = base+"kraken_results/{file}.report",
            kraken2_out = base+"kraken_results/{file}.out",
            read1=base+'not_classified_fastq/'+"{file}"+'_1.fastq.gz',
            read2=base+'not_classified_fastq/'+"{file}"+'_2.fastq.gz'
    params: 
            kraken2_db = "/mnt/disk1/DATABASES/kraken2/pro_and_eu",
            sample = lambda wildcards: wildcards.file
    priority: 4
    conda:
        "envs/kraken2.yaml"
    threads: 10
    shell: 
        f"""
        kraken2 --threads {{threads}} --confidence 0.7 --db {{params.kraken2_db}} {{input.read1}} {{input.read2}} --use-names --report {{output.kraken2_report}} --output {{output.kraken2_out}} --unclassified-out {base}/not_classified_fastq/{{params.sample}}#.fastq.gz.tmp --paired
        gzip -c {{output.read1}}.tmp > {{output.read1}}; rm {{output.read1}}.tmp
        gzip -c {{output.read2}}.tmp > {{output.read2}}; rm {{output.read2}}.tmp
        """



rule spades:
    input:
        read1=base+'not_classified_fastq/'+"{file}"+'_1.fastq.gz',
        read2=base+'not_classified_fastq/'+"{file}"+'_2.fastq.gz'
    output:
        f"{base}/spades/{{file}}/contigs.fasta"
    conda:
        'envs/spades.yaml'
    threads: 20
    priority: 5
    params:
        sample = lambda wildcards: wildcards.file
    shell:
        f"""
        spades.py --meta -1 {{input.read1}} -2 {{input.read2}} -o {base}/spades/{{params.sample}} -t {{threads}} 
        """
    
    
rule deduplicate_fastq:
    input:
        read1=base+'not_classified_fastq/'+"{file}"+'_1.fastq.gz',
        read2=base+'not_classified_fastq/'+"{file}"+'_2.fastq.gz'
    output:
        read1=base+'deduplicated_fastq/'+"{file}"+'_1.fastq.gz',
        read2=base+'deduplicated_fastq/'+"{file}"+'_2.fastq.gz'
    threads: 10
    priority: 8
    params: 
        quality = 15
    conda:
        "envs/fastp.yaml"
    shell:
        """
        fastp -w {threads} -i {input.read1} -I {input.read2} -o {output.read1} -O {output.read2} -q {params.quality} --dedup
        """

rule map_extracted_fastq:
    input: 
        read1=base+'deduplicated_fastq/'+"{file}"+'_1.fastq.gz',
        read2=base+'deduplicated_fastq/'+"{file}"+'_2.fastq.gz',
        reference = genome_reference,
        index_ref = genome_reference + '.amb'
    output: base + 'unclassified_sorted_bam/' + '{file}' + '.sorted.bam'
    threads: 10
    priority: 10
    conda:
        "envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} {input.reference} {input.read1} {input.read2} | samtools view -q 20 -h - | samtools sort -@ {threads} -o {output} - 
        samtools index {output}
        """
        
rule calculate_coverage_and_quality:
    input: base + 'unclassified_sorted_bam/' + '{file}' + '.sorted.bam'
    output: 
        coverage=base + 'calculated_coverage_and_quality/' + '{file}_coverage' + '.txt',
        quality=base + 'calculated_coverage_and_quality/' + '{file}_quality' + '.txt'
    priority:12
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
        14
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
        16
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
    conda:
        'envs/bamsnap.yaml'
    params:
        reference=genome_reference
    
    shell:
        """
        bash scripts/plot_coverage.sh -f={params.reference} -b={input.bam} -l={input.coverage} -o={output}
        """
        

rule hmm_scan:
    input: 
        read1=f"{base}/spades/{{file}}/contigs.fasta"
    params:
        reference = HMM_folder,
    threads: 10
    priority:
        38
    conda: 'envs/hmm_scan.yaml'
    output:
        read1=f'{base}/hmm_results/{{file}}.csv'
    shell:
        f"""
        python scripts/analyze_seq_hmm.py -f {{input.read1}} -o {{output.read1}} -m {{params.reference}} -t {{threads}}  --batch 10000 
        """
         
    
        
        
rule hmm_report:
    input: 
        read1=f'{base}/hmm_results/{{file}}.csv'
    params:
        reference = HMM_info,
    conda: 'envs/plot_hmm.yaml'
    priority:
        40
    output:
        report=f'{base}/hmm_reports/{{file}}.csv',
        hmm_list=f'{base}/phylo/{{file}}/hmm_list.csv'
    shell:
        f"""
        python scripts/generate_hmm_report.py -i {{input.read1}} -o {{output.report}} -m {{params.reference}} -t {{output.hmm_list}}
        """
rule generalise_hmm_list:
    input: expand(f'{base}/phylo/{{file}}/hmm_list.csv', file=files)
    output: base + 'phylo/' + 'ictv_list.csv'
    shell:
        """
        cat {input} > {output}
        """
rule analyse_ictv:
    input: 
        read1=genome_reference,
        hmm_list = base + 'phylo/' + 'ictv_list.csv'
    params:
        reference = HMM_folder
    threads: 10
    priority:
        38
    conda: 'envs/hmm_scan.yaml'
    output:
        read1=f'{base}/phylo/ictv_result.csv'
    shell:
        f"""
        python scripts/analyze_seq_hmm.py -f {{input.read1}} -o {{output.read1}} -m {{params.reference}} -t {{threads}}  --batch 10000 -l {{input.hmm_list}}
        """

rule ictv_report:
    input: 
        read1=f'{base}/phylo/ictv_result.csv'
    params:
        reference = HMM_info,
    conda: 'envs/plot_hmm.yaml'
    priority:
        40
    output:
        report=f'{base}/phylo/ictv_report.csv'
    shell:
        f"""
        python scripts/generate_hmm_report.py -i {{input.read1}} -o {{output.report}} -m {{params.reference}} 
        """
        
rule collect_msa:
    input:
        genome_reference = genome_reference,
        contig_fasta = f"{base}/spades/{{file}}/contigs.fasta",
        ictv_report = f'{base}/phylo/ictv_report.csv',
        contig_report = f'{base}/hmm_reports/{{file}}.csv'
    output:
        directory(f"{base}/phylo/msas/{{file}}/")
    conda: 'envs/phylo.yaml'
    shell:
        """
        python scripts/collect_msa.py -g {input.genome_reference} -c {input.contig_fasta} -i {input.ictv_report} -r {input.contig_report} -o {output}
        """



        
rule perform_msa:
    input: directory(f"{base}/phylo/msas/{{file}}/"),
    output: directory(f"{base}/phylo/msas_performed/{{file}}/")
    threads: 2
    conda: 'envs/phylo.yaml'
    shell:
        "bash scripts/muscle_script.sh {input} {output} {threads}"
    
rule create_tree:
    input: directory(f"{base}/phylo/msas_performed/{{file}}/"),
    output: directory(f"{base}/phylo/trees/{{file}}/")
    threads: 5
    conda: 'envs/phylo.yaml'
    shell:
        "bash scripts/create_tree.sh {input} {output} {threads}"        
 
rule draw_tree:
    input: directory(f"{base}/phylo/trees/{{file}}/")
    output: directory(f"{base}/phylo/trees_drawings/{{file}}/")
    conda: 'envs/phylo.yaml'
    shell:
        "python scripts/draw_trees.py -i {input} -o {output}"

rule plot_hmm:
    input:
        report=f'{base}/hmm_reports/{{file}}.csv'
    output: directory(f'{base}/hmm_plots/{{file}}/')
    priority:
        41
    conda:
        'envs/plot_hmm.yaml'
    params:
        reference=genome_reference
    
    shell:
        """
        python scripts/plot_hmm.py -i {input.report} -o {output}
        """
        

rule make_html:
    input:
        coverage=base + 'ictv_coverage/' + '{file}' + '.csv',
        drawings=base+'drawings/'+ '{file}',
        hmm=f'{base}/hmm_reports/{{file}}.csv',
        hmm_plots = f'{base}/hmm_plots/{{file}}/',
        trees = directory(f"{base}/phylo/trees_drawings/{{file}}/")
    output: f'{base}/htmls/{{file}}.html'
    priority:
        45
    shell:
        """
        python scripts/make_html.py -c {input.coverage} -d {input.drawings} -o {output} -m {input.hmm} -a {input.hmm_plots} -t {input.trees}
        """           
    
        
rule make_general_files:
    input: 
        coverages = base + 'ictv_coverage/',
        hmms = base + 'hmm_reports/'
    conda: 'envs/plot_hmm.yaml'
    output: 
        hmm_general = f'{base}/hmm_general.csv',
        coverage_general = f'{base}/coverage_general.csv',
        hmm_general_pic = f'{base}/hmm_general.png',
        coverage_general_pic = f'{base}/coverage-general.csv'
    shell:
        """
        python scripts/make_general_plots.py -i {input.coverages} -m {input.hmms} -p {output.coverage_general_pic} -c {output.coverage_general} -t {output.hmm_general} -u {output.hmm_general_pic}

        """
        
rule general_html:
    input: 
        hmm_general = f'{base}/hmm_general.csv',
        coverage_general = f'{base}/coverage_general.csv',
        hmm_general_pic = f'{base}/hmm_general.png',
        coverage_general_pic = f'{base}/coverage-general.csv'
    output: 
        f"{base}/htmls/general_html.html"
    conda: 'envs/plot_hmm.yaml'
    shell:
        """
        python scripts/make_general_html.py -c {input.coverage_general} -e {input.coverage_general_pic} -m {input.hmm_general} -p {input.hmm_general_pic} -o {output} 
        """