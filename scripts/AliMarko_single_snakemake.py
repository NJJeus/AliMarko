import os
import re
import glob
import random


workdir: config['workdir'] 
basedir = config['output_directory'] # A folder where output files is written
base = basedir
input_folder = config['input_directory'] # A folder with input fastq files 
genome_reference = config['genome_reference'] # A fasta file with reference sequnces for alignment
HMM_folder = config['hmm_folder']# A folder with HMM for analyzis
HMM_info = config['hmm_info'] # A folder with taxonomy information of HMM
kraken_database = config['kraken_database']

suffix = '.fastq.gz' # An ending and extension of FASTQ files. They may be comressed with gz or not

files, = glob_wildcards(input_folder+"{file}"+suffix)

want_all = (expand(basedir + 'htmls/{file}.html', file=files))

rule all:
    input: basedir + 'htmls/general_html.html', want_all
    
rule index_reference:
    input: genome_reference
    output : genome_reference + '.amb'
    conda:
         'envs/bwa.yaml'
    shell:
        """
        bwa index {input}
        
        """ 
    
rule kraken2:
    input: read1 = input_folder+"{file}"+ suffix
    output: kraken2_report = basedir + "kraken_results/{file}.report.gz",
            kraken2_out = temp(basedir + "kraken_results/{file}.out"),
            read1 = temp(basedir + 'not_classified_fastq/{file}.fastq.gz')
    params: 
            kraken2_db = kraken_database,
            sample = lambda wildcards: wildcards.file
    conda:
        'envs/kraken2.yaml'
    threads: 10
    priority:
        4
    shell: 
        f"""
        kraken2 --threads {{threads}} --confidence 0.9 --db {{params.kraken2_db}} {{input.read1}} --use-names --report {{output.kraken2_report}}.tmp --output {{output.kraken2_out}} --unclassified-out {base}/not_classified_fastq/{{params.sample}}.fastq.gz.tmp 
        gzip -c {{output.read1}}.tmp > {{output.read1}}; rm {{output.read1}}.tmp
        gzip -c {{output.kraken2_report}}.tmp > {{output.kraken2_report}}; rm {{output.kraken2_report}}.tmp
        
        """
        
rule deduplicate_fastq:
    input:
        read1 = basedir + 'not_classified_fastq/{file}.fastq.gz'
    output:
        read1 = basedir + 'deduplicated_fastq/{file}.fastq.gz'
    threads: 10
    priority: 13
    conda:
        'envs/fastp.yaml'
    shell:
        """
        fastp -w {threads} -i {input.read1}  -o {output.read1} --dedup
        
        """

rule megahit:
    input:
        read1 = temp(basedir + 'not_classified_fastq/{file}.fastq.gz')
    output:
        basedir + 'megahit/{file}/contigs.fasta'
    conda:
        'envs/megahit.yaml'
    threads: 10
    params:
        sample = lambda wildcards: wildcards.file
    shell:
        f"""
        rm -r {base}/megahit/{{params.sample}}
        megahit.py  -s {{input.read1}} -o {base}/megahit/{{params.sample}} -t {{threads}}
        
        """        

rule map_extracted_fastq:
    input: 
        read1 = basedir + 'deduplicated_fastq/{file}.fastq.gz',
        reference = genome_reference,
        index_ref = genome_reference + '.amb'
    output: basedir + 'unclassified_sorted_bam/{file}.sorted.bam'
    threads: 10
    priority:
        5
    conda:
        'envs/bwa.yaml'
    shell:
        """
        bwa mem -t {threads} {input.reference} {input.read1} | samtools sort -@ {threads} -o {output}.tmp -
        python scripts/filter_bam.py -i {output}.tmp -o {output}
        samtools index {output}
        
        """
        
rule calculate_coverage_and_quality:
    input: base + 'unclassified_sorted_bam/{file}.sorted.bam'
    output: 
        coverage = base + 'calculated_coverage_and_quality/{file}_coverage.txt',
        quality = base + 'calculated_coverage_and_quality/{file}_quality.txt'
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
        coverage = basedir + 'calculated_coverage_and_quality/{file}_coverage.txt',
        bam = basedir + 'unclassified_sorted_bam/{file}.sorted.bam'
    output:
        base + 'snps/{file}.csv'
    conda:
        'envs/freebayes.yaml'
    params:
        reference = genome_reference,
        threshold = 10
    priority:
        29
    shell:
        """
        bash scripts/calculate_variance.sh -f={params.reference} -b={input.bam} -l={input.coverage} -o={output} -t={params.threshold}
        
        """
         
rule convert_coverage:
    input: 
        coverage = basedir + 'calculated_coverage_and_quality/{file}_coverage.txt',
        snps = basedir + 'snps/{file}.csv'
    output: 
        main = basedir + 'ictv_coverage/{file}.csv',
        tmp = basedir + 'tmp/cov_tmp/{file}.txt'
    conda: 'envs/scripts.yaml'
    priority:
        32
    shell:
        """
        python scripts/convert_ictv.py -c {input.coverage} -o {output.main} -t ictv_tables -s {input.snps} -m {output.tmp}
        
        """
        
rule plot_coverage:
    input:
        bam = basedir + 'unclassified_sorted_bam/{file}.sorted.bam',
        coverage = basedir + 'tmp/cov_tmp/{file}.txt'
    output: directory(basedir + 'drawings/{file}/')
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
        
rule hmm_scan:
    input: 
        read1 = basedir + 'megahit/{file}/contigs.fasta'
    params:
        reference = HMM_folder
    threads: 10
    conda: 'envs/hmm_scan.yaml'
    output:
        read1 = basedir + 'hmm_results/{file}.csv',
    shell:
        """
        python scripts/analyze_seq_hmm.py -f {input.read1} -o {output.read1} -m {params.reference} -t {threads}  --batch 10000
        
        """
            
rule hmm_report:
    input: 
        read1 = basedir + 'hmm_results/{file}.csv'
    params:
        reference = HMM_info,
    conda: 'envs/hmm_scan.yaml'
    output:
        report = basedir + 'hmm_reports/{file}.csv',
        hmm_list = basedir + 'phylo/{file}/hmm_list.csv'
    shell:
        """
        python scripts/generate_hmm_report.py -i {input.read1} -o {output.report} -m {params.reference} -t {output.hmm_list}
        
        """         

rule generalise_hmm_list:
    input: expand(basedir + 'phylo/{file}/hmm_list.csv', file=files)
    output: basedir + 'phylo/ictv_list.csv'
    shell:
        """
        cat {input} > {output}
        """        
 
rule analyse_ictv:
    input: 
        read1 = genome_reference,
        hmm_list = base + 'phylo/ictv_list.csv'
    params:
        reference = HMM_folder
    threads: 10
    priority:
        38
    conda: 'envs/hmm_scan.yaml'
    output:
        read1 = basedir + 'phylo/ictv_result.csv'
    shell:
        """
        python scripts/analyze_seq_hmm.py -f {input.read1} -o {output.read1} -m {params.reference} -t {threads}  --batch 10000 -l {input.hmm_list}
        
        """

rule ictv_report:
    input: 
        read1 = basedir + 'phylo/ictv_result.csv',
    params:
        reference = HMM_info,
    conda: 'envs/plot_hmm.yaml'
    priority:
        40
    output:
        report = base + 'phylo/ictv_report.csv'
    shell:
        """
        python scripts/generate_hmm_report.py -i {input.read1} -o {output.report} -m {params.reference} 
        
        """
        
rule collect_msa:
    input:
        genome_reference = genome_reference,
        contig_fasta = basedir + 'megahit/{file}/contigs.fasta',
        ictv_report = basedir + 'phylo/ictv_report.csv',
        contig_report = basedir + 'hmm_reports/{file}.csv'
    output:
        directory(basedir + 'phylo/msas/{file}/')
    conda: 'envs/phylo.yaml'
    shell:
        """
        python scripts/collect_msa.py -g {input.genome_reference} -c {input.contig_fasta} -i {input.ictv_report} -r {input.contig_report} -o {output}
        
        """
        
rule perform_msa:
    input: basedir + 'phylo/msas/{file}/',
    output: directory(basedir + 'phylo/msas_performed/{file}/')
    threads: 2
    conda: 'envs/phylo.yaml'
    shell:
        "bash scripts/muscle_script.sh {input} {output} {threads}"
    
rule create_tree:
    input: basedir + 'phylo/msas_performed/{file}/',
    output: directory(basedir + 'phylo/trees/{file}/')
    threads: 5
    conda: 'envs/phylo.yaml'
    shell:
        "bash scripts/create_tree.sh {input} {output} {threads}"        
 
rule draw_tree:
    input: basedir + 'phylo/trees/{file}/'
    output: directory(basedir + 'phylo/trees_drawings/{file}/')
    conda: 'envs/phylo.yaml'
    shell:
        "python scripts/draw_trees.py -i {input} -o {output}"
        
rule plot_hmm:
    input:
        report= basedir + 'hmm_reports/{file}.csv'
    output: directory(basedir + 'hmm_plots/{file}/')
    priority:
        41
    conda:
        'envs/plot_hmm.yaml'
    params:
        reference = genome_reference
    
    shell:
        """
        python scripts/plot_hmm.py -i {input.report} -o {output}
        
        """
        
rule make_html:
    input:
        coverage = basedir + 'ictv_coverage/{file}.csv',
        drawings = basedir +'drawings/{file}',
        hmm = basedir + 'hmm_reports/{file}.csv',
        hmm_plots = basedir + 'hmm_plots/{file}/',
        trees = directory(basedir + 'phylo/trees_drawings/{file}/')
    output: basedir + 'htmls/{file}.html'
    priority:
        45
    shell:
        """
        python scripts/make_html.py -c {input.coverage} -d {input.drawings} -o {output} -m {input.hmm} -a {input.hmm_plots} -t {input.trees}
        
        """           
    
rule make_general_files:
    input: 
        coverages = expand(basedir + 'ictv_coverage/{file}.csv', file=files),
        hmms = expand(basedir + 'hmm_reports/{file}.csv', file=files)
    conda: 'envs/plot_hmm.yaml'
    output: 
        hmm_general = basedir + 'hmm_general.csv',
        coverage_general = basedir + 'coverage_general.csv',
        hmm_general_pic = basedir + 'hmm_general.png',
        coverage_general_pic = basedir + 'coverage_general.png'
    params:
        reference = genome_reference,
        coverage_dir = basedir + 'ictv_coverage/',
        hmm_dir = basedir + 'hmm_reports/'
    shell:
        """
        python scripts/make_general_plots.py -i {params.coverage_dir} -m {params.hmm_dir} -p {output.coverage_general_pic} -c {output.coverage_general} -t {output.hmm_general} -u {output.hmm_general_pic}

        """
        
rule general_html:
    input: 
        hmm_general = basedir + 'hmm_general.csv',
        coverage_general = basedir + 'coverage_general.csv',
        hmm_general_pic = basedir + 'hmm_general.png',
        coverage_general_pic = basedir + 'coverage_general.png'
    output: 
        basedir + 'htmls/general_html.html'
    conda: 'envs/plot_hmm.yaml'
    shell:
        """
        python scripts/make_general_html.py -c {input.coverage_general} -e {input.coverage_general_pic} -m {input.hmm_general} -p {input.hmm_general_pic} -o {output} 
        
        """        
