import os
from glob import glob

# --------------------------
# Configuration
# --------------------------

# Load configuration with clear variable names
WORKDIR = config["workdir"]
OUTPUT_DIR = config["output_directory"]
INPUT_DIR = config["input_directory"]
GENOME_REF = config["genome_reference"]
HMM_DIR = config["hmm_folder"]
HMM_INFO = config["hmm_info"]
KRAKEN_DB = config["kraken_database"]
READ_SUFFIX = config["suffix"]  # Single suffix for single-end reads
BLAST_DB = f"{GENOME_REF}.db"


def get_samples():
    """Get sample names from input files, excluding special files - SINGLE-END VERSION"""
    files = glob(os.path.join(INPUT_DIR, "*" + READ_SUFFIX))
    samples = [os.path.basename(f).replace(READ_SUFFIX, "") for f in files]
    # Filter out any non-sample files (like reports, etc.)
    return [s for s in samples if not s.startswith(('report', 'general', 'summary'))]

# --------------------------
# Sample Detection
# --------------------------

# Get sample names from input files (single-end version)
SAMPLES, = glob_wildcards(os.path.join(INPUT_DIR, "{sample}" + READ_SUFFIX))
# --------------------------
# Rule Definitions
# --------------------------

rule all:
    input:
        # Individual sample reports first
        expand(os.path.join(OUTPUT_DIR, "htmls/{sample}.html"), sample=SAMPLES),
        # Then general summary files
        os.path.join(OUTPUT_DIR, "htmls/general_html.html")

# --------------------------
# Reference Preparation
# --------------------------

rule index_reference:
    input:
        GENOME_REF
    output:
        GENOME_REF + ".amb"
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa index {input}"

rule create_blast_db:
    input:
        GENOME_REF
    output:
        f"{BLAST_DB}.ndb"
    conda:
        "envs/phylo.yaml"
    shell:
        """
        blast_db="{output}"
        blast_db=${{blast_db%.ndb}}
        makeblastdb -in {input} -dbtype nucl -out $blast_db
        """

# --------------------------
# Read Processing
# --------------------------

rule kraken2_classification:
    input:
        reads = os.path.join(INPUT_DIR, "{sample}" + READ_SUFFIX)
    output:
        report = os.path.join(OUTPUT_DIR, "kraken_results/{sample}.report"),
        out = os.path.join(OUTPUT_DIR, "kraken_results/{sample}.out"),
        uc_reads = os.path.join(OUTPUT_DIR, "not_classified_fastq/{sample}.fastq.gz")
    params:
        db = KRAKEN_DB
    threads: 10
    priority: 4
    conda: "envs/kraken2.yaml"
    shell:
        """
        kraken2 --threads {threads} --confidence 0.7 --db {params.db} \
            {input.reads} --use-names --report {output.report} \
            --output {output.out} --unclassified-out {output.uc_reads}.tmp
        gzip -c {output.uc_reads}.tmp > {output.uc_reads}
        rm {output.uc_reads}.tmp
        """

rule deduplicate_reads:
    input:
        reads = rules.kraken2_classification.output.uc_reads
    output:
        dedup_reads = os.path.join(OUTPUT_DIR, "deduplicated_fastq/{sample}.fastq.gz")
    params:
        quality = 15
    threads: 10
    priority: 8
    conda: "envs/fastp.yaml"
    shell:
        """
        fastp -w {threads} -i {input.reads} \
              -o {output.dedup_reads} \
              -q {params.quality} --dedup
        """

# --------------------------
# Assembly
# --------------------------

rule megahit_assembly:
    input:
        reads = rules.deduplicate_reads.output.dedup_reads
    output:
        contigs = os.path.join(OUTPUT_DIR, "megahit/{sample}/final.contigs.fa")
    threads: 20
    priority: 5
    conda: 'envs/megahit.yaml'
    shell:
        """
        rm -rf {OUTPUT_DIR}/megahit/{wildcards.sample}
        megahit -r {input.reads} \
                -o {OUTPUT_DIR}/megahit/{wildcards.sample} \
                -t {threads}
        """

# --------------------------
# Mapping and Variant Calling
# --------------------------

rule map_reads:
    input:
        reads = rules.deduplicate_reads.output.dedup_reads,
        ref = GENOME_REF,
        ref_index = rules.index_reference.output
    output:
        bam = os.path.join(OUTPUT_DIR, "unclassified_sorted_bam/{sample}.sorted.bam")
    threads: 10
    priority: 10
    conda: "envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.reads} \
            | samtools view -q 20 -h - \
            | samtools sort -@ {threads} -o {output.bam}.tmp -
        python scripts/filter_bam.py -i {output.bam}.tmp -o {output.bam}
        rm {output.bam}.tmp
        samtools index {output.bam}
        """

rule calculate_coverage:
    input:
        rules.map_reads.output.bam
    output:
        coverage = os.path.join(OUTPUT_DIR, "calculated_coverage_and_quality/{sample}_coverage.txt"),
        quality = os.path.join(OUTPUT_DIR, "calculated_coverage_and_quality/{sample}_quality.txt")
    priority: 12
    conda:
        "envs/bwa.yaml"
    shell:
        """
        samtools coverage {input} > {output.coverage}
        samtools view {input} | awk '{{print $3 "\t" $5}}' > {output.quality}
        sed -i '/^\*/d' {output.quality}
        """


rule call_variants:
    input:
        coverage = rules.calculate_coverage.output.coverage,
        bam = rules.map_reads.output.bam
    output:
        snps = os.path.join(OUTPUT_DIR, "snps/{sample}.csv")
    params:
        ref = GENOME_REF,
        threshold = 5
    priority: 14
    conda:
        "envs/freebayes.yaml"
    shell:
        """
        bash scripts/calculate_variance.sh \
            -f={params.ref} \
            -b={input.bam} \
            -l={input.coverage} \
            -o={output.snps} \
            -t={params.threshold}
        """

rule convert_coverage:
    input:
        coverage = rules.calculate_coverage.output.coverage,
        snps = rules.call_variants.output.snps
    output:
        main = os.path.join(OUTPUT_DIR, "ictv_coverage/{sample}.csv"),
        tmp = os.path.join(OUTPUT_DIR, "tmp/cov_tmp/{sample}.txt")
    conda: 
        'envs/scripts.yaml'
    priority: 16
    shell:
        """
        python scripts/convert_ictv.py \
            -c {input.coverage} \
            -o {output.main} \
            -t ictv_tables \
            -s {input.snps} \
            -m {output.tmp}
        """

# --------------------------
# HMM Analysis
# --------------------------

rule hmm_scan:
    input:
        rules.megahit_assembly.output.contigs
    output:
        hits = os.path.join(OUTPUT_DIR, "hmm_results/{sample}.csv")
    params:
        hmm_dir = HMM_DIR
    threads: 10
    priority: 38
    conda:
        'envs/hmm_scan.yaml'
    shell:
        """
        python scripts/analyze_seq_hmm.py \
            -f {input} \
            -o {output.hits} \
            -m {params.hmm_dir} \
            -t {threads} \
            --batch 10000
        """

rule generate_hmm_report:
    input:
        hits = rules.hmm_scan.output.hits,
        contigs = rules.megahit_assembly.output.contigs
    output:
        report = os.path.join(OUTPUT_DIR, "hmm_reports/{sample}.csv"),
        matched_contigs = os.path.join(OUTPUT_DIR, "hmm_reports/{sample}_matched_contigs.fasta"),
        hmm_list = os.path.join(OUTPUT_DIR, "phylo/{sample}/hmm_list.csv")
    params:
        hmm_info = HMM_INFO
    conda:
        'envs/plot_hmm.yaml'
    priority: 40
    shell:
        """
        python scripts/generate_hmm_report.py \
            -i {input.hits} \
            -c {input.contigs} \
            -o {output.report} \
            -m {params.hmm_info} \
            -t {output.hmm_list} \
            -n {output.matched_contigs}
        """

# --------------------------
# Phylogenetic Analysis
# --------------------------

rule collect_msa:
    input:
        ref = GENOME_REF,
        contigs = rules.megahit_assembly.output.contigs,
        ictv_report = "ictv_tables/HMM_reports.csv",
        hmm_report = rules.generate_hmm_report.output.report
    output:
        directory(os.path.join(OUTPUT_DIR, "phylo/msas/{sample}/"))
    conda:
        'envs/phylo.yaml'
    shell:
        """
        python scripts/collect_msa.py \
            -g {input.ref} \
            -c {input.contigs} \
            -i {input.ictv_report} \
            -r {input.hmm_report} \
            -o {output}
        """

rule perform_msa:
    input:
        rules.collect_msa.output
    output:
        directory(os.path.join(OUTPUT_DIR, "phylo/msas_performed/{sample}/"))
    threads: 2
    conda:
        'envs/phylo.yaml'
    shell:
        "bash scripts/muscle_script.sh {input} {output} {threads}"

rule build_tree:
    input:
        rules.perform_msa.output
    output:
        directory(os.path.join(OUTPUT_DIR, "phylo/trees/{sample}/"))
    threads: 5
    conda:
        'envs/phylo.yaml'
    shell:
        "bash scripts/create_tree.sh {input} {output} {threads}"

rule visualize_tree:
    input:
        rules.build_tree.output
    output:
        directory(os.path.join(OUTPUT_DIR, "phylo/trees_drawings/{sample}/"))
    conda:
        'envs/phylo.yaml'
    shell:
        "python scripts/draw_trees.py -i {input} -o {output}"
rule blast_matched_contigs:
    input:
        contigs = rules.generate_hmm_report.output.matched_contigs,
        blast_db = rules.create_blast_db.output
    output:
        os.path.join(OUTPUT_DIR, "hmm_reports/{sample}_blastn_results.tsv")
    conda:
        'envs/phylo.yaml'
    shell:
        """
        blast_db="{input.blast_db}"
        blast_db=${{blast_db%.ndb}}
        blastn -query {input.contigs} \
               -db $blast_db \
               -out {output} \
               -outfmt '6 qseqid sseqid stitle salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore'
        """

# --------------------------
# Visualization and Reporting
# --------------------------


rule plot_coverage:
    input:
        bam = rules.map_reads.output.bam,
        coverage = rules.convert_coverage.output.tmp
    output:
        directory(os.path.join(OUTPUT_DIR, "drawings/{sample}"))
    params:
        reference = GENOME_REF
    priority: 37
    conda:
        'envs/bamsnap.yaml'
    shell:
        """
        bash scripts/plot_coverage.sh \
            -f={params.reference} \
            -b={input.bam} \
            -l={input.coverage} \
            -o={output}
        """
rule plot_hmm_results:
    input:
        rules.generate_hmm_report.output.report
    output:
        directory(os.path.join(OUTPUT_DIR, "hmm_plots/{sample}"))
    params:
        reference = GENOME_REF
    priority: 41
    conda:
        'envs/plot_hmm.yaml'
    shell:
        """
        python scripts/plot_hmm.py \
            -i {input} \
            -o {output}
        """

rule generate_html_report:
    input:
        coverage = os.path.join(OUTPUT_DIR, "ictv_coverage/{sample}.csv"),
        drawings = rules.plot_coverage.output,
        hmm_report = rules.generate_hmm_report.output.report,
        hmm_plots = rules.plot_hmm_results.output,
        trees = rules.visualize_tree.output,
        blast_results = rules.blast_matched_contigs.output,
        fasta = rules.generate_hmm_report.output.matched_contigs
    output:
        os.path.join(OUTPUT_DIR, "htmls/{sample}.html")
    conda:
        'envs/plot_hmm.yaml'
    priority: 45
    shell:
        """
        python scripts/make_html.py \
            -c {input.coverage} \
            -d {input.drawings} \
            -o {output} \
            -m {input.hmm_report} \
            -a {input.hmm_plots} \
            -t {input.trees} \
            -b {input.blast_results} \
            -f {input.fasta}
        """

rule generate_summary_reports:
    input:
        # Wait for all sample reports to be generated first
        sample_reports = expand(os.path.join(OUTPUT_DIR, "htmls/{sample}.html"), sample=SAMPLES),
        # Then collect the data files
        coverages = expand(os.path.join(OUTPUT_DIR, "ictv_coverage/{sample}.csv"), sample=SAMPLES),
        hmm_reports = expand(os.path.join(OUTPUT_DIR, "hmm_reports/{sample}.csv"), sample=SAMPLES)
    output:
        hmm_summary = os.path.join(OUTPUT_DIR, "hmm_general.csv"),
        coverage_summary = os.path.join(OUTPUT_DIR, "coverage_general.csv"),
        hmm_plot = os.path.join(OUTPUT_DIR, "hmm_general.jpg"),
        coverage_plot = os.path.join(OUTPUT_DIR, "coverage_general.jpg")
    params:
        coverage_dir = os.path.join(OUTPUT_DIR, "ictv_coverage/"),
        hmm_dir = os.path.join(OUTPUT_DIR, "hmm_reports/")
    priority: 100  # Give this low priority to run last
    conda:
        'envs/plot_hmm.yaml'
    shell:
        """
        python scripts/make_general_plots.py \
            -i {params.coverage_dir} \
            -m {params.hmm_dir} \
            -p {output.coverage_plot} \
            -c {output.coverage_summary} \
            -t {output.hmm_summary} \
            -u {output.hmm_plot}
        """

rule generate_summary_html:
    input:
        hmm_summary = rules.generate_summary_reports.output.hmm_summary,
        coverage_summary = rules.generate_summary_reports.output.coverage_summary,
        hmm_plot = rules.generate_summary_reports.output.hmm_plot,
        coverage_plot = rules.generate_summary_reports.output.coverage_plot,
        # Ensure all sample reports exist first
        _sample_reports = expand(os.path.join(OUTPUT_DIR, "htmls/{sample}.html"), sample=SAMPLES)
    output:
        os.path.join(OUTPUT_DIR, "htmls/general_html.html")
    conda:
        'envs/plot_hmm.yaml'
    shell:
        """
        python scripts/make_general_html.py \
            -c {input.coverage_summary} \
            -e {input.coverage_plot} \
            -m {input.hmm_summary} \
            -p {input.hmm_plot} \
            -o {output}
        """
