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
# HMM Analysis
# --------------------------

rule hmm_scan:
    input:
        os.path.join(INPUT_DIR, "{sample}" + READ_SUFFIX)
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
        contigs = os.path.join(INPUT_DIR, "{sample}" + READ_SUFFIX)
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
        contigs = os.path.join(INPUT_DIR, "{sample}" + READ_SUFFIX),
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
        hmm_reports = expand(os.path.join(OUTPUT_DIR, "hmm_reports/{sample}.csv"), sample=SAMPLES)
    output:
        hmm_summary = os.path.join(OUTPUT_DIR, "hmm_general.csv"),
        hmm_plot = os.path.join(OUTPUT_DIR, "hmm_general.jpg"),
    params:
        hmm_dir = os.path.join(OUTPUT_DIR, "hmm_reports/")
    priority: 100  # Give this low priority to run last
    conda:
        'envs/plot_hmm.yaml'
    shell:
        """
        python scripts/make_general_plots.py \
            -m {params.hmm_dir} \
            -t {output.hmm_summary} \
            -u {output.hmm_plot}
        """

rule generate_summary_html:
    input:
        hmm_summary = rules.generate_summary_reports.output.hmm_summary,
        hmm_plot = rules.generate_summary_reports.output.hmm_plot,
        # Ensure all sample reports exist first
        _sample_reports = expand(os.path.join(OUTPUT_DIR, "htmls/{sample}.html"), sample=SAMPLES)
    output:
        os.path.join(OUTPUT_DIR, "htmls/general_html.html")
    conda:
        'envs/plot_hmm.yaml'
    shell:
        """
        python scripts/make_general_html.py \
            -m {input.hmm_summary} \
            -p {input.hmm_plot} \
            -o {output}
        """
