configfile: "config.yaml"

import os

# Define directories based on config
OUTPUT_DIR    = config["output_dir"]
DATA_DIR      = os.path.join(OUTPUT_DIR, "data")
QUAST_DIR     = os.path.join(OUTPUT_DIR, "quast")
ANALYSIS_DIR  = os.path.join(OUTPUT_DIR, "analysis")
PAINTING_DIR  = os.path.join(ANALYSIS_DIR, "painting_probes")
ONT_MAPPING_DIR  = os.path.join(ANALYSIS_DIR, "mapping_ONT_reads")
CLIPPING_DIR  = os.path.join(ANALYSIS_DIR, "clipping_info")
CONTIG_PAIRS_DIR = os.path.join(ANALYSIS_DIR, "contig_pairs")
TIDECLUSTER_DIR = os.path.join(ANALYSIS_DIR, "tidecluster")
MITO_DB = config["mitochondrial_db"]
PLASTID_DB = config["plastid_db"]

# get limits for size fractions, if not provided, use defaults 80k and 120k
min_limit = config.get("min_limit", 80000)
max_limit = config.get("max_limit", 120000)
mid_limit = max_limit - 1
trim_start = config.get("trim_start", 0)
trim_end = config.get("trim_end", 0)


if config.get("assembly_fai_ref", False):
    make_doplots = 'true'
    FAI_REF = config["assembly_fai_ref"]
    BLAST_REF = config["blast_probes_ref"]
else:
    make_doplots = 'false'



# Create required output directories if they don't exist
def create_dirs(*dirs):
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

create_dirs(OUTPUT_DIR, QUAST_DIR, ANALYSIS_DIR, DATA_DIR,
            PAINTING_DIR, ONT_MAPPING_DIR, CLIPPING_DIR, CONTIG_PAIRS_DIR, TIDECLUSTER_DIR)
print("Output directory structure created:", os.listdir(OUTPUT_DIR))

# Final targets: final assembly fasta and a Quast marker file
rule all:
    input:
        os.path.join(OUTPUT_DIR, "hifiasm_assembly.bp.p_ctg.gfa.fasta"),
        os.path.join(QUAST_DIR, "quast_done.txt"),
        os.path.join(PAINTING_DIR, "hifiasm_assembly.bp.p_ctg_x_oligos_CAMv2r2.blast_out"),
        os.path.join(ONT_MAPPING_DIR, "hifiasm_assembly.bp.p_ctg.gfa_mapped_longest_fraction.sorted.bam"),
        os.path.join(ONT_MAPPING_DIR, "hifiasm_assembly.bp.p_ctg.gfa_mapped_midsized_fraction.sorted.bam"),
        os.path.join(CLIPPING_DIR, "hifiasm_assembly.bp.p_ctg.gfa_mapped_longest_fraction_clip_info.tsv"),
        os.path.join(CONTIG_PAIRS_DIR, "contig_pairs.csv"),
        os.path.join(TIDECLUSTER_DIR, "tc_done.txt"),
        os.path.join(CLIPPING_DIR, "hifiasm_assembly.bp.p_ctg.gfa_mapped_longest_fraction_zero_coverage.bed")


####################################################################
# 1. Convert Oxford Nanopore FASTQ to FASTA using seqkit
####################################################################

rule convert_fastq_to_fasta:
    input:
        fastq = config["oxford_nanopore_reads"]
    output:
        fasta = os.path.join(DATA_DIR, "reads.fasta")
    conda:
        "envs/pysam.yaml"
    threads: 1
    shell:
        """
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        seqkit fq2fa {input.fastq} > {output.fasta}
        """

####################################################################
# 1a. Create different read length fractions (not used for assembly)
####################################################################

rule fraction_reads:
    input:
        fastq = os.path.join(DATA_DIR, "reads.filtered.fastq"),
        fasta = os.path.join(DATA_DIR, "reads.filtered.fasta")
    output:
        fasta_long   = os.path.join(DATA_DIR, "reads_longest_fraction.fasta"),
        fastq_long   = os.path.join(DATA_DIR, "reads_longest_fraction.fastq"),
        fasta_med    = os.path.join(DATA_DIR, "reads_midsized_fraction.fasta"),
        fastq_med    = os.path.join(DATA_DIR, "reads_midsized_fraction.fastq")
    conda:
        "envs/pysam.yaml"
    threads: 4
    shell:
        """
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        seqkit seq -m {max_limit} {input.fasta} > {output.fasta_long} &
        seqkit seq -m {max_limit} {input.fastq} > {output.fastq_long} &
        seqkit seq -m {min_limit} -M {mid_limit} {input.fasta} > {output.fasta_med} &
        seqkit seq -m {min_limit} -M {mid_limit} {input.fastq} > {output.fastq_med}
        wait
        """

####################################################################
# 2. Detect contamination (plastid and mitochondrial) using BLAST
####################################################################

rule detect_contamination_plastid:
    input:
        fasta = os.path.join(DATA_DIR, "reads.fasta")
    output:
        tsv = os.path.join(DATA_DIR, "contam_plastid.tsv")
    conda:
        "envs/pysam.yaml"
    threads: workflow.cores
    shell:
        """
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        # copy fasta file to the data directory
        cp {PLASTID_DB} {DATA_DIR}/plastid.fasta
        makeblastdb -in {DATA_DIR}/plastid.fasta -dbtype nucl
        detect_contamination.py -i {input.fasta} -o {output.tsv} --min_coverage 90 --min_identity 90 --num_cpu {threads} -d {DATA_DIR}/plastid.fasta --word_size 21 --raw_blast {DATA_DIR}/plastid.blast
        """

rule detect_contamination_mito:
    input:
        fasta = os.path.join(DATA_DIR, "reads.fasta")
    output:
        tsv = os.path.join(DATA_DIR, "contam_mito.tsv")
    conda:
        "envs/pysam.yaml"
    threads: workflow.cores
    shell:
        """
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        # copy fasta file to the data directory
        cp {MITO_DB} {DATA_DIR}/mitochondrial.fasta
        makeblastdb -in {DATA_DIR}/mitochondrial.fasta -dbtype nucl
        detect_contamination.py -i {input.fasta} -o {output.tsv} --min_coverage 90 --min_identity 90 --num_cpu {threads} -d {DATA_DIR}/mitochondrial.fasta --word_size 21 --raw_blast {DATA_DIR}/mitochondrial.blast
        """

####################################################################
# 3. Remove contaminant reads from the FASTQ using filter_fastq_fasta.py
####################################################################

rule filter_fastq:
    input:
        fastq   = config["oxford_nanopore_reads"],
        fasta   = os.path.join(DATA_DIR, "reads.fasta"),
        plastid = os.path.join(DATA_DIR, "contam_plastid.tsv"),
        mito    = os.path.join(DATA_DIR, "contam_mito.tsv")
    output:
        fastq   = os.path.join(DATA_DIR, "reads.filtered.fastq"),
        fasta   = os.path.join(DATA_DIR, "reads.filtered.fasta")
    conda:
        "envs/pysam.yaml"
    threads: 1
    shell:
        """
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        filter_fastq_fasta.py -i {input.fastq} -r {input.plastid} {input.mito} -o {output.fastq} \
        --format fastq --trim_start {trim_start} --trim_end {trim_end} 
        filter_fastq_fasta.py -i {input.fasta} -r {input.plastid} {input.mito} -o {output.fasta} \
         --format fasta --trim_start {trim_start} --trim_end {trim_end}
        """

####################################################################
# 4. Run genome assembly using hifiasm (outputs in hifiasm_outputs folder)
####################################################################

rule hifiasm_assembly:
    input:
        fastq = os.path.join(DATA_DIR, "reads.filtered.fastq")
    output:
        gfa = os.path.join(OUTPUT_DIR, "hifiasm_assembly.bp.p_ctg.gfa")
    threads: workflow.cores
    shell:
        """
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        hifiasm -o {OUTPUT_DIR}/hifiasm_assembly -t {threads} --ont {input.fastq}
        """

####################################################################
# 5. Convert the hifiasm .gfa file to FASTA using awk
####################################################################

rule gfa_to_fasta:
    input:
        gfa = os.path.join(OUTPUT_DIR, "hifiasm_assembly.bp.p_ctg.gfa")
    output:
        fasta = os.path.join(OUTPUT_DIR, "hifiasm_assembly.bp.p_ctg.gfa.fasta")
    threads: 1
    shell:
        """
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        awk '/^S/{{print ">"$2"\\n"$3}}' {input.gfa} | fold > {output.fasta}
        """


####################################################################
# 6. Run Quast for quality assessment (outputs in quast folder)
####################################################################

rule run_quast:
    input:
        fasta = os.path.join(OUTPUT_DIR, "hifiasm_assembly.bp.p_ctg.gfa.fasta")
    output:
        marker = os.path.join(QUAST_DIR, "quast_done.txt")
    params:
        quast_dir = QUAST_DIR
    conda:
        "envs/quast.yaml"
    threads: 1
    shell:
        """
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        quast {input.fasta} -o {params.quast_dir}
        touch {output.marker}
        """

rule make_index_and_blast_database:
    input:
        fasta = os.path.join(OUTPUT_DIR, "hifiasm_assembly.bp.p_ctg.gfa.fasta")
    output:
        index = os.path.join(OUTPUT_DIR, "hifiasm_assembly.bp.p_ctg.gfa.fasta.fai"),
        db = os.path.join(OUTPUT_DIR, "hifiasm_assembly.bp.p_ctg.gfa.fasta.nhr")
    conda:
        "envs/basic_tools.yaml"
    threads: 1
    shell:
        """
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        samtools faidx {input.fasta}
        makeblastdb -in {input.fasta} -dbtype nucl
        """


rule map_probes:
    input:
        assembly_query=os.path.join(OUTPUT_DIR,"hifiasm_assembly.bp.p_ctg.gfa.fasta"),
        assembly_fai=os.path.join(OUTPUT_DIR,"hifiasm_assembly.bp.p_ctg.gfa.fasta.fai"),
        probes=config["painting_probes"]
    output:
        blast_probes_query=os.path.join(PAINTING_DIR,"hifiasm_assembly.bp.p_ctg_x_oligos_CAMv2r2.blast_out"),
    params:
        dotplot=os.path.join(PAINTING_DIR,"hifiasm_assembly.bp.p_ctg_x_oligos_CAMv2r2_vs_ref.png")
    threads: workflow.cores

    conda:
        "envs/basic_tools.yaml"
    shell:
        """
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
    
        # Run blastn: map the painting probes against the assembly
        blastn -task blastn -db {input.assembly_query} -query {input.probes} -evalue 1e-10 \
          -num_alignments 1 -outfmt 6 -num_threads {threads} -out {output.blast_probes_query}
        touch {output.blast_probes_query}
        
        # run next step only if make_doplots is true
        if [ "{make_doplots}" = "true" ]; then
            make_dotplot_from_probes2.R -q {output.blast_probes_query} -Q {input.assembly_fai} \
              -S {FAI_REF} -s {BLAST_REF} -o {params.dotplot}
        fi
        """


rule map_nanopore_reads_to_assembly:
    input:
        assembly=os.path.join(OUTPUT_DIR,"hifiasm_assembly.bp.p_ctg.gfa.fasta"),
        fastq_long=os.path.join(DATA_DIR,"reads_longest_fraction.fastq"),
        filtq_med=os.path.join(DATA_DIR,"reads_midsized_fraction.fastq")
    output:
        bam_long=os.path.join(ONT_MAPPING_DIR, "hifiasm_assembly.bp.p_ctg.gfa_mapped_longest_fraction.sorted.bam"),
        bam_med=os.path.join(ONT_MAPPING_DIR, "hifiasm_assembly.bp.p_ctg.gfa_mapped_midsized_fraction.sorted.bam")
    conda:
        "envs/basic_tools.yaml"
    threads: workflow.cores
    log: os.path.join(OUTPUT_DIR, "logs/map_nanopore_reads_to_assembly.log")
    shell:
        """
        # dorado is in path
        echo "lgging test"
        dorado aligner {input.assembly} {input.fastq_long}  2> {log} 1> {output.bam_long}.unsorted.bam
        samtools sort -@ {threads} -o {output.bam_long} {output.bam_long}.unsorted.bam
        samtools index -c {output.bam_long}
        rm {output.bam_long}.unsorted.bam
        
        dorado aligner {input.assembly} {input.filtq_med}  2> {log} 1> {output.bam_med}.unsorted.bam
        samtools sort -@ {threads} -o {output.bam_med} {output.bam_med}.unsorted.bam
        samtools index -c {output.bam_med}
        rm {output.bam_med}.unsorted.bam
        
        """


rule detect_clipping:
    input:
        bam_long=os.path.join(ONT_MAPPING_DIR, "hifiasm_assembly.bp.p_ctg.gfa_mapped_longest_fraction.sorted.bam"),
    output:
        clipping_long=os.path.join(CLIPPING_DIR, "hifiasm_assembly.bp.p_ctg.gfa_mapped_longest_fraction_clip_info.tsv"),
    conda:
        "envs/pysam.yaml"
    params:
         prefix = os.path.join(CLIPPING_DIR, "hifiasm_assembly.bp.p_ctg.gfa_mapped_longest_fraction")
    threads: 1
    log: os.path.join(OUTPUT_DIR, "logs/detect_clipping.log")
    shell:
        """
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        detect_clipping.py --bam {input.bam_long} --min_clip 10000 --add_coverage --out {params.prefix} > {log} 2>&1
        """


rule get_zero_coverage:
    input:
        bam_long=os.path.join(ONT_MAPPING_DIR, "hifiasm_assembly.bp.p_ctg.gfa_mapped_longest_fraction.sorted.bam"),
        bam_mid=os.path.join(ONT_MAPPING_DIR, "hifiasm_assembly.bp.p_ctg.gfa_mapped_midsized_fraction.sorted.bam")

    output:
        zero_coverage_bed_long=os.path.join(CLIPPING_DIR, "hifiasm_assembly.bp.p_ctg.gfa_mapped_longest_fraction_zero_coverage.bed"),
        zero_coverage_bed_mid=os.path.join(CLIPPING_DIR, "hifiasm_assembly.bp.p_ctg.gfa_mapped_midsized_fraction_zero_coverage.bed")
    conda:
        "envs/pysam.yaml"
    threads: 2
    log: "logs/get_zero_coverage.log"
    shell:
        """
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        zero_coverage.py --bam {input.bam_long} --out {output.zero_coverage_bed_long} &
        zero_coverage.py --bam {input.bam_mid} --out {output.zero_coverage_bed_mid} &
        wait
        """


rule extract_contig_ends:
    input:
        assembly=os.path.join(OUTPUT_DIR,"hifiasm_assembly.bp.p_ctg.gfa.fasta"),
        bam=os.path.join(ONT_MAPPING_DIR, "hifiasm_assembly.bp.p_ctg.gfa_mapped_longest_fraction.sorted.bam"),
        reads=os.path.join(DATA_DIR,"reads_longest_fraction.fasta")
    output:
        contig_pairs=os.path.join(CONTIG_PAIRS_DIR, "contig_pairs.csv")
    conda:
        "envs/pysam.yaml"
    params:
        output_dir = CONTIG_PAIRS_DIR
    threads: workflow.cores
    log: "logs/extract_contig_ends.log"

    shell:
        """
        scripts_dir=$(realpath scripts)
        export PATH=$scripts_dir:$PATH
        # get full path to extract_contig_ends2.py
        extract_script=$(which extract_contig_ends2.py)
        echo $extract_script

        get_contig_pairs_from_bam.py --bam {input.bam} --out {output.contig_pairs}
        run_extract_contig_ends_batch.py --pairs_file {output.contig_pairs} --bam {input.bam} --reads_fasta {input.reads} \
        --out_dir {params.output_dir} --genome_fasta {input.assembly} --min_aln_length 10000 --extract_script $extract_script \
        --cpu {threads}
        
        """



rule run_tidecluster:
    input:
        assembly=os.path.join(OUTPUT_DIR,"hifiasm_assembly.bp.p_ctg.gfa.fasta"),
        library=config['tandem_repeats']
    output:
        tidecluster=os.path.join(TIDECLUSTER_DIR, "tc_done.txt")
    params:
        prefix = os.path.join(TIDECLUSTER_DIR, "tc")
    conda:
        "envs/tidecluster.yaml"
    threads:
        workflow.cores
    shell:
        """
        TideCluster.py run_all -f {input.assembly} -pr {params.prefix} -l {input.library} -c {threads}
        touch {output.tidecluster}  
        """
