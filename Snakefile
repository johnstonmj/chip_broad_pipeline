# # Run
# snakemake -nq --jobs 200 --use-conda --cluster-config cluster.yaml --cluster "bsub -J {cluster.jobname} -n {cluster.numcpu} -R {cluster.span} -R {cluster.memory} -M {cluster.maxmem} -We {cluster.wall_est} -W {cluster.wall_max} -o {cluster.output} -e {cluster.error} < " all

# snakemake --dag all | dot -Tsvg > chip_broad_dag.svg

# Path to config file
configfile: "config.yaml"

# Uses 'chain' from itertools extract all elements from 'list of lists'
# 'set' removes duplicates before converting back to 'list'
all_samples = list(set(chain.from_iterable(config["ip_input_pairs"])))
ip_samples = [ x[0] for x in config["ip_input_pairs"] ]
input_samples = [ x[1] for x in config["ip_input_pairs"] ]

# Basename for the last part of the filename
# Extract genome_name from zeroth part of split filename
genome_name = os.path.splitext(os.path.basename(config["genome_path"]))[0]

# Don't submit full server jobs for these rules
localrules: all

# Collector rule
rule all:
    input: 
        expand("fastqc_out/{sample}.{read}/{sample}.{read}_fastqc.html", sample = all_samples, read = ['R1','R2']),
        # expand("bwa_alignment/{sample}." + genome_name + ".bam.bai", sample = all_samples),
        # Note the 'zip'ping of these two lists rather than the default 'product'
        expand("macs_peaks/{treat_samp}_WRT_{ctrl_samp}_FE.tdf", zip, treat_samp = ip_samples, ctrl_samp = input_samples),
        # expand("pearson/{sample}." + genome_name + ".norm.bw", sample = all_samples),
        # "pearson/multiBigwig_norm_counts.tsv",
        "pearson/pearsonCorrelation_ChIP.tsv",


# FastQC Quality check
#########
rule fastqc:
    message: "Running FastQC sequencing quality stats for: {wildcards.sample}"
    input:
        "input_fastqs/{sample}.{read}.fastq"
    output:
        "fastqc_out/{sample}.{read}/{sample}.{read}_fastqc.html"
    conda:
        "envs/chip.yaml"
    params:
        tmp_dir = config["tmp_dir"]
    threads:
        4
    shell:
        """
        # Make ouput directory because fastqc won't on its own
        # -p prevents errors if already exists
        mkdir -p fastqc_out/{wildcards.sample}.{wildcards.read}

        # Run fastqc
        fastqc --outdir fastqc_out/{wildcards.sample}.{wildcards.read} --dir {params.tmp_dir} --threads {threads} {input}
        """

# Align with bwa mem
# Filter and sort bam file
#########
rule bwa_align:
    message: "Aligning reads, filtering reads, sorting, and indexing for: {wildcards.sample}"
    input:
        fq1 = "input_fastqs/{sample}.R1.fastq",
        fq2 = "input_fastqs/{sample}.R2.fastq",
        ref = config["genome_path"],
    output:
        sam = "bwa_alignment/{sample}." + genome_name + ".sam",
        bam = "bwa_alignment/{sample}." + genome_name + ".bam",
        bai = "bwa_alignment/{sample}." + genome_name + ".bam.bai",
    conda:
        "envs/chip.yaml"
    params:
        blacklist = config["blacklist"],
        min_aln_q = config["min_aln_q"]
    threads:
        16
    shell:
        """
        # bwa mem align
        bwa mem -t {threads} {input.ref} {input.fq1} {input.fq2} > {output.sam}

        # Read sam file
        samtools view --threads {threads} -hS {output.sam} |\
        # Filter out unwanted chromosome names
        grep -E -v "chrUn|random|_alt" |\
        # Filter out low quality
        samtools view --threads {threads} -bS -q {params.min_aln_q} - |\
        # Filter out blacklisted regions
        bedtools intersect -v -abam stdin -b {params.blacklist} |\
        # Sort by position prior to indexing
        samtools sort --threads {threads} - > {output.bam}

        # Index
        samtools index {output.bam}
        """

# Call broad peaks with MACS2
#########
rule macs2_call_broadpeaks:
    message: "Calling broad peaks using MACS2: {wildcards.treat_samp}_WRT_{wildcards.ctrl_samp}"
    input:
        ctrl_bam = "bwa_alignment/{ctrl_samp}." + genome_name + ".bam",
        treat_bam = "bwa_alignment/{treat_samp}." + genome_name + ".bam",
    output:
        bed = "macs_peaks/{treat_samp}_WRT_{ctrl_samp}_peaks.broadPeak.bed",
        broadpeaks = "macs_peaks/{treat_samp}_WRT_{ctrl_samp}_peaks.broadPeak",
        control = "macs_peaks/{treat_samp}_WRT_{ctrl_samp}_control_lambda.bdg",
        treat = "macs_peaks/{treat_samp}_WRT_{ctrl_samp}_treat_pileup.bdg",
    conda:
        "envs/chip.yaml"
    params:
        tmp_dir = config["tmp_dir"],
        genome_size = config["genome_size"],
        q_val = config["q_val"]
    shell:
        """
        # MACS2 arguments
        # -t and -c: "treat" and "control" bams
        # --broad: Broadpeaks instead of narrow.
        # -f: Input file format
        # -n: Output prefix
        # -g: Genome size = hs for homo sapiens
        # -B / --bdg: Save bedgraph of pileups
        # --SPMR: signal per million reads. 
        # --SPMR Recommended for displaying normalized pileup tracks across many datasets. Requires -B to be set.
        # --SPMR Not recommended for bdgcmp and bdgpeakcall as results will differ 
        macs2 callpeak -t "{input.treat_bam}" -c "{input.ctrl_bam}" --broad -f "BAMPE" \
        -n "macs_peaks/{wildcards.treat_samp}_WRT_{wildcards.ctrl_samp}" \
        --tempdir "{params.tmp_dir}" -g "{params.genome_size}" -q "{params.q_val}" -B --SPMR

        # Make additional simple bed file from broadpeak results
        cut -f 1-3 "{output.broadpeaks}" > "{output.bed}"
        """


# Generate fold enrichment traces with MACS2 and igvtools
#########
rule foldEnrich:
    message: "Generating fold enrichment traces for: {wildcards.treat_samp}_WRT_{wildcards.ctrl_samp}"
    input:
        control = "macs_peaks/{treat_samp}_WRT_{ctrl_samp}_control_lambda.bdg",
        treat = "macs_peaks/{treat_samp}_WRT_{ctrl_samp}_treat_pileup.bdg",
    output:
        FEbdg = "macs_peaks/{treat_samp}_WRT_{ctrl_samp}_FE.bdg",
        FEbdg_sorted = temp("macs_peaks/{treat_samp}_WRT_{ctrl_samp}_FE.sorted.bdg"),
        FEtdf = "macs_peaks/{treat_samp}_WRT_{ctrl_samp}_FE.tdf",
    conda:
        "envs/chip.yaml"
    params:
        tmp_dir = config["tmp_dir"],
        chrom_sizes = config["chrom_sizes"],
    shell:
        """
        # MACS2 arguments
        # -t and -c: "treat" and "control" bams
        # -o: Output filename
        # -m: Method of comparison. FE = fold enrichment
        macs2 bdgcmp -t "{input.treat}" -c "{input.control}" -o "{output.FEbdg}" -m FE

        # Convert to binary format for efficient IGV loading
        # -z: Sets maximum zoom level computed.
        igvtools sort -t "{params.tmp_dir}" "{output.FEbdg}" "{output.FEbdg_sorted}"
        igvtools toTDF -z 10 "{output.FEbdg_sorted}" "{output.FEtdf}" "{params.chrom_sizes}"
        """

# Normalize coverage of each sample with DeepTools
###########
rule coverage_normalize:
    message: "Using DeepTools to normalize ChIP signal for comparison between samples."
    input:
        bam = "bwa_alignment/{sample}." + genome_name + ".bam",
    output:
        norm_bw = "pearson/{sample}." + genome_name + ".norm.bw",
    conda:
        "envs/chip.yaml"
    threads:
        4
    params:
        eff_genome = config["eff_genome"],
    shell:
        """
        # bamCoverage to normalize
        # Effective genome size taken from https://www.nature.com/articles/nbt.1518/tables/1  @70bp seq length
        # RPGC = reads per genomic content (1x normalization)
        bamCoverage --bam {input.bam} \
            --outFileName {output.norm_bw} \
            --effectiveGenomeSize {params.eff_genome} \
            --normalizeUsing RPGC \
            --extendReads \
            --numberOfProcessors {threads}
        """

# Assemble normalized count table using DeepTools
##########
rule multiBigwig_norm_counts:
    message: "Using DeepTools to create a summary table containing all samples."
    input:
        all_broadpeaks = expand("macs_peaks/{treat_samp}_WRT_{ctrl_samp}_peaks.broadPeak", zip, treat_samp = ip_samples, ctrl_samp = input_samples),
        all_norm = expand("pearson/{sample}." + genome_name + ".norm.bw", sample = all_samples),
    output:
        union_bed = "pearson/union_broadpeaks.bed",
        npz = "pearson/multiBigwigSummary.npz",
        counts = "pearson/multiBigwig_norm_counts.tsv",
    conda:
        "envs/chip.yaml"
    params:
        labels = expand("{sample}", sample = all_samples),
    threads:
        4
    shell:
        """
        # Generate consensus file of all regions identified in any sample
        cat {input.all_broadpeaks} | bedtools sort -i stdin | bedtools merge -i stdin > {output.union_bed}

        # Summarize multipe bigwig files together in one table
        multiBigwigSummary BED-file --bwfiles {input.all_norm} \
            --outFileName {output.npz} \
            --BED {output.union_bed} \
            --labels {params.labels} \
            --numberOfProcessors {threads} \
            --outRawCounts {output.counts}
        """

# Generate Pearson Plot
# Uses R script accessing snakemake S4 object within R
##################
rule pearson:
    message: "Using R to calculate and plot Pearson correlation between samples."
    input:
        counts = "pearson/multiBigwig_norm_counts.tsv",
    output:
        pdf = "pearson/pearsonCorrelation_ChIP.pdf",
        tsv = "pearson/pearsonCorrelation_ChIP.tsv",
    conda:
        "envs/r-pearson.yaml"
    threads:
        4
    script:
        "scripts/pearson.R"