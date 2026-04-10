# Modular ChIP-seq workflow
# Author: Himanshu Bhandary
# Contact: 2032ushimanshu@gmail.com

import csv
import subprocess
from pathlib import Path

configfile: "config.yaml"

subprocess.run(
    ["python3", "rules/scripts/validate_config.py", "config.yaml"],
    check=True,
)

SAMPLES_TSV = Path(config["global"]["samples"])
with SAMPLES_TSV.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

SAMPLES = [row["sample"] for row in rows]
FASTQ_R1 = {row["sample"]: row["fastq_r1"] for row in rows}
FASTQ_R2 = {row["sample"]: row["fastq_r2"] for row in rows}
CONTROLS = {row["sample"]: row["control"] for row in rows if row["control"] != "NONE"}
TREATMENTS = list(CONTROLS.keys())

if not SAMPLES:
    raise ValueError(f"No samples found in sample sheet: {SAMPLES_TSV}")


# Workflow modules
include: "rules/fastp.smk"
include: "rules/fastqc.smk"
include: "rules/bowtie2.smk"
include: "rules/samtools_sort.smk"
include: "rules/calculate_mito_reads.smk"
include: "rules/remove_mito_reads.smk"
include: "rules/samtools_fixmate.smk"
include: "rules/samtools_markdup.smk"
include: "rules/samtools_index_after_markdup.smk"
include: "rules/samtools_view.smk"
include: "rules/samtools_index_post_filter.smk"
include: "rules/samtools_stats.smk"
include: "rules/fragment_size_analysis.smk"
include: "rules/picard_alignment_metrics.smk"
include: "rules/picard_insert_size_metrics.smk"
include: "rules/bedtools_genomecov.smk"
include: "rules/sorted_bedgraph.smk"
include: "rules/bigwig_conversion.smk"
include: "rules/correlation_analysis.smk"
include: "rules/normalize_coverage.smk"
include: "rules/macs2_peak_calling.smk"
include: "rules/blacklist_region_filter.smk"
include: "rules/heatmap.smk"
include: "rules/frip_calculation.smk"
include: "rules/peak_annotation.smk"
include: "rules/motif_analysis.smk"
include: "rules/phantompeakqualtool.smk"
include: "rules/preseq.smk"
include: "rules/qualimap_bamqc.smk"
include: "rules/multiqc.smk"
include: "rules/qc_gate.smk"
include: "rules/bamCompare.smk"


# Pipeline targets grouped by stage
PREPROCESSING_TARGETS = (
    expand("results/fastp/{sample}_R1_trimmed.fastq.gz", sample=SAMPLES)
    + expand("results/fastp/{sample}_R2_trimmed.fastq.gz", sample=SAMPLES)
    + expand("results/fastp/{sample}.html", sample=SAMPLES)
    + expand("results/fastp/{sample}.json", sample=SAMPLES)
    + expand("results/fastqc/{sample}_R1_trimmed_fastqc.html", sample=SAMPLES)
    + expand("results/fastqc/{sample}_R1_trimmed_fastqc.zip", sample=SAMPLES)
    + expand("results/fastqc/{sample}_R2_trimmed_fastqc.html", sample=SAMPLES)
    + expand("results/fastqc/{sample}_R2_trimmed_fastqc.zip", sample=SAMPLES)
)

ALIGNMENT_TARGETS = expand("results/bowtie2/{sample}.bam", sample=SAMPLES)

POST_ALIGNMENT_TARGETS = (
    expand("results/samtools_sort/{sample}.sorted.bam", sample=SAMPLES)
    + expand("results/mito_ChIP/{sample}_mito_stats.txt", sample=SAMPLES)
    + expand("results/remove_mito_reads/{sample}_noMT.sorted.bam", sample=SAMPLES)
    + expand("results/samtools_fixmate/{sample}_noMT.sorted.fixmate.bam", sample=SAMPLES)
    + expand("results/samtools_markdup/{sample}_noMT.sorted.dedup.bam", sample=SAMPLES)
    + expand("results/samtools_index/post_markdup/{sample}_noMT.sorted.dedup.bam.bai", sample=SAMPLES)
    + expand("results/samtools_view/{sample}.filtered.bam", sample=SAMPLES)
    + expand("results/samtools_view/{sample}.filtered.bam.bai", sample=SAMPLES)
    + expand("results/samtools_stats/{sample}_postFiltering.stats.txt", sample=SAMPLES)
    + expand("results/fragment_size_analysis/{sample}_fragment_sizes.txt", sample=SAMPLES)
    + expand("results/fragment_size_analysis/{sample}_fragment.png", sample=SAMPLES)
    + expand("results/fragment_size_analysis/{sample}_fragment_stats.txt", sample=SAMPLES)
)

PICARD_METRICS_TARGETS = (
    expand("results/picard/CollectAlignmentSummaryMetrics/{sample}.alignment_metrics.txt", sample=SAMPLES)
    + expand("results/picard/CollectInsertSizeMetrics/{sample}.insert_metrics.txt", sample=SAMPLES)
    + expand("results/picard/CollectInsertSizeMetrics/{sample}.insert_histogram.pdf", sample=SAMPLES)
)

COVERAGE_TARGETS = (
    expand("results/bedtools_genomecov/{sample}.bedGraph", sample=SAMPLES)
    + expand("results/sorted_bedgraph_file/{sample}.sorted.bedGraph", sample=SAMPLES)
    + expand("results/bigwig/{sample}.bw", sample=SAMPLES)
    + expand("results/normalized_coverage/{sample}_CPM.bw", sample=SAMPLES)
)

CORRELATION_TARGETS = [
    "results/correlation_analysis/matrix.npz",
    "results/correlation_analysis/matrix.tab",
    "results/correlation_analysis/correlation_heatmap.png",
    "results/correlation_analysis/correlation_values.tab",
]

PEAK_CALLING_TARGETS = (
    expand("results/macs2_peakcall/{sample}_peaks.narrowPeak", sample=SAMPLES)
    + expand("results/filtered_peaks/{sample}_filtered_peaks.bed", sample=SAMPLES)
)

HEATMAP_AND_FRIP_TARGETS = (
    expand("results/heatmap/matrix/{sample}_matrix.gz", sample=SAMPLES)
    + expand("results/heatmap/{sample}_regions.bed", sample=SAMPLES)
    + expand("results/heatmap/plot/{sample}_tss_heatmap.pdf", sample=SAMPLES)
    + expand("results/frip_calculation/{sample}_frip.txt", sample=SAMPLES)
)

ANNOTATION_AND_MOTIF_TARGETS = (
    expand("results/peak_annotation/{sample}_peak_annotation.txt", sample=SAMPLES)
    + ["results/motif_analysis"]
)

QC_TARGETS = (
    expand("results/phantompeakqualtools/{sample}_qc.txt", sample=SAMPLES)
    + expand("results/phantompeakqualtools/{sample}_qc.pdf", sample=SAMPLES)
    + expand("results/preseq/{sample}.ccurve.txt", sample=SAMPLES)
    + expand("results/qualimap/{sample}_qualimap_report", sample=SAMPLES)
    + expand("results/qc_gate/{sample}_qc_pass.txt", sample=SAMPLES)
)

BAMCOMPARE_TARGETS = expand("results/bamCompare/{sample}_vs_control.bw", sample=TREATMENTS)

MULTIQC_TARGETS = ["results/multiqc"]


rule all:
    input:
        (
            PREPROCESSING_TARGETS
            + ALIGNMENT_TARGETS
            + POST_ALIGNMENT_TARGETS
            + PICARD_METRICS_TARGETS
            + COVERAGE_TARGETS
            + CORRELATION_TARGETS
            + PEAK_CALLING_TARGETS
            + HEATMAP_AND_FRIP_TARGETS
            + ANNOTATION_AND_MOTIF_TARGETS
            + QC_TARGETS
            + BAMCOMPARE_TARGETS
            + MULTIQC_TARGETS
        )
