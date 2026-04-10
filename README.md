# Modular ChIP-seq Pipeline

A scalable, reproducible, and modular [Snakemake](https://snakemake.readthedocs.io/) pipeline for end-to-end processing of paired-end ChIP-seq data. Starting from raw FASTQ files, the pipeline performs quality control, alignment, duplicate removal, peak calling, coverage generation, control-aware signal normalization, and downstream analysis.

---

## Features

- **Reproducible** — All tools run inside strictly versioned [Galaxy Project Singularity containers](https://depot.galaxyproject.org/singularity/) or Conda environments. Results are byte-for-byte reproducible.
- **Modular** — Each analysis step lives in its own `rules/*.smk` file. Add, remove, or swap steps without touching the master `Snakefile`.
- **Config-driven** — A single `config.yaml` controls all parameters. No hardcoded values inside rules.
- **Validated** — `rules/scripts/validate_config.py` checks the config and sample sheet at startup, catching errors before any job is submitted.
- **Control-aware** — Matched control samples can be used for MACS2 peak calling and `bamCompare` signal tracks.
- **QC-gated** — Key QC metrics can fail the workflow if thresholds are not met.
- **Cluster-ready** — All rules declare `threads`, `mem_mb`, and `time` resources, making it compatible with SLURM, SGE, and PBS via Snakemake profiles.

---

## Pipeline Overview

```
Raw FASTQ
   │
   ├─ 01 fastp          — Adapter trimming & quality filtering
   ├─ 02 FastQC         — Per-sample QC report (post-trim)
   │
   ├─ 03 Bowtie2        — Paired-end alignment (--very-sensitive)
   ├─ 04 samtools sort  — Coordinate-sorted BAM
   │
   ├─ 05 Mito filter    — Count & remove mitochondrial reads
   ├─ 06 samtools fixmate
   ├─ 07 samtools markdup — Duplicate removal
   ├─ 08 samtools index
   ├─ 09 samtools view  — MAPQ ≥ 30, flag filter (3844)
   ├─ 10 samtools index
   ├─ 11 samtools stats — Post-filtering alignment statistics
   ├─ 12 Fragment size analysis
   │
   ├─ 13 Picard AlignmentSummaryMetrics
   ├─ 14 Picard InsertSizeMetrics
   │
   ├─ 15 bedtools genomecov — Genome coverage (BedGraph)
   ├─ 16 sort BedGraph
   ├─ 17 bedGraphToBigWig — BigWig for genome browsers
   ├─ 18 deepTools bamCompare — log2(ChIP / control) signal track
   │
   ├─ 19 deepTools multiBigwigSummary — Correlation matrix
   │
   ├─ 20 MACS2          — Narrow peak calling
   ├─ 21 Blacklist filter — Remove ENCODE blacklist regions
   │
   ├─ 22 deepTools computeMatrix / plotHeatmap
   ├─ 23 FRiP calculation
   ├─ 24 QC gate         — Threshold-based pass/fail checks
   │
   ├─ 25 ChIPseeker      — Peak annotation
   ├─ 26 HOMER           — De novo motif analysis
   │
   ├─ 27 PhantomPeakQualTools — NSC / RSC strand cross-correlation
   ├─ 28 Preseq          — Library complexity estimation
   ├─ 29 Qualimap bamqc  — Comprehensive BAM QC
   │
   └─ 30 MultiQC         — Aggregated QC report
```

---

## Requirements

| Software | Version | Purpose |
|----------|---------|---------|
| Snakemake | ≥ 8.0 | Workflow manager |
| Singularity / Apptainer | ≥ 3.8 | Container runtime |
| Conda / Mamba | any | Environment manager |
| Python | ≥ 3.10 | Config validation script |

---

## Installation

```bash
# Create and activate a Snakemake environment
conda create -n snakemake -c bioconda -c conda-forge snakemake singularity graphviz -y
conda activate snakemake

# Clone the pipeline
git clone https://github.com/<your-org>/chipseq_pipeline.git
cd chipseq_pipeline
```

---

## Quick Start

### 1. Prepare reference files

Place the following in `data/reference/`:

| File | Description |
|------|-------------|
| `genome.fa` | Reference genome FASTA |
| `genome.chrom.sizes` | Chromosome sizes (from `samtools faidx` + `cut -f1,2`) |
| `ENCODE_blacklist.bed` | ENCODE blacklist (e.g., [hg38](https://github.com/Boyle-Lab/Blacklist)) |
| `annotation.gtf` | Gene annotation GTF (Ensembl or GENCODE) |

Build the Bowtie2 index:

```bash
bowtie2-build data/reference/genome.fa data/02_alignment/bowtie2/index/genome
```

### 2. Define your samples

Create `data/fastp/samples.tsv` (tab-separated):

```tsv
sample	condition	replicate	fastq_r1	fastq_r2	control
Sample1	Control	1	/path/to/Sample1_R1.fastq.gz	/path/to/Sample1_R2.fastq.gz	NONE
Sample2	Control	2	/path/to/Sample2_R1.fastq.gz	/path/to/Sample2_R2.fastq.gz	NONE
Sample3	Treatment	1	/path/to/Sample3_R1.fastq.gz	/path/to/Sample3_R2.fastq.gz	Sample1
Sample4	Treatment	2	/path/to/Sample4_R1.fastq.gz	/path/to/Sample4_R2.fastq.gz	Sample2
```

> **Note**: Use absolute paths for FASTQ files.
>
> Use `NONE` in the `control` column for control-only rows.

### 3. Configure parameters

Edit `config.yaml` — at minimum, verify:
- `global.bowtie_index` — path to your Bowtie2 index prefix
- `global.genome_fa` — path to reference FASTA
- `global.genome_chrom_sizes` — chrom sizes file
- `mito_ChIP_calculate.params.mito_chr` — `"MT"` (Ensembl) or `"chrM"` (UCSC)
- `macs_peakcall.params.genome_size` — `"hs"` (human), `"mm"` (mouse), etc.
- `qc_gate.params.*` — FRiP, NSC, RSC, mapping-rate, and duplicate-rate thresholds

### 4. Run the pipeline

```bash
# Dry run (check DAG without executing)
snakemake --dry-run --cores 8

# Run locally with Singularity
snakemake --use-singularity --cores 8

# Run locally with Conda
snakemake --use-conda --cores 8

# Run on SLURM cluster
snakemake --use-singularity --profile slurm --jobs 50
```

---

## Output Structure

```
results/
├── fastp/                      # Trimmed FASTQs + fastp reports
├── fastqc/                     # FastQC HTML + ZIP reports
├── bowtie2/                    # Raw aligned BAMs
├── samtools_sort/              # Coordinate-sorted BAMs
├── mito_ChIP/                  # Mitochondrial read statistics
├── remove_mito_reads/          # MT-filtered BAMs
├── samtools_fixmate/           # Fixmate BAMs
├── samtools_markdup/           # Deduplicated BAMs
├── samtools_index/             # BAM indices
├── samtools_view/              # MAPQ/flag filtered BAMs
├── samtools_stats/             # Alignment statistics
├── fragment_size_analysis/     # Fragment size distributions & plots
├── picard/
│   ├── CollectAlignmentSummaryMetrics/
│   └── CollectInsertSizeMetrics/
├── bedtools_genomecov/         # Raw BedGraphs
├── sorted_bedgraph_file/       # Sorted BedGraphs
├── bigwig/                     # Raw BigWig tracks
├── normalized_coverage/        # CPM-normalized BigWigs
├── bamCompare/                 # ChIP vs control log2 signal tracks
├── correlation_analysis/       # Sample correlation matrix & heatmap
├── macs2_peakcall/             # narrowPeak files
├── filtered_peaks/             # Blacklist-filtered peaks
├── heatmap/                    # TSS heatmaps (matrix + PDF)
├── frip_calculation/           # FRiP scores per sample
├── qc_gate/                    # QC pass/fail markers
├── peak_annotation/            # ChIPseeker annotation tables
├── motif_analysis/             # HOMER motif discovery output
├── phantompeakqualtools/       # NSC/RSC metrics
├── preseq/                     # Library complexity curves
├── qualimap/                   # BAM QC reports
└── multiqc/                    # Aggregated MultiQC report ← start here
```

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fastp.params.trim_front1/2` | 5 | Bases trimmed from 5' end (R1/R2) |
| `fastp.params.length_required` | 30 | Minimum read length after trimming |
| `bowtie2.sensitive` | `--very-sensitive` | Alignment sensitivity preset |
| `mito_ChIP_calculate.params.mito_chr` | `MT` | Mitochondrial chromosome name |
| `samtools_view.params.MAPQ` | 30 | Minimum mapping quality |
| `samtools_view.params.flags` | 3844 | SAM flags to exclude |
| `macs_peakcall.params.genome_size` | `hs` | Effective genome size |
| `macs_peakcall.params.qvalue` | 0.01 | MACS2 peak calling q-value threshold |
| `macs_peakcall.params.format` | `BAMPE` | Paired-end BAM input format |
| `qc_gate.params.min_frip` | 0.05 | Minimum FRiP threshold |
| `qc_gate.params.min_nsc` | 1.05 | Minimum NSC threshold |
| `qc_gate.params.min_rsc` | 0.8 | Minimum RSC threshold |
| `qc_gate.params.min_mapping_rate` | 90.0 | Minimum mapping-rate threshold |
| `qc_gate.params.max_duplicate_rate` | 20.0 | Maximum duplicate-rate threshold |

---

## Software Versions (Singularity Containers)

| Tool | Version | Container |
|------|---------|-----------|
| fastp | 1.1.0 | `fastp:1.1.0--heae3180_0` |
| FastQC | — | `fastqc` |
| Bowtie2 | 2.5.4 | `bowtie2:2.5.4--he96a11b_5` |
| samtools | — | `samtools` |
| Picard | — | `picard` |
| bedtools | 2.31.1 | `bedtools:2.31.1--h13024bc_3` |
| MACS2 | 2.2.9.1 | `macs2:2.2.9.1--py39hbcbf7aa_2` |
| deepTools | — | `deeptools` |
| HOMER | 4.11 | `homer:4.11--pl5262h4ac6f70_9` |
| ChIPseeker | 1.46.1 | `bioconductor-chipseeker:1.46.1--r45hdfd78af_0` |
| PhantomPeakQualTools | 1.2.2 | `phantompeakqualtools:1.2.2--0` |
| MultiQC | — | `multiqc` |

All containers are pulled from the [Galaxy Project Singularity depot](https://depot.galaxyproject.org/singularity/).

---

## Authors

**Himanshu Bhandary**
Email: [2032ushimanshu@gmail.com](mailto:2032ushimanshu@gmail.com)

---

## License

MIT License. See [`LICENSE`](LICENSE) for full terms.

---

## Citation

> Bhandary H. *et al.* (2026). Modular ChIP-seq Pipeline [Software]. GitHub. https://github.com/<your-org>/chipseq_pipeline

*(Update with journal citation once published.)*
