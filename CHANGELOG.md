# Changelog

All notable changes to this project will be documented in this file.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) — versioning follows [Semantic Versioning](https://semver.org/).

---

## [Unreleased]

---

## [1.0.0] — 2026-04-03

### Added
- Initial release of the modular ChIP-seq Snakemake pipeline.
- End-to-end workflow: QC → alignment → filtering → peak calling → downstream analysis → MultiQC report.
- `rules/scripts/validate_config.py` — config and sample sheet validation at startup.
- Singularity containers for all rules (Galaxy Project biocontainers depot).
- Conda environment YAMLs as fallback for all rules.
- `benchmarks/` and `logs/` directories automatically populated per rule per sample.
- Resource declarations (`mem_mb`, `time`, `threads`) on all rules for HPC compatibility.

### Rules included (v1.0.0)
`fastp_01` · `fastqc_02` · `bowtie2_03` · `samtools_sort_04` ·
`calculate_mito_reads_05` · `remove_mito_reads_06` ·
`samtools_fixmate_pre_08` · `samtools_markdup_08` ·
`samtools_index_after_markdup_09` · `samtools_view_10` ·
`samtools_index_post_filter_11` · `samtools_stats_12` ·
`fragment_size_analysis_13` · `picard_CollectAlignmentSummaryMetrics_14` ·
`picard_CollectInsertSizeMetrics_15` · `bedtools_genomecov_17` ·
`sorted_bedgraph_18` · `bigwig_conversion_19` · `correlation_analysis_20` ·
`normalize_coverage_21` · `macs2_peak_calling_22` ·
`blacklist_region_filter_23` · `heatmap_24` · `frip_calculation_25` ·
`peak_annotation_26` · `motif_analysis_27` · `phantompeakqualtool_28` ·
`preseq_29` · `qualimap_bamqc_30` · `multiqc_31`

### Fixed
- Resolved `KeyError` in `macs2_peak_calling_22.smk`: `config['macs_peakcall']['input']['shifted_bam']` → `filtered_bam`.
- Removed undefined `params.motif_db` directive from `motif_analysis_27.smk`.
- Uncommented `phantompeakqualtools` outputs in `rule all`.

### Removed
- `samtools_index_07.smk` — merged into downstream rules.
- `tss_enrichment_16.smk` and associated script — removed; not part of ChIP-seq standard workflow.
- Orphan files: `preseq_30.smk`, `qualimap_bamqc_31.smk`.

### Changed
- Renamed working directory from `atacseq/` to `chipseq_pipeline/` to correctly reflect pipeline type.
- Cleaned `envs/` folder; removed stale environment files.
- Fixed numbering mismatches across rule filenames in `rules/`.

---

[Unreleased]: https://github.com/<your-org>/chipseq_pipeline/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/<your-org>/chipseq_pipeline/releases/tag/v1.0.0
