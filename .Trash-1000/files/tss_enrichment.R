suppressPackageStartupMessages({
    library(ATACseqQC)
    library(GenomicFeatures)
    library(GenomicAlignments)
    library(Rsamtools)
    library(ChIPpeakAnno)
})

# Get command line arguments from snakemake
bamfile <- snakemake@input[["shifted_bam"]]
out_text <- snakemake@output[["text"]]
out_pdf <- snakemake@output[["pdf"]]
sample_name <- snakemake@wildcards[["sample"]]
txdb_pkg <- snakemake@params[["txdb"]]
upstream <- as.numeric(snakemake@params[["upstream"]])
downstream <- as.numeric(snakemake@params[["downstream"]])

cat("===========================================\n")
cat("TSS Enrichment Analysis\n")
cat("===========================================\n")
cat("Sample:", sample_name, "\n")
cat("BAM file:", bamfile, "\n")
cat("TxDb:", txdb_pkg, "\n")
cat("===========================================\n\n")

# Load the specified TxDb package
library(txdb_pkg, character.only = TRUE)
txdb <- get(txdb_pkg)

# Check BAM chromosome names
cat("Checking BAM chromosome naming style...\n")
bam_chroms <- scanBamHeader(bamfile)[[1]]$targets
bam_chr_names <- names(bam_chroms)
cat("BAM chromosomes (first 5):", paste(head(bam_chr_names, 5), collapse=", "), "\n")

# Determine BAM naming style
has_chr_prefix <- any(grepl("^chr", bam_chr_names))
cat("BAM uses", ifelse(has_chr_prefix, "UCSC", "NCBI/Ensembl"), "style naming\n\n")

# Load transcripts and create TSS regions
cat("Loading transcript database...\n")
txs <- transcripts(txdb)

# Create TSS regions BEFORE converting chromosome names
cat("Creating TSS regions...\n")
tss_regions <- promoters(txs, upstream=upstream, downstream=downstream)

# Convert TSS regions to match BAM chromosome naming
if(has_chr_prefix) {
    cat("Converting TSS regions to UCSC style (with 'chr' prefix)...\n")
    seqlevelsStyle(tss_regions) <- "UCSC"
    standard_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
} else {
    cat("Converting TSS regions to NCBI style (no 'chr' prefix)...\n")
    seqlevelsStyle(tss_regions) <- "NCBI"
    standard_chroms <- c(as.character(1:22), "X", "Y", "MT")
}

# Keep only standard chromosomes that exist in the annotation
available_chroms <- intersect(seqlevels(tss_regions), standard_chroms)
tss_regions <- keepSeqlevels(tss_regions, available_chroms, pruning.mode="coarse")

cat("TSS regions after filtering:", length(tss_regions), "regions\n")
cat("Using chromosomes:", paste(head(seqlevels(tss_regions), 10), collapse=", "), "\n\n")

# Read BAM file
cat("Reading BAM file...\n")
bamParam <- ScanBamParam(
    flag = scanBamFlag(isProperPair = TRUE, 
                      isUnmappedQuery = FALSE,
                      isSecondaryAlignment = FALSE,
                      isSupplementaryAlignment = FALSE),
    what = c("qname", "flag", "mapq")
)

gal <- readGAlignments(bamfile, use.names = TRUE, param = bamParam)
cat("Read", format(length(gal), big.mark=","), "alignments\n")

# Check overlap between BAM and TSS regions
cat("\nChecking chromosome overlap...\n")
gal_chroms <- unique(as.character(seqnames(gal)))
tss_chroms <- unique(as.character(seqnames(tss_regions)))
common_chroms <- intersect(gal_chroms, tss_chroms)

cat("BAM chromosomes:", length(gal_chroms), "\n")
cat("TSS chromosomes:", length(tss_chroms), "\n")
cat("Common chromosomes:", length(common_chroms), "\n")

if(length(common_chroms) == 0) {
    stop("ERROR: No common chromosomes between BAM and TSS regions!\n",
         "BAM has: ", paste(head(gal_chroms, 5), collapse=", "), "\n",
         "TSS has: ", paste(head(tss_chroms, 5), collapse=", "))
}

cat("Common chromosomes:", paste(head(common_chroms, 10), collapse=", "), "\n\n")

# Keep only common chromosomes in both objects
gal <- keepSeqlevels(gal, common_chroms, pruning.mode="coarse")
tss_regions <- keepSeqlevels(tss_regions, common_chroms, pruning.mode="coarse")

cat("After filtering to common chromosomes:\n")
cat("  Alignments:", format(length(gal), big.mark=","), "\n")
cat("  TSS regions:", format(length(tss_regions), big.mark=","), "\n\n")

# Calculate TSS enrichment score
cat("Calculating TSS enrichment score...\n")

tsse_score <- tryCatch({
    tsse_result <- TSSEscore(gal, txs = tss_regions)
    score <- if(is.list(tsse_result)) tsse_result$TSSEscore else tsse_result
    cat("TSS Enrichment Score:", round(score, 4), "\n\n")
    score
}, error = function(e) {
    cat("ERROR in TSSEscore:", conditionMessage(e), "\n")
    cat("Attempting alternative calculation...\n")
    
    # Alternative: calculate enrichment manually
    reads_gr <- as(gal, "GRanges")
    overlaps <- countOverlaps(tss_regions, reads_gr)
    
    tss_coverage <- sum(overlaps > 0) / length(tss_regions)
    tss_avg_depth <- mean(overlaps[overlaps > 0])
    
    genome_size <- sum(as.numeric(bam_chroms[common_chroms]))
    genome_avg <- length(reads_gr) / (genome_size / median(width(tss_regions)))
    
    score <- tss_avg_depth / max(genome_avg, 1)
    cat("Manual TSS Enrichment Score:", round(score, 4), "\n\n")
    score
})

# Determine quality
if(tsse_score > 7) {
    quality <- "Excellent"
    color <- "darkgreen"
} else if(tsse_score > 5) {
    quality <- "Good" 
    color <- "green"
} else if(tsse_score > 3) {
    quality <- "Acceptable"
    color <- "orange"
} else {
    quality <- "Poor"
    color <- "red"
}

cat("Quality Assessment:", quality, "\n")
cat("===========================================\n\n")

# Create publication-quality multi-panel plot
cat("Generating comprehensive TSS enrichment plots...\n")
pdf(out_pdf, width=14, height=10)

# Set up layout for multiple plots
layout(matrix(c(1,1,2,3), 2, 2, byrow=TRUE), heights=c(1.2, 1))
par(mar=c(5, 5, 4, 2))

# PANEL 1: TSS Enrichment Profile (Main Plot)
cat("Panel 1: Computing TSS enrichment profile...\n")
tryCatch({
    reads <- coverage(gal)
    
    # Calculate signal around TSS - sample if too many regions
    n_sample <- min(length(tss_regions), 10000)
    tss_sample <- tss_regions
    if(length(tss_regions) > n_sample) {
        set.seed(42)
        tss_sample <- tss_regions[sample(length(tss_regions), n_sample)]
        cat("  Sampling", n_sample, "TSS regions for visualization\n")
    }
    
    cat("DEBUG: Starting featureAlignedSignal with", length(tss_sample), "TSS regions and", length(reads), "reads\n")
    
    sigs <- featureAlignedSignal(reads, 
                                 feature.gr = tss_sample,
                                 upstream = upstream,
                                 downstream = downstream,
                                 n.tile = 100)
    
    # Average signal
    avg_signal <- colMeans(sigs, na.rm = TRUE)
    
    # Normalize to baseline
    baseline <- mean(c(avg_signal[1:10], avg_signal[91:100]), na.rm = TRUE)
    if(baseline > 0 && !is.na(baseline)) {
        avg_signal <- avg_signal / baseline
    }
    
    # Create main plot with enhanced styling
    plot(avg_signal, 
         type = "l", 
         lwd = 4, 
         col = "#2E86AB",
         xlab = "Distance from TSS (bp)",
         ylab = "Normalized Read Density",
         main = paste0(sample_name, " - TSS Enrichment Profile\nScore: ", 
                      round(tsse_score, 3), " (", quality, ")"),
         xaxt = "n",
         las = 1,
         ylim = c(0, max(avg_signal, na.rm=TRUE) * 1.15),
         cex.lab = 1.4,
         cex.axis = 1.2,
         cex.main = 1.5)
    
    # Add shaded confidence region (optional)
    se_signal <- apply(sigs, 2, sd, na.rm=TRUE) / sqrt(nrow(sigs))
    polygon(c(1:100, 100:1),
            c(avg_signal + se_signal, rev(avg_signal - se_signal)),
            col = rgb(0.18, 0.53, 0.67, 0.2), border = NA)
    
    # Redraw main line
    lines(avg_signal, lwd = 4, col = "#2E86AB")
    
    # Add custom x-axis
    axis(1, at = seq(1, 100, by = 25),
         labels = round(seq(-upstream, downstream, length.out = 5)),
         cex.axis = 1.2)
    
    # Add vertical line at TSS
    abline(v = 50, lty = 2, col = "red", lwd = 2)
    
    # Add horizontal baseline
    if(baseline > 0) {
        abline(h = 1, lty = 3, col = "gray40", lwd = 1.5)
    }
    
    # Add grid
    grid(col = "gray70", lty = "dotted")
    
    # PANEL 2: Heatmap of individual TSS regions
    cat("Panel 2: Creating TSS heatmap...\n")
    par(mar=c(5, 5, 4, 2))
    
    # Sample and sort TSS regions by signal strength
    n_heatmap <- min(500, nrow(sigs))
    signal_strength <- rowMeans(sigs[, 45:55], na.rm=TRUE)
    top_idx <- order(signal_strength, decreasing=TRUE)[1:n_heatmap]
    
    heatmap_data <- sigs[top_idx, ]
    
    # Create color palette
    colors <- colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100)
    
    # Plot heatmap
    image(t(heatmap_data[nrow(heatmap_data):1, ]),
          col = colors,
          xlab = "Distance from TSS (bp)",
          ylab = paste0("TSS Regions (top ", n_heatmap, " by signal)"),
          main = "TSS Signal Heatmap",
          axes = FALSE,
          cex.lab = 1.3,
          cex.main = 1.3)
    
    # Add axes
    axis(1, at = seq(0, 1, by = 0.25),
         labels = round(seq(-upstream, downstream, length.out = 5)),
         cex.axis = 1.1)
    axis(2, las = 1, cex.axis = 1.1)
    
    # Add TSS line
    abline(v = 0.5, col = "white", lwd = 2, lty = 2)
    
    # Add color scale legend
    legend_breaks <- seq(min(heatmap_data, na.rm=TRUE), 
                        max(heatmap_data, na.rm=TRUE), 
                        length.out = 5)
    legend("topright", 
           legend = round(legend_breaks, 2),
           fill = colorRampPalette(colors)(5),
           title = "Signal",
           cex = 0.9,
           bg = "white")
    
    # PANEL 3: Summary statistics
    cat("Panel 3: Creating summary plot...\n")
    par(mar=c(5, 5, 4, 2))
    
    # Calculate signal distribution
    peak_signal <- avg_signal[45:55]
    flanking_signal <- c(avg_signal[1:20], avg_signal[80:100])
    
    boxplot(list("TSS Peak\n(-250 to +250bp)" = peak_signal,
                 "Flanking Regions\n(±1500-2000bp)" = flanking_signal),
            col = c("#2E86AB", "#95D5D8"),
            main = "Signal Distribution",
            ylab = "Normalized Read Density",
            las = 1,
            cex.lab = 1.3,
            cex.axis = 1.2,
            cex.main = 1.3)
    
    # Add enrichment score
    text(1.5, max(c(peak_signal, flanking_signal)) * 0.9,
         paste0("Enrichment: ", round(tsse_score, 2), "x\n",
                "Quality: ", quality),
         cex = 1.4, col = color, font = 2)
    
    grid(col = "gray70", lty = "dotted")
    
}, error = function(e) {
    cat("Error in creating plots:", conditionMessage(e), "\n")
    
    # Fallback: simple plot
    plot(1, type = "n", 
         xlim = c(-upstream, downstream),
         ylim = c(0, max(2, tsse_score)),
         xlab = "Distance from TSS (bp)",
         ylab = "Enrichment Score",
         main = paste0(sample_name, "\nTSS Score: ", round(tsse_score, 3)),
         las = 1,
         cex.lab = 1.3,
         cex.main = 1.5)
    
    text(0, tsse_score/2, 
         paste0("TSS Enrichment Score:\n", round(tsse_score, 3)),
         cex = 2.5, col = color, font = 2)
    
    text(0, tsse_score * 0.75,
         paste0("Quality: ", quality),
         cex = 2, col = color)
    
    abline(v = 0, lty = 2, col = "red", lwd = 2)
    grid(col = "gray70")
})

dev.off()
cat("Plots saved to:", out_pdf, "\n\n")

# Save results
cat("Saving results...\n")
result_df <- data.frame(
    Sample = sample_name, 
    TSS_Enrichment = round(tsse_score, 4),
    Total_Alignments = length(gal),
    TSS_Regions = length(tss_regions),
    Common_Chromosomes = length(common_chroms),
    Quality = quality,
    stringsAsFactors = FALSE
)

write.table(result_df, file=out_text, sep="\t", 
            quote=FALSE, row.names=FALSE)

cat("Results saved to:", out_text, "\n")
cat("===========================================\n")
cat("TSS Enrichment Analysis Complete!\n")
cat("Final Score:", round(tsse_score, 4), "-", quality, "\n")
cat("===========================================\n")
