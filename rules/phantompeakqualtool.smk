rule phantompeakqualtools:
    input:
        indexed_dedup=f"{config['samtools_index_post_markdup']['output']['index']}/{{sample}}_noMT.sorted.dedup.bam.bai",
        dedup_bam=f"{config['samtools_markdup']['output']['markdup_bam']}/{{sample}}_noMT.sorted.dedup.bam"
    output:
        qc_metrics=f"results/phantompeakqualtools/{{sample}}_qc.txt",
        qc_pdf=f"results/phantompeakqualtools/{{sample}}_qc.pdf"
    resources:
        mem_mb=16000,
        time="04:00:00"
    benchmark: "benchmarks/phantompeakqualtools/{sample}.txt"
    log: "logs/phantompeakqualtools/{sample}.err"
    conda: "envs/phantompeakqualtools.yaml"
    container: "https://depot.galaxyproject.org/singularity/phantompeakqualtools:1.2.2--0"
    threads: 8
    message:
        "[PhantomPeakQualTools] Sample: {wildcards.sample} | Quality control."
    shell:
        """
        Rscript run_spp.R -c={input.dedup_bam} -savp={output.qc_pdf} -out={output.qc_metrics} 2> {log}
        """
