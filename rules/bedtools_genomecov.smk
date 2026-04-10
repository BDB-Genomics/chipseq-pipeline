rule bedtools_genomecov:
    input:
        shifted_bam=lambda wildcards: f"{config['bedtools_genomecov']['input']['filtered_bam']}/{wildcards.sample}.filtered.bam"
 
    output:
        bedgraph=f"{config['bedtools_genomecov']['output']['bedgraph']}/{{sample}}.bedGraph"

    params:
        extra=config['bedtools_genomecov']['params']['extra']

    resources:
        mem_mb=config['bedtools_genomecov']['resources']['mem_mb'], 
        time=config['bedtools_genomecov']['resources']['time']

    benchmark: "benchmarks/bedtools_genomecov/{sample}.txt"
    log: "logs/bedtools_genomecov/{sample}.err"
    conda: "envs/bedtools.yaml"
    container: "https://depot.galaxyproject.org/singularity/bedtools:2.31.1--h13024bc_3"
    threads: config['bedtools_genomecov']['threads']

    message:
        "[bedtools genomecov] sample: {wildcards.sample} | BAM : {input.shifted_bam}| Output: {output.bedgraph}..."

    shell:
        """
        bedtools genomecov \
          -ibam {input.shifted_bam} \
          {params.extra} \
          > {output.bedgraph} \
          2> {log}
        """
