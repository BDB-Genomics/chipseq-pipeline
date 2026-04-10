rule calculate_mito_reads:
    input:
        sorted_bam=lambda wildcards: f"{config['mito_ChIP_calculate']['input']['sorted_bam']}/{wildcards.sample}.sorted.bam"
        
    output:
        mito_stats=f"{config['mito_ChIP_calculate']['output']['mito_stats']}/{{sample}}_mito_stats.txt"
        
    params:
        mito_chr=config['mito_ChIP_calculate']['params']['mito_chr']
    
    resources:
        mem_mb=config['mito_ChIP_calculate']['resources']['mem_mb'],
        time=config['mito_ChIP_calculate']['resources']['time']
            
    benchmark: "benchmarks/mito_ChIP_calculate/{sample}.txt"        
    log: "logs/mito_ChIP_calculate/{sample}.err"        
    conda: "envs/samtools.yaml"
    container: "https://depot.galaxyproject.org/singularity/samtools:1.21--h96c455f_1"
    threads: config['mito_ChIP_calculate']['threads']
        
    message:
        "[MITOCHONDRIAL READS] SAMPLES: {wildcards.sample}|INPUT: {input.sorted_bam}|OUTPUT: {output.mito_stats}|PATTERN: {params.mito_chr}"
        
    shell:
        """
        #Index BAM if not already indexed
            if [ ! -f {input.sorted_bam}.bai ]; then
                samtools index {input.sorted_bam}
            fi
        
            # Total mapped reads (excluding unmapped)
            total=$(samtools view -c -F 4 {input.sorted_bam})
        
            # Mitochondrial reads
            mito=0
            if [ "${{total}}" -ne 0 ]; then
                mito=$(samtools view -c {input.sorted_bam} {params.mito_chr})
            fi
        
            # Calculate fraction
            fraction=0
            if [ "${{total}}" -gt 0 ]; then
                fraction=$(echo "scale=6; ${{mito}} / ${{total}}" | bc -l)
            fi
        
            echo "Total Reads: ${{total}}" > {output.mito_stats}
            echo "Mito Reads: ${{mito}}" >> {output.mito_stats}
            echo "Mito Fraction: ${{fraction}}" >> {output.mito_stats}

          """
