rule deeptools_bamCompare:
    input:
        treatment_bam=lambda wildcards: f"{config['bamCompare']['input']['filtered_bam']}/{wildcards.sample}.filtered.bam",
        control_bam=lambda wildcards: f"{config['bamCompare']['input']['filtered_bam']}/{CONTROLS[wildcards.sample]}.filtered.bam"
        
    output:
        normalized_bw=f"{config['bamCompare']['output']['normalized_coverage']}/{{sample}}_vs_control.bw"
        
    params:
        operation=config['bamCompare']['params']['operation'],
        binSize=config['bamCompare']['params']['binSize'],
        format=config['bamCompare']['params']['outFileFormat']
        
    resources:
        mem_mb=config['bamCompare']['resources']['mem_mb'],
        time=config['bamCompare']['resources']['time']
        
    benchmark: f"benchmarks/bamCompare/{{sample}}.txt"
    log: f"logs/bamCompare/{{sample}}.log"
    conda: "envs/deeptools.yaml"
    container: "https://depot.galaxyproject.org/singularity/deeptools:3.5.5--pyhdfd78af_0"
    threads: config['bamCompare']['threads']

    message: "[bamCompare] Normalizing {wildcards.sample} against {input.control_bam} | Operation: {params.operation}"
    
    shell:
        """
        bamCompare \
            -b1 {input.treatment_bam} \
            -b2 {input.control_bam} \
            -o {output.normalized_bw} \
            --operation {params.operation} \
            --binSize {params.binSize} \
            --outFileFormat {params.format} \
            --numberOfProcessors {threads} \
            2> {log}
        """
