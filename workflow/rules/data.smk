rule bwakit_genref:
    output:
        "data/references/{name}.fa"
    conda:
        "../envs/bwakit.yaml"
    shell:
        """
        run-gen-ref {wildcards.name}
        mv {wildcards.name}.fa* $(dirname {output[0]})
        """

rule make_chromosomes_bed:
    input:
        "data/references/{reference}.fa.fai"
    output:
        "data/references/{reference}.chromosomes.bed"
    shell:
        "head -25 {input} | awk -v OFS='\t' '{{print $1,0,$2}}' > {output}"

rule link_germline_fastq:
    input:
        lambda wildcards: config["germline_fastqs"][0 if int(wildcards.strand) == 1 else 1]
    output:
        "data/reads/raw/germline/{sample}_R{strand}.fastq.gz"
    shell:
        "ln -sr {input} {output}"
