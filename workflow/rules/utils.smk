rule samtools_index_bam:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input}"

rule tabix_vcf:
    input:
        "{prefix}.vcf.gz"
    output:
        "{prefix}.vcf.gz.tbi"
    conda:
        "../envs/samtools.yaml"
    shell:
        "tabix {input}"

rule sort_bam:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.sorted.bam"
    conda:
        "../envs/samtools.yaml"
    threads: int(config["threads"])
    shell:
        "samtools sort -@{threads} -o {output} {input}"

rule sort_and_index_bam:
    input:
        "{prefix}.bam"
    output:
        bam="{prefix}.sorted.bam",
        bai="{prefix}.sorted.bam.bai"
    params:
        bai="{prefix}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    threads: int(config["threads"])
    shell:
        """
        if [ -f "{params.bai}" ]; then
            ln -sr {input} {output.bam}
            ln -sr {params.bai} {output.bai}
        else
            samtools sort -@{threads} -o {output.bam} {input}
            samtools index -@{threads} {output.bam}
        fi
        """

ruleorder: sort_and_index_bam > sort_bam > samtools_index_bam

rule bam_to_fastqs:
    input:
        f"data/reads/mapped/somatic/merged/{{sample}}.{config['reference']}.{config['mapper']}.{config['caller']}.bam"
    output:
        expand("data/reads/raw/somatic/{{sample}}_{end}.fastq.gz",
               end=["R1", "R2", "single", "unpaired"])
    conda:
        "../envs/samtools.yaml"
    threads: int(config["threads"])
    shell:
        "samtools sort -n -@{threads} {input} | \
          samtools fastq -N -@{threads} \
            -1 {output[0]} -2 {output[1]} \
            -s {output[2]} -0 {output[3]}"

rule complete_fastqs:
    input:
        somatic_fq1="data/reads/raw/somatic/{sample}_{tumour}_R1.fastq.gz",
        somatic_fq2="data/reads/raw/somatic/{sample}_{tumour}_R2.fastq.gz",
        somatic_single="data/reads/raw/somatic/{sample}_{tumour}_single.fastq.gz",
        somatic_unpaired="data/reads/raw/somatic/{sample}_{tumour}_unpaired.fastq.gz",
        germline_fq1="data/reads/raw/germline/{sample}_R1.fastq.gz",
        germline_fq2="data/reads/raw/germline/{sample}_R2.fastq.gz"
    output:
        fq1="results/{sample}_{tumour}_R1.fastq.gz",
        fq2="results/{sample}_{tumour}_R2.fastq.gz"
    conda:
        "../envs/py39.yaml"
    script:
        "../scripts/complete_fastqs.py"

rule link_variants:
    input:
        "data/variants/somatic/{sample}.{reference}.vcf.gz",
        "data/variants/somatic/{sample}.{reference}.vcf.gz.tbi"
    output:
        "results/{sample}.{reference}.vcf.gz",
        "results/{sample}.{reference}.vcf.gz.tbi"
    shell:
        """
        ln -sr {input[0]} {output[0]}
        ln -sr {input[1]} {output[1]}
        """

localrules: link_variants
ruleorder: link_variants > tabix_vcf

rule link_bam:
    input:
        "data/reads/mapped/somatic/{sample}.{reference}.{mapper}.bam",
        "data/reads/mapped/somatic/{sample}.{reference}.{mapper}.bam.bai"
    output:
        "results/{sample}.{reference}.{mapper}.bam",
        "results/{sample}.{reference}.{mapper}.bam.bai"
    shell:
        """
        ln -sr {input[0]} {output[0]}
        ln -sr {input[1]} {output[1]}
        """

localrules: link_bam
ruleorder: link_bam > samtools_index_bam
