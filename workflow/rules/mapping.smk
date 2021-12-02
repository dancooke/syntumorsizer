rule samtools_index_fasta:
    input:
        "{fasta}"
    output:
         "{fasta}.fai"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input}"

localrules: samtools_index_fasta

rule bwa_index:
    input:
        fa="data/references/{reference}.fa",
        fai="data/references/{reference}.fa.fai"
    output:
        "data/references/{reference}.fa.amb",
        "data/references/{reference}.fa.ann",
        "data/references/{reference}.fa.bwt",
        "data/references/{reference}.fa.pac",
        "data/references/{reference}.fa.sa"
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa index {input.fa}"

rule bwa_map:
    input:
        rules.bwa_index.output,
        fa="data/references/{reference}.fa",
        fai="data/references/{reference}.fa.fai",
        fq1="data/reads/raw/{kind}/{sample}_R1.fastq.gz",
        fq2="data/reads/raw/{kind}/{sample}_R2.fastq.gz"
    output:
        "data/reads/mapped/{kind}/{sample}.{reference}.bwa-mem.bam"
    params:
        rg=r"@RG\tID:0\tSM:{sample}\tLB:syntumour\tPU:illumina",
        sort_threads=4,
        sort_memory_per_thread="4G"
    log:
        "logs/bwa/{kind}/{sample}.{reference}.log"
    threads: int(config["threads"])
    conda:
        "../envs/bwa.yaml"
    shell:
        "(bwa mem -t {threads} -R '{params.rg}' {input.fa} {input.fq1} {input.fq2} | \
          samtools view -bh | \
          samtools sort -@ {params.sort_threads} -m {params.sort_memory_per_thread} -o {output}) \
         2> {log}"
