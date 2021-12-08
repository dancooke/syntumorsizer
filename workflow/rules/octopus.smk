def _get_contig_ploidies(wildcards, input):
    if config["sex"] == "male":
        return ["X=1", "Y=1", "chrX=1", "chrY=1"]
    else:
        return ["X=2", "Y=0", "chrX=2", "chrY=0"]

def _get_calling_bed(wildcards):
    if "regions" in config:
        return config["regions"]
    else:
        return f"data/references/{wildcards.reference}.chromosomes.bed"

rule octopus_call:
    input:
        reference="data/references/{reference}.fa",
        bam="data/reads/mapped/germline/{sample}.{reference}.{mapper}.bam",
        bai="data/reads/mapped/germline/{sample}.{reference}.{mapper}.bam.bai",
        bed=_get_calling_bed
    output:
        vcf="data/variants/germline/{sample}.{reference}.{mapper}.octopus.vcf.gz",
        vcf_index="data/variants/germline/{sample}.{reference}.{mapper}.octopus.vcf.gz.tbi"
    params:
        contig_ploidies=_get_contig_ploidies,
        forest="/opt/octopus/resources/forests/germline.v0.7.4.forest.gz",
        other=config["octopus"] if "octopus" in config else ""
    log:
        "logs/octopus/{sample}.{reference}.{mapper}.log"
    threads: int(config["threads"])
    container:
        "docker://dancooke/octopus:develop"
    shell:
        "(octopus \
         -R {input.reference} \
         -I {input.bam} \
         -t {input.bed} \
         -o {output} \
         --contig-ploidies {params.contig_ploidies} \
         --forest {params.forest} \
         --threads {threads} \
         {params.other} \
         )2> {log}"

rule octopus_realign:
    input:
        reference="data/references/{reference}.fa",
        bam="data/reads/mapped/germline/{sample}.{reference}.{mapper}.bam",
        bai="data/reads/mapped/germline/{sample}.{reference}.{mapper}.bam.bai",
        bed=_get_calling_bed,
        source_vcf=rules.octopus_call.output.vcf
    output:
        vcf="data/variants/germline/{sample}.{reference}.{mapper}.octopus-realign.vcf.gz",
        vcf_index="data/variants/germline/{sample}.{reference}.{mapper}.octopus-realign.vcf.gz.tbi",
        bam="data/reads/mapped/germline/realigned/haplotagged/{sample}.{reference}.{mapper}.octopus.bam"
    params:
        contig_ploidies=_get_contig_ploidies
    log:
        "logs/octopus/{sample}.{reference}.{mapper}.realign.log"
    threads: int(config["threads"])
    container:
        "docker://dancooke/octopus:develop"
    shell:
        "(octopus \
         -R {input.reference} \
         -I {input.bam} \
         -t {input.bed} \
         -o {output} \
         --contig-ploidies {params.contig_ploidies} \
         --threads {threads} \
         --disable-denovo-variant-discovery \
         --source-candidates {input.source_vcf} \
         --max-haplotypes 500 \
         --lagging-level OPTIMISTIC \
         --backtrack-level AGGRESSIVE \
         --min-protected-haplotype-posterior 1e-5 \
         --bamout {output.bam} \
         --bamout-type FULL \
         )2> {log}"

ruleorder: octopus_call > tabix_vcf
ruleorder: octopus_realign > tabix_vcf
