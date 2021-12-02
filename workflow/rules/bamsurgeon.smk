rule split_bam:
    input:
        "data/reads/mapped/germline/realigned/haplotagged/{sample}.{reference}.{mapper}.{caller}.bam"
    output:
        expand("data/reads/mapped/germline/realigned/split/{{sample}}.{{reference}}.{{mapper}}.{{caller}}.hap_{haplotype}.bam", 
               haplotype=[0, 1])
    conda:
        "../envs/py39.yaml"
    script:
        "../scripts/split_haplotagged_bam.py"

rule make_bamsurgeon_bed:
    input:
        "data/variants/somatic/{sample}.{reference}.vcf.gz"
    output:
        "data/variants/somatic/bamsurgeon/{sample}.{reference}.hap_{haplotype}.{kind}.bed"
    conda:
        "../envs/py39.yaml"
    script:
        "../scripts/vcf_to_bamsurgeon_bed.py"

rule bamsurgeon_addsnv:
    input:
        reference="data/references/{reference}.fa",
        bam="data/reads/mapped/germline/realigned/split/{sample}.{reference}.{mapper}.{caller}.hap_{haplotype}.sorted.bam",
        bai="data/reads/mapped/germline/realigned/split/{sample}.{reference}.{mapper}.{caller}.hap_{haplotype}.sorted.bam.bai",
        snv_bed="data/variants/somatic/bamsurgeon/{sample}_{tumour}.{reference}.hap_{haplotype}.snv.bed"
    output:
        bam="data/reads/mapped/somatic/split/{sample}_{tumour}.{reference}.{mapper}.{caller}.hap_{haplotype}.snv.bam",
        vcf=temp("{sample}_{tumour}.{reference}.{mapper}.{caller}.hap_{haplotype}.snv.addsnv.{sample}_{tumour}.{reference}.hap_{haplotype}.snv.vcf"),
        logs=temp(directory("addsnv_logs_{sample}_{tumour}.{reference}.{mapper}.{caller}.hap_{haplotype}.snv.bam"))
    wildcard_constraints:
        haplotype="[0-9]+"
    params:
        path_hack="/opt/conda/envs/bamsurgeon/bin:/opt/bamsurgeon/bin",
        tmp_dir="data/reads/mapped/somatic/split/{sample}_{tumour}.{reference}.{mapper}.{caller}.hap_{haplotype}.addsnv",
        seed=int(config["seed"])
    threads: int(config["threads"])
    log:
        "logs/bamsurgeon/addsnv/{sample}_{tumour}.{reference}.{mapper}.{caller}.{haplotype}.log"
    container:
         "docker://dancooke/bamsurgeon"
    shell:
        """
        export PATH={params.path_hack}:$PATH
        (addsnv.py \
            -r {input.reference} \
            -f {input.bam} \
            -v {input.snv_bed} \
            -o {output.bam} \
            -p {threads} \
            -m 1 \
            --mindepth 1 \
            --force \
            --ignoresnps \
            --ignorepileup \
            --aligner ssw \
            --haplosize 200 \
            --seed {params.seed} \
            --tmpdir {params.tmp_dir} \
        )&> {log}
        ln -sr {input.bam} {output.bam} || true
        rm -rf {params.tmp_dir}
        touch {output.vcf}
        """

rule bamsurgeon_addindel:
    input:
        reference="data/references/{reference}.fa",
        bam="data/reads/mapped/somatic/split/{sample}.{reference}.{mapper}.{caller}.hap_{haplotype}.snv.sorted.bam",
        bai="data/reads/mapped/somatic/split/{sample}.{reference}.{mapper}.{caller}.hap_{haplotype}.snv.sorted.bam.bai",
        indel_bed="data/variants/somatic/bamsurgeon/{sample}.{reference}.hap_{haplotype}.indel.bed"
    output:
        bam="data/reads/mapped/somatic/split/{sample}.{reference}.{mapper}.{caller}.hap_{haplotype}.bam",
        vcf=temp("{sample}.{reference}.{mapper}.{caller}.hap_{haplotype}.addindel.{sample}.{reference}.hap_{haplotype}.indel.vcf"),
        logs=temp(directory("addindel_logs_{sample}.{reference}.{mapper}.{caller}.hap_{haplotype}.bam"))
    wildcard_constraints:
        haplotype="[0-9]+"
    params:
        path_hack="/opt/conda/envs/bamsurgeon/bin:/opt/bamsurgeon/bin",
        tmp_dir="data/reads/mapped/somatic/split/{sample}.{reference}.{mapper}.{caller}.hap_{haplotype}.addindel",
        seed=int(config["seed"])
    threads: int(config["threads"])
    log:
        "logs/bamsurgeon/addindel/{sample}.{reference}.{mapper}.{caller}.hap_{haplotype}.log"
    container:
         "docker://dancooke/bamsurgeon"
    shell:
        """
        export PATH={params.path_hack}:$PATH
        (addindel.py \
            -r {input.reference} \
            -f {input.bam} \
            -v {input.indel_bed} \
            -o {output.bam} \
            -p {threads} \
            -m 1 \
            --mindepth 1 \
            --force \
            --ignorepileup \
            --aligner ssw \
            --seed {params.seed} \
            --tmpdir {params.tmp_dir} \
        )&> {log}
        ln -sr {input.bam} {output.bam} || true
        rm -rf {params.tmp_dir}
        touch {output.vcf}
        """

rule merge_bams:
    input:
        expand("data/reads/mapped/somatic/split/{{sample}}.{{reference}}.{{mapper}}.{{caller}}.hap_{haplotype}.bam",
               haplotype=[0, 1])
    output:
        "data/reads/mapped/somatic/merged/{sample}.{reference}.{mapper}.{caller}.bam"
    conda:
        "../envs/samtools.yaml"
    threads: int(config["threads"])
    shell:
        """
        samtools merge -@{threads} -o {output} {input}
        rm -f add*.*.muts.bam
        """ # Cleaning up files here from bamsurgen as names are random!
