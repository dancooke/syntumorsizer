rule generate_somatic_mutations:
    input:
        vcf=config["source_vcf"],
        bed="data/references/{reference}.chromosomes.bed"
    output:
        "data/variants/somatic/{sample}_{tumour}.{reference}.vcf.gz"
    params:
        sex=config["sex"],
        tumour_purity=float(config["tumour_purity"]),
        subclone_ccfs=list(map(float, config["subclone_ccfs"])),
        mutation_rates=list(map(float, config["mutation_rates"])),
        callable_bed=config["regions"] if "regions" in config else None,
        seed=int(config["seed"])
    conda:
        "../envs/py39.yaml"
    script:
        "../scripts/simulate_tumour.py"
