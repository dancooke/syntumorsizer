import random
import csv
import pysam as ps

def sum_region_sizes(bed_path):
    with open(bed_path) as bed:
        bedreader = csv.reader(bed, delimiter="\t")
        return sum(int(rec[2]) - int(rec[1]) for rec in bedreader)

def count_records(vcf_path):
    vcf = ps.VariantFile(vcf_path)
    return sum(1 for _ in vcf)

def get_germline_copy_numbers(chrom, sex):
    if "Y" in chrom:
        return (1, 0) if sex == "male" else (0, 0)
    elif "X" in chrom:
        return (1, 0) if sex == "male" else (1, 1)
    else:
        return (1, 1)

def choose_clone(clone_probs):
    return random.choices(range(len(clone_probs)), weights=clone_probs)[0]

def choose_spike_haplotype(maternal_cn, paternal_cn):
    maternal_prob = float(maternal_cn) / (maternal_cn + paternal_cn)
    return 0 if random.random() < maternal_prob else 1

def simulate_tumour(source_vcf, output_vcf, mutation_rates, calling_bed, tumour_purity, subclone_ccfs, sex):
    assert len(mutation_rates) == len(subclone_ccfs) + 1
    assert sex in ["male", "female"]

    callable_genome_size = sum_region_sizes(calling_bed)
    
    in_vcf = ps.VariantFile(source_vcf)
    in_vcf.header.add_meta("INFO", items=[("ID","CLONE"), ("Number",1), ("Type","Integer"), ("Description","Spike-in clone")])
    in_vcf.header.add_meta("INFO", items=[("ID","HAP"), ("Number",1), ("Type","Integer"), ("Description","Spike-in haplotype")])
    in_vcf.header.add_meta("INFO", items=[("ID","VAF"), ("Number",1), ("Type","Float"), ("Description","Somatic mutation VAF")])
    in_vcf.header.add_meta("INFO", items=[("ID","CN"), ("Number",1), ("Type","Integer"), ("Description","Total copy number")])
    in_vcf.header.add_meta("INFO", items=[("ID","HAP_CN"), ("Number",1), ("Type","Integer"), ("Description","Copy number of the spike-in in haplotype")])

    num_source_variants = count_records(source_vcf)
    assert num_source_variants > 0
    selection_rate = sum(mutation_rates) * callable_genome_size / num_source_variants
    if selection_rate >= 1.0:
        selection_rate = None
    clone_probs = [r / sum(mutation_rates) for r in mutation_rates]
    ccfs = [1.0] + subclone_ccfs

    with ps.VariantFile(output_vcf, "wz", header=in_vcf.header) as out_vcf:
        for rec in in_vcf:
            if selection_rate is None or random.random() < selection_rate:
                clone = choose_clone(clone_probs)
                germline_copy_numbers = get_germline_copy_numbers(rec.chrom, sex)
                ploidy = sum(germline_copy_numbers)
                hap = choose_spike_haplotype(*germline_copy_numbers)
                hap_cn = germline_copy_numbers[hap]
                allele_cn = random.randint(1, hap_cn)
                vaf = ccfs[clone] * tumour_purity * allele_cn / ploidy
                rec.info["CLONE"] = clone
                rec.info["VAF"] = vaf
                rec.info["HAP"] = hap
                rec.info["CN"] = ploidy
                rec.info["HAP_CN"] = hap_cn
                out_vcf.write(rec)


random.seed(snakemake.params.seed)

simulate_tumour(snakemake.input.vcf, \
                snakemake.output[0], \
                snakemake.params.mutation_rates, \
                snakemake.params.callable_bed if snakemake.params.callable_bed is not None else snakemake.input.bed, \
                snakemake.params.tumour_purity, \
                snakemake.params.subclone_ccfs, \
                snakemake.params.sex)
