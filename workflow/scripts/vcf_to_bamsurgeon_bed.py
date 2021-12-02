import csv
import pysam as ps

def is_ins(rec):
    return len(rec.alts[0]) > len(rec.ref)

def is_del(rec):
    return len(rec.alts[0]) < len(rec.ref)

def is_indel(rec):
    return is_ins(rec) or is_del(rec)

def make_bamurgeon_record(rec, cn):
    result = [rec.chrom]
    vaf = cn * rec.info["VAF"]
    if is_indel(rec):
        if is_del(rec):
            result += [rec.pos, rec.stop]
            result += [vaf, 'DEL']
        else:
            result += [rec.stop, rec.stop]
            result += [vaf, 'INS', rec.alts[0]]
    else:
        result += [rec.pos, rec.pos]
        result += [vaf, rec.alts[0]]
    return result

def vcf_to_bamsurgeon_bed(vcf_path, bed_path, kind, haplotype):
    def _get_kind(rec):
        return "indel" if is_indel(rec) else "snv"

    with ps.VariantFile(vcf_path) as vcf, open(bed_path, "w") as bed:
        bedwriter = csv.writer(bed, delimiter="\t")
        for rec in vcf:
            if _get_kind(rec) == kind:
                hap, cn = rec.info["HAP"], rec.info["CN"]
                if hap == haplotype:
                    bs_rec = make_bamurgeon_record(rec, cn)
                    bedwriter.writerow(bs_rec)

vcf_to_bamsurgeon_bed(snakemake.input[0], snakemake.output[0], \
                      snakemake.wildcards.kind, \
                      int(snakemake.wildcards.haplotype))
