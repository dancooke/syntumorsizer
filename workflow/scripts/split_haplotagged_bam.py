from pathlib import Path
import numpy as np
import pysam as ps

def get_haplotype_ids(read):
    try:
        return tuple([int(id) for id in read.get_tag('HP').split(',')])
    except KeyError:
        return None

def split_bam(in_bam_path, out_bam_paths):
    in_bam = ps.AlignmentFile(in_bam_path)
    out_bams = [ps.AlignmentFile(bam, 'wb', template=in_bam) for bam in out_bam_paths]

    for read in in_bam:
        haplotype_ids = get_haplotype_ids(read)
        if haplotype_ids is None:
            haplotype = np.random.choice(len(out_bams))
        else:
            haplotype = np.random.choice(haplotype_ids)
        out_bams[haplotype].write(read)

split_bam(Path(snakemake.input[0]), [Path(bam) for bam in snakemake.output])
