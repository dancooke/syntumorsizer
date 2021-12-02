import gzip
from collections import defaultdict
from pathlib import Path
import pysam as ps

def get_read_number(read):
    return int(read.name[-1])

def complete_fastqs(somatic_fq1_path, somatic_fq2_path, somatic_single_fq_path, somatic_unpaired_fq_path, \
                    germline_fq1_path, germline_fq2_path, \
                    output_fq1_path, output_fq2_path):
    singleton_reads = defaultdict(list), defaultdict(list)
    with ps.FastxFile(somatic_single_fq_path) as somatic_single_fq:
        for read in somatic_single_fq:
            read_num = get_read_number(read)
            if read_num is not None:
                read.name = read.name[-2] # remove pair number
                singleton_reads[read_num - 1][read.name].append(read)

    with gzip.open(output_fq1_path, "wt") as output_fq1, gzip.open(output_fq2_path, "wt") as output_fq2:
        # copy paired somatic read, removing redundant /x tag
        with ps.FastxFile(somatic_fq1_path) as somatic_fq1, ps.FastxFile(somatic_fq2_path) as somatic_fq2:
            for read1, read2 in zip(somatic_fq1, somatic_fq2):
                read1.name = read1.name[:-2]
                output_fq1.write(str(read1) + "\n")
                read2.name = read2.name[:-2]
                output_fq2.write(str(read2) + "\n")

        # complete singltons
        with ps.FastxFile(germline_fq1_path) as germline_fq1, ps.FastxFile(germline_fq2_path) as germline_fq2:
            for read1, read2 in zip(germline_fq1, germline_fq2):
                if read1.name in singleton_reads[0]:
                    output_fq1.write(str(singleton_reads[0].get(read1.name)) + "\n")
                    output_fq2.write(str(read2) + "\n")
                elif read2.name in singleton_reads[1]:
                    output_fq1.write(str(read1) + "\n")
                    output_fq2.write(str(singleton_reads[1].get(read2.name)) + "\n")

complete_fastqs(Path(snakemake.input.somatic_fq1), Path(snakemake.input.somatic_fq2), \
                Path(snakemake.input.somatic_single), Path(snakemake.input.somatic_unpaired), \
                Path(snakemake.input.germline_fq1), Path(snakemake.input.germline_fq2), \
                Path(snakemake.output.fq1), Path(snakemake.output.fq2))
