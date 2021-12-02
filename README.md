# Syntumorsizer

Synumorsizer is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow for making synthetic tumour sequencing data. It implements the procedure described in [Cooke et al.](https://www.nature.com/articles/s41587-021-00861-3) that haplotypes and realigns germline reads before somatic mutation spike-in, resulting in more realistic haplotype structure than if naive spike-in is used.

The basic input of the workflow is:

* Germline sequencing reads of the the sample to be mutated, in FASTQ format.
* Somatic variants to spike into the germline sequencing reads, in VCF.

and the output is:

* Sequencing reads with somatic mutations spiked in, in FASTQ or BAM format.

## Running

To run the workflow, you'll need [Snakemake](https://snakemake.readthedocs.io/en/stable/), [Conda](https://docs.conda.io/en/latest/), and [Singularity](https://sylabs.io/guides/latest/user-guide/) installed.

### Configuration

You'll need to complete a YAML file specifying the inputs and parameters. An example is included in the  `config` directory.

### Example

To run the bundled example locally:

```shell
$ snakemake --configfile config/example.yaml --cores 16 --use-conda --use-singularity
```

Once the workfow is complete, you should find a directory `results` including:

* `NA24385_COAD_R{1,2}.fastq.gz`: raw spiked sequencing reads.
* `NA24385_COAD.hs37d5.bwa-mem.bam`: remapped spiked sequencing reads.
* `NA24385_COAD.hs37d5.vcf.gz`: spiked-in somatic mutations - a subset of the input variants with added annotations.

Note that the input variants will be sampled to satisfy the specified somatic mutation rate(s), so if you want all the input variants to be included then just make sure the `mutation_rates` are sufficiently high.

Take a look at the Snakemake docs to learn other ways of executing the workflow, such as running on a [cluster](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).

## Limitations

* Only paired-end Illuina quality reads are supported.
* Only [BWA-MEM](https://github.com/lh3/bwa) and [Octopus](https://github.com/luntergroup/octopus) are options for read mapping and haplotyping, respectively.
* Tumour subclonal structure (aka tumour phylogeny) is not modelled (only local haplotype structure is respected).
* No CNVs - only small variants (SNVs/MNVs/indels).
* Human samples only.
