# Manage Genomes
This snakemake workflow is written to perform multiple steps to prepare the genomes for further downstream analysis in various other workflows.

The steps performed in this workflow include:

1. Download genome from Ensembl (or Refseq: under development)
2. Download gtf from Ensembl (or Refseq: under development)
3. Create indices for bowtie1, bowtie2, bwa, star, cellran ger, cell ranger vdj
4. Generate various derivatives of gtf: bed, refflat, igv, and gene-transcript relationship file

