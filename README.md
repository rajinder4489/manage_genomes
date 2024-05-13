# Manage Genomes
This snakemake workflow is written to perform multiple steps to prepare the genomes for further downstream analysis in various other workflows.

The steps performed in this workflow include:

1. Download genome from Ensembl (or Refseq: under development)
2. Download gtf from Ensembl (or Refseq: under development)
3. Create indices for bowtie1, bowtie2, bwa, star, cellranger, cellranger vdj
4. Generate various derivatives of gtf: bed, refflat, igv, and gene-transcript relationship file

To do:

1. Handle multiple options for assembly, species, release and type : Done
2. Split all rules and scripts in individual files: Done
3. Build indices for xengsort for multiple species together: 
4. Get repeats: from ensembl under mysql folder
5. vep files
6. Download data from Refseq or other sources
7. ERCC: may be add this manually, if too much of work
8. Restructure paths for download and file creation: Done


Refer this path for info /dcgc/support/pipeline/fastqc 