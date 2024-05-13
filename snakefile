# Download the genome (fasta) and gtf files from the Ensembl ftp site (https://ftp.ensembl.org/pub/) 
# For assembly GRch38, works for release >= 47 (all available species)
# For assembly GRch37, works for release >= 83 (only homo_sapiens is availbale)
# Builds genome indices for various mapping tools such as STAR, bwa, bowtie, ...
# Convert gtf to refflat, bed, igv and other formats


import os
import subprocess
from os.path import exists, abspath
from bs4 import BeautifulSoup
#import urllib.request
import re
import random
import shutil
from collections import Counter


# config file
configfile: "config.yaml"

# Set seed for reproducibility
random.seed(config["seed"])

# Base url for fetching data
source = config["source"]
base_url = config[source]["base_url"]


# Import all scripts
# Define a function to filter script files
def get_script_files(directory):
    script_files = []
    for filename in os.listdir(directory):
        if filename.endswith((".py", ".R", ".pl")):     # ".sh",
            script_files.append(os.path.join(directory, filename))
    return script_files

# Fetch all script files from the scripts directory
script_files = get_script_files("scripts/")

# Include all script files in the workflow
for script in script_files:
    include: script


###########################################
## Parameters from the config and checks ##
###########################################

# resources local base path ############

## Downloads
if(config[source]["fasta"]["download"]):
    fasta_download_path = path_checker(name = "FASTA download", specific_path = config[source]["fasta"]["local_path"], general_path = config["local"]["path_genome"])
else:
    fasta_download_path = os.getcwd()

if(config[source]["annotation"]["download"]):
    annotation_download_path = path_checker(name = "ANNOTATION download", specific_path = config[source]["annotation"]["local_path"], general_path = config["local"]["path_genome"])
else:
    annotation_download_path = os.getcwd()

if(config[source]["repeats"]["download"]):
    ensembl_repeats_download_path = path_checker(name = "ENSEMBL REPEATS download", specific_path = config[source]["annotation"]["local_path"], general_path = config["local"]["path_genome"])
else:
    ensembl_repeats_download_path = os.getcwd()

## Indices
### bowtie1
if(config["build_indices"]["bowtie1"]["run"]):
    bowtie1_indices_path = path_checker(name = "BOWTIE1 index", specific_path = config["build_indices"]["bowtie1"]["local_path"], general_path = config["local"]["path_indices"])
else:
    bowtie1_indices_path = os.getcwd()

bowtie1_indices_path = os.path.join(bowtie1_indices_path, "indices", "bowtie1")

### bowtie2
if(config["build_indices"]["bowtie2"]["run"]):
    bowtie2_indices_path = path_checker(name = "BOWTIE2 index", specific_path = config["build_indices"]["bowtie2"]["local_path"], general_path = config["local"]["path_indices"])
else:
    bowtie2_indices_path = os.getcwd()

bowtie2_indices_path = os.path.join(bowtie2_indices_path, "indices", "bowtie2")

### bwa
if(config["build_indices"]["bwa"]["run"]):
    bwa_indices_path = path_checker(name = "BWA index", specific_path = config["build_indices"]["bwa"]["local_path"], general_path = config["local"]["path_indices"])
else:
    bwa_indices_path = os.getcwd()

bwa_indices_path = os.path.join(bwa_indices_path, "indices", "bwa")

### cellranger
if(config["build_indices"]["cellranger"]["run"]):
    cellranger_indices_path = path_checker(name = "CELLRANGER index", specific_path = config["build_indices"]["cellranger"]["local_path"], general_path = config["local"]["path_indices"])
else:
    cellranger_indices_path = os.getcwd()

cellranger_indices_path = os.path.join(cellranger_indices_path, "indices", "cellranger")

### kallisto
if(config["build_indices"]["kallisto"]["run"]):
    kallisto_indices_path = path_checker(name = "KALLISTO index", specific_path = config["build_indices"]["kallisto"]["local_path"], general_path = config["local"]["path_indices"])
else:
    kallisto_indices_path = os.getcwd()

kallisto_indices_path = os.path.join(kallisto_indices_path, "indices", "kallisto")

### star
if(config["build_indices"]["star"]["run"]):
    star_indices_path = path_checker(name = "STAR index", specific_path = config["build_indices"]["star"]["local_path"], general_path = config["local"]["path_indices"])
else:
    star_indices_path = os.getcwd()

star_indices_path = os.path.join(star_indices_path, "indices", "star")

### xengsort
if(config["build_indices"]["xengsort"]["run"]):
    xengsort_indices_path = path_checker(name = "XENGSORT index", specific_path = config["build_indices"]["xengsort"]["local_path"], general_path = config["local"]["path_indices"])
else:
    xengsort_indices_path = os.getcwd()

xengsort_indices_path = os.path.join(xengsort_indices_path, "indices", "xengsort")


## GTF derivatives


# Check what combinations exist ##############

if source == "ensembl":
    working_dict = combination_checks_ensembl(config)
elif source == "refseq":
    working_dict = combination_checks_refseq(config)


# Check what files are available for the valid combinations and patterns ##############

if source == "ensembl":
    downloadable_files = downloadable_ensembl(config, working_dict)
elif source == "refseq":
    downloadable_files = downloadable_refseq(config, working_dict)

print(downloadable_files)
# targets ############

targets = []
files_to_decompress_fa = []
files_to_decompress_anno = []

for key, values in downloadable_files.items():
    
    assembly = values["assembly"]
    release = values["release"]
    species = values["species"]
    seqtype = values["seqtype"]
    genome_files = values["fasta"]
    anno_files = values["annotation"]
        
    if(config["ensembl"]["fasta"]["download"]):
        targets.append(os.path.join(fasta_download_path, species, f"{assembly}_{release}_{seqtype}", "genome.fa"))

    
    if(config["ensembl"]["annotation"]["download"]):
        targets.append(os.path.join(fasta_download_path, species, f"{assembly}_{release}_annotation", "gtf.gtf"))


    if(source == "ensembl" and 
        config["build_indices"]["bowtie1"]["run"]):
        targets.append(os.path.join(bowtie1_indices_path, species, f"{assembly}_{release}_{seqtype}", f"{species}.4.ebwt"))
        targets.append(os.path.join(bowtie1_indices_path, species, f"{assembly}_{release}_{seqtype}", f"{species}.rev.2.ebwt"))


    if(source == "ensembl" and 
        config["build_indices"]["bowtie2"]["run"]):
        targets.append(os.path.join(bowtie2_indices_path, species, f"{assembly}_{release}_{seqtype}", f"{species}.4.bt2"))
        targets.append(os.path.join(bowtie2_indices_path, species, f"{assembly}_{release}_{seqtype}", f"{species}.rev.2.bt2"))


    if(source == "ensembl" and 
        config["build_indices"]["bwa"]["run"]):
        targets.append(os.path.join(bwa_indices_path, f"{species}", f"{assembly}_{release}_{seqtype}", f"{species}.sa"))


    if(source == "ensembl" and 
        config["build_indices"]["star"]["run"]):
        targets.append(os.path.join(star_indices_path, f"{species}", f"{assembly}_{release}_{seqtype}", "SA"))


    if(source == "ensembl" and 
        config["build_indices"]["kallisto"]["run"]):
        targets.append(os.path.join(kallisto_indices_path, f"{species}", f"{assembly}_{release}_{seqtype}", f"{species}.idx"))

    
    if(source == "ensembl" and 
        config["build_indices"]["cellranger"]["run"]):
        targets.append(os.path.join(cellranger_indices_path, f"{species}", f"{assembly}_{release}_{seqtype}", "star", "SA"))
    
    
    if(config["annotation_files"]["gene_transcript_relation"]["create"]):
        targets.append(os.path.join(annotation_download_path, f"{species}", f"{assembly}_{release}_annotation", "gene-transcript.txt"))

if(source == "ensembl" and 
    config["build_indices"]["xengsort"]["run"]):
    
    if(len(downloadable_files) != 2):
        error(f"For xengsort index building, two genomes are required. Found {len(downloadable_files)}")
    else:
        species_list = []
        pathname_list = []

        for key in downloadable_files:
            species_list.append(downloadable_files[key]['species'])
            
            #s_parts = downloadable_files[key]['species'].split("_")
            #pathname_list.append("".join(word[0] for word in s_parts))

            pathname_list.append(downloadable_files[key]['species'])
            pathname_list.append(downloadable_files[key]['assembly'])

        species_list = [downloadable_files[key]['species'] for key in downloadable_files]
        species_counts = Counter(species_list)
        
        if(len(species_counts) != 2):
            error(f"For xengsort index building, two species are required. Found {len(species_counts)}")
        else:
            pathname = '_'.join(pathname_list)
            targets.append(os.path.join(xengsort_indices_path, pathname, "xengsort"))


targets = list(set(targets))
print('\n\n'.join(map(str, targets)) + '\n\n')


# rules ############
rule all:
    input: targets


# Import all rule files from the rules directory ##########

# Downloads #
include: "rules/download.smk"


# Indices #
include: "rules/bowtie1_index.smk"
include: "rules/bowtie2_index.smk"
include: "rules/bwa_index.smk"
include: "rules/cellranger_index.smk"
#include: "rules/cellranger_vdj_index.smk"
include: "rules/kallisto_index.smk"
include: "rules/star_index.smk"
#include: "rules/xengsort_index.smk"


# GTF derivatives #
#include: "rules/bed.smk"
include: "rules/gene_transcript.smk"
#include: "rules/igv.smk"
#include: "rules/refflat.smk"
