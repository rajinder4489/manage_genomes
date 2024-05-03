# Download the genome (fasta) and gtf files from the Ensembl ftp site (https://ftp.ensembl.org/pub/) 
# For assembly GRch38, works for release >= 47 (all available species)
# For assembly GRch37, works for release >= 83 (only homo_sapiens is availbale)
# Builds genome indices for various mapping tools such as STAR, bwa, bowtie, ...
# Convert gtf to refflat, bed, igv and other formats


import os
import subprocess
from os.path import exists, abspath
from bs4 import BeautifulSoup
import urllib.request
import re
import glob
import random
import shutil


# config file
configfile: "config.yaml"

# Set seed for reproducibility
random.seed(config["seed"])

# Base url for fetching data
source = config["source"]
base_url = config[source]["base_url"]

# Build indices
indices_run = [value["run"] for value in config["build_indices"].values()]

# GTF derivatives
gtf_derivatives_create = [value["create"] for value in config["annotation_files"].values()]


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

# Downloads
if(config[source]["fasta"]["download"]):
    fasta_download_path = path_checker(name = "FASTA download", specific_path = config[source]["fasta"]["local_path"], general_path = config["local"]["path_genome"])

if(config[source]["annotation"]["download"]):
    annotation_download_path = path_checker(name = "ANNOTATION download", specific_path = config[source]["annotation"]["local_path"], general_path = config["local"]["path_genome"])

if(config[source]["ensembl_repeats"]["download"]):
    ensembl_repeats_download_path = path_checker(name = "ENSEMBL REPEATS download", specific_path = config[source]["annotation"]["local_path"], general_path = config["local"]["path_genome"])


# Indices
if(config["build_indices"]["bowtie1"]["run"]):
    bowtie1_indices_path = path_checker(name = "BOWTIE1 index", specific_path = config["build_indices"]["bowtie1"]["local_path"], general_path = config["local"]["path_indices"])
    bowtie1_indices_path = os.path.join(bowtie1_indices_path, "indices", "bowtie1")

if(config["build_indices"]["bowtie2"]["run"]):
    bowtie2_indices_path = path_checker(name = "BOWTIE2 index", specific_path = config["build_indices"]["bowtie2"]["local_path"], general_path = config["local"]["path_indices"])
    bowtie2_indices_path = os.path.join(bowtie2_indices_path, "indices", "bowtie2")

if(config["build_indices"]["bwa"]["run"]):
    bwa_indices_path = path_checker(name = "BWA index", specific_path = config["build_indices"]["bwa"]["local_path"], general_path = config["local"]["path_indices"])
    bwa_indices_path = os.path.join(bwa_indices_path, "indices", "bowtie2")

if(config["build_indices"]["cellranger"]["run"]):
    cellranger_indices_path = path_checker(name = "CELLRANGER index", specific_path = config["build_indices"]["cellranger"]["local_path"], general_path = config["local"]["path_indices"])
    cellranger_indices_path = os.path.join(cellranger_indices_path, "indices", "cellranger")

if(config["build_indices"]["cellranger_vdj"]["run"]):
    cellranger_vdj_indices_path = path_checker(name = "CELLRANGER VDJ index", specific_path = config["build_indices"]["cellranger_vdj"]["local_path"], general_path = config["local"]["path_indices"])
    cellranger_vdj_indices_path = os.path.join(cellranger_vdj_indices_path, "indices", "cellranger_vdj")

if(config["build_indices"]["kallisto"]["run"]):
    kallisto_indices_path = path_checker(name = "KALLISTO index", specific_path = config["build_indices"]["kallisto"]["local_path"], general_path = config["local"]["path_indices"])
    kallisto_indices_path = os.path.join(kallisto_indices_path, "indices", "kallisto")

if(config["build_indices"]["star"]["run"]):
    star_indices_path = path_checker(name = "STAR index", specific_path = config["build_indices"]["star"]["local_path"], general_path = config["local"]["path_indices"])
    star_indices_path = os.path.join(star_indices_path, "indices", "star")

if(config["build_indices"]["xengsort"]["run"]):
    xengsort_indices_path = path_checker(name = "XENGSORT index", specific_path = config["build_indices"]["xengsort"]["local_path"], general_path = config["local"]["path_indices"])
    xengsort_indices_path = os.path.join(xengsort_indices_path, "indices", "xengsort")


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


# targets ############

targets = []
files_to_decompress_fa = []
files_to_decompress_anno = []

if(config["ensembl"]["fasta"]["download"] or 
    config["ensembl"]["annotation"]["download"]):
    
    for key, values in downloadable_files.items():
        
        assembly = values["assembly"]
        release = values["release"]
        species = values["species"]
        seqtype = values["seqtype"]
        genome_files = values["fasta"]
        anno_files = values["annotation"]

        if(config["ensembl"]["fasta"]["download"]):
            for file in genome_files:
                if(any(indices_run)):
                    files_to_decompress_fa.append(os.path.join(fasta_download_path, species, f"{assembly}_{release}_{seqtype}", file))
                    file = os.path.splitext(file)[0]
                
                targets.append(os.path.join(fasta_download_path, species, f"{assembly}_{release}_{seqtype}", file))

        if(config["ensembl"]["annotation"]["download"]):
            for file in anno_files:
                if(any(gtf_derivatives_create)):
                    
                    files_to_decompress_anno.append(os.path.join(annotation_download_path, species, f"{assembly}_{release}_annotation", file))
                    file = os.path.splitext(file)[0]
                
                targets.append(os.path.join(annotation_download_path, species, f"{assembly}_{release}_annotation", file))


if(source == "ensembl" or 
    config["build_indices"]["bowtie1"]["run"]):
    
    for key, values in downloadable_files.items():
        
        assembly = values["assembly"]
        release = values["release"]
        species = values["species"]
        seqtype = values["seqtype"]
        genome_files = values["fasta"]
        anno_files = values["annotation"]


targets = list(set(targets))
#print('\n\n'.join(map(str, targets)) + '\n\n')


# rules ############
rule all:
    input: targets


# Import all rule files from the rules directory ##########

# Downloads #
include: "rules/download.smk"

if(any(indices_run)):
    decompress(files_to_decompress_fa)

if(any(gtf_derivatives_create)):
    decompress(files_to_decompress_anno)

# Indices #
#include: "rules/bowtie1_index.smk"
#include: "rules/bowtie2_index.smk"
#include: "rules/bwa_index.smk"
#include: "rules/cellranger_index.smk"
#include: "rules/cellranger_vdj_index.smk"
#include: "rules/kallisto_index.smk"
#include: "rules/star_index.smk"
#include: "rules/xengsort_index.smk"


# GTF derivatives #
#include: "rules/bed.smk"
#include: "rules/gene_transcript.smk"
#include: "rules/igv.smk"
#include: "rules/refflat.smk"
