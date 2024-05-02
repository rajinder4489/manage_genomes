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

# resources base path ############

if(config[source]["fasta"]["download"] or 
    config[source]["ensembl_repeats"]["download"]):

    if config[source]["fasta"]["local_path"]:
        fasta_download_path = config[source]["fasta"]["local_path"]
    else:
        fasta_download_path = config["local"]["path_genome"]
    
    
    if not fasta_download_path:
        warnings.warn("Fasta download directory not specified in config via '{source}:fasta:local_path' or 'local:path_genome'!. Setting to current directory")
        fasta_download_path = "."
    if not exists(fasta_download_path):
        raise WorkflowError("Could not find the genome download directory '" + fasta_download_path + "' specified via  or '{source}:fasta:local_path' or 'local:path_genome'!")
    
    fasta_download_path.rstrip(os.sep)
    fasta_download_path = abspath(fasta_download_path)


if(config[source]["annotation"]["download"]):

    if config[source]["annotation"]["local_path"]:
        annotation_download_path = config[source]["annotation"]["local_path"]
    else:
        annotation_download_path = config["local"]["path_genome"]
    
    
    if not annotation_download_path:
        warnings.warn("Annotation download directory not specified in config via '{source}:annotation:local_path' or 'local:path_genome'!. Setting to current directory")
        annotation_download_path = "."
    if not exists(annotation_download_path):
        raise WorkflowError("Could not find the annotation download directory '" + annotation_download_path + "' specified via '{source}:annotation:local_path' or 'local:path_genome!")
    
    annotation_download_path.rstrip(os.sep)
    annotation_download_path = abspath(annotation_download_path)


if(config["build_indices"]["bowtie1"]["run"] or 
    config["build_indices"]["bowtie2"]["run"] or 
    config["build_indices"]["bwa"]["run"] or 
    config["build_indices"]["star"]["run"] or 
    config["build_indices"]["cellranger"]["run"] or 
    config["build_indices"]["cellranger_vdj"]["run"] or
    config["build_indices"]["kallisto"]["run"] or 
    config["build_indices"]["xengsort"]["run"]):

    indices_build_path = config["local"]["path_indices"]

    if not indices_build_path:
        warnings.warn("Indices build directory not specified in config via 'local:path_indices'!; setting to current directory for now. If local_path defined for the specific tools, it will be reset to it later.")
        indices_build_path = "."
    if not exists(indices_build_path):
        raise WorkflowError("Could not find the indices build directory '" + indices_build_path + "' specified via 'local:path_indices'!")

    indices_build_path.rstrip(os.sep)
    indices_build_path = abspath(indices_build_path)


# Check what combinations exist
if source == "ensembl":
    working_dict = combination_checks_ensembl(config)
elif source == "refseq":
    working_dict = combination_checks_refseq(config)


# Check what files are available for the valid combinations and patterns
if source == "ensembl":
    downloadable_files = downloadable_ensembl(config, working_dict)
elif source == "refseq":
    downloadable_files = downloadable_refseq(config, working_dict)


# rules ############
targets = []

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
                targets.append(os.path.join(fasta_download_path, species, assembly, release, seqtype, file))

        if(config["ensembl"]["annotation"]["download"]):
            for file in anno_files:
                targets.append(os.path.join(annotation_download_path, species, assembly, release, "annotation", file))


#print('\n\n'.join(map(str, targets)) + '\n\n')


# rules ############
rule all:
    input: targets


# Import all rule files from the rules directory
include: "rules/download.smk"