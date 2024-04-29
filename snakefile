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

# print("\n\n" + "Checking " + source + " " + base_url)

# Import all rule files from the rules directory
#include: "rules/bed.smk"


# Import all scripts
# Define a function to filter script files
def get_script_files(directory):
    script_files = []
    for filename in os.listdir(directory):
        if filename.endswith((".py", ".R", ".sh", ".pl")):
            script_files.append(os.path.join(directory, filename))
    return script_files

# Fetch all script files from the scripts directory
script_files = get_script_files("scripts/")

# Include all script files in the workflow
#for script in script_files:
#    include: script


###########################################
## Parameters from the config and checks ##
###########################################

# resources base path ############

if(config[source]["fasta"]["download"] or 
    config[source]["annotation"]["download"] or 
    config[source]["ensembl_repeats"]["download"]):

    genome_download_path = config["local"]["path_genome"]
    if not exists(genome_download_path):
        raise WorkflowError("Could not find the genome download directory '" + genome_download_path + "' specified via 'genome:local_path'!")
    #genome_download_path = "."
    genome_download_path.rstrip(os.sep)
    genome_download_path = abspath(genome_download_path)


if(config["build_indices"]["bowtie1"]["run"] or 
    config["build_indices"]["bowtie2"]["run"] or 
    config["build_indices"]["bwa"]["run"] or 
    config["build_indices"]["star"]["run"] or 
    config["build_indices"]["cellranger"]["run"] or 
    config["build_indices"]["cellranger_vdj"]["run"] or
    config["build_indices"]["kallisto"]["run"] or 
    config["build_indices"]["xengsort"]["run"]):
    
    indices_build_path = config["local"]["path_indices"]
    if not exists(indices_build_path):
        raise WorkflowError("Could not find the indices build directory '" + indices_build_path + "' specified via 'indices:local_path'!")
    # indices_build_path = "."
    indices_build_path.rstrip(os.sep)
    indices_build_path = abspath(indices_build_path)


include: "scripts/check_exists.py"
if source == "ensembl":
    working_combinations = combination_checks_ensembl(config)


print(working_combinations)