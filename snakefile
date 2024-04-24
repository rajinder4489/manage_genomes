# Download the genome (fasta) and gtf files from the Ensembl ftp site (https://ftp.ensembl.org/pub/) 
# For assembly GRch38, works for release >= 47 (all available species)
# For assembly GRch37, works for release >= 83 (only homo_sapiens is availbale)
# Builds genome indices for various mapping tools such as STAR, bwa, bowtie, ...
# Convert gtf to refflat, bed and other formats


import os
import subprocess
from os.path import exists, abspath
from bs4 import BeautifulSoup
import urllib.request
import re
import glob
import random
import shutil


#### Colors for the print
# ANSI color escape codes
RED_START = "\033[91m"
GREEN_START = "\033[92m"
BLUE_START = "\033[94m"
COLOR_END = "\033[0m"  # Reset color to default

configfile: "config.yaml"

# set seed for reproducibility
random.seed(config["seed"])
# base url for fetching data
BASE_URL = "https://ftp.ensembl.org/pub"


#########################################
# Parameters from the config and checks #
#########################################

# resources base path ############
GENOME_DOWNLOAD_PATH = config["genome"]["local_path"]
if not exists(GENOME_DOWNLOAD_PATH):
    raise WorkflowError("Could not find the genome download directory '" + GENOME_DOWNLOAD_PATH + "' specified via 'genome:local_path'!")
GENOME_DOWNLOAD_PATH = "."
GENOME_DOWNLOAD_PATH.rstrip(os.sep)
GENOME_DOWNLOAD_PATH = abspath(GENOME_DOWNLOAD_PATH)


INDICES_BUILD_PATH = config["indices"]["local_path"]
if not exists(INDICES_BUILD_PATH):
    raise WorkflowError("Could not find the indices build directory '" + INDICES_BUILD_PATH + "' specified via 'indices:local_path'!")
INDICES_BUILD_PATH = "."
INDICES_BUILD_PATH.rstrip(os.sep)
INDICES_BUILD_PATH = abspath(INDICES_BUILD_PATH)


# assembly ############
ALL_ASSEMBLIES = ['grch37', 'grch38']
ASSEMBLY = config["genome"]["assembly"] # config.get('assembly', '')

if not ASSEMBLY:
    raise WorkflowError("Please provide the assembly via 'assembly' in the config.yaml!")

if ASSEMBLY not in ALL_ASSEMBLIES:
    terminal_width = shutil.get_terminal_size().columns
    assembly_per_line = terminal_width // max(len(assembly) for assembly in ALL_ASSEMBLIES)
    print(f"\n\nAvailable assemblies:")
    for i in range(0, len(ALL_ASSEMBLIES), assembly_per_line):
        assembly = "  ".join(ALL_ASSEMBLIES[i:i+assembly_per_line])
        colored_assembly = f"{GREEN_START}{assembly}{COLOR_END}"
        print(colored_assembly)
    raise ValueError(f"The provided assembly '{ASSEMBLY}' is not available in Ensembl")

if ASSEMBLY == 'grch38':
    ASSEMBLYPATH = ''
else:
    ASSEMBLYPATH = 'grch37'


# release ############
url = f"{BASE_URL}/{ASSEMBLYPATH}/"
with urllib.request.urlopen(url) as response:
    html_content = response.read()

soup = BeautifulSoup(html_content, "html.parser")
ALL_RELEASES = [link.text.strip("/") for link in soup.find_all("a") if "release" in link.text]

RELEASE = config["genome"]["release"]
if not RELEASE:
    raise WorkflowError("Please provide the release via 'release' in the config.yaml!")

if RELEASE not in ALL_RELEASES:
    terminal_width = shutil.get_terminal_size().columns
    release_per_line = terminal_width // max(len(release) for release in ALL_RELEASES)
    print(f"\n\nAvailable release for assembly '{ASSEMBLY}':")
    for i in range(0, len(ALL_RELEASES), release_per_line):
        release = "  ".join(ALL_RELEASES[i:i+release_per_line])
        colored_release = f"{BLUE_START}{release}{COLOR_END}"
        print(colored_release)
    raise ValueError(f"The provided release '{RELEASE}' is not available in Ensembl.")


# species ############
url = f"{BASE_URL}/{ASSEMBLYPATH}/{RELEASE}/fasta/"
with urllib.request.urlopen(url) as response:
    html_content = response.read()

soup = BeautifulSoup(html_content, "html.parser")
ALL_SPECIES = [link.text.strip("/") for link in soup.find_all("a") if "/" in link.text and link.text.endswith("/")]

SPECIES = config["genome"]["species"]
if not SPECIES:
    raise WorkflowError("Please provide the species via 'species' in the config.yaml!")

if SPECIES not in ALL_SPECIES:
    terminal_width = shutil.get_terminal_size().columns
    species_per_line = terminal_width // max(len(species) for species in ALL_SPECIES)
    print(f"\n\nAvailable species for assembly '{ASSEMBLY}' and release '{RELEASE}':")
    for i in range(0, len(ALL_SPECIES), species_per_line):
        species = "   ".join(ALL_SPECIES[i:i+species_per_line])
        colored_species = f"{RED_START}{species}{COLOR_END}"
        print(colored_species)
    raise ValueError(f"The provided species '{SPECIES}' is not available in Ensembl.")


# seq type ###########
SEQTYPE = config["genome"]["seq_type"]
if not SEQTYPE:
    raise WorkflowError("Please provide the seq type via 'seq_type' in the config.yaml!")


# get fasta files ###########
url = f"{BASE_URL}/{ASSEMBLYPATH}/{RELEASE}/fasta/{SPECIES}/{SEQTYPE}/"
with urllib.request.urlopen(url) as response:
    html_content = response.read()

soup = BeautifulSoup(html_content, "html.parser")
ALL_FILES = list([link['href'] for link in soup.find_all('a') if link.get('href') and not link['href'].endswith('/')])

# patterns fasta file ########
fasta_patterns_include = config["genome"]["fasta"]["file_patterns_include"]
fasta_patterns_exclude = config["genome"]["fasta"]["file_patterns_exclude"]

fasta_to_keep = re.compile('|'.join(fasta_patterns_include))
fasta_to_remove = re.compile('|'.join(fasta_patterns_exclude))

fasta_download = [s for s in ALL_FILES if re.search(fasta_to_keep, s)]
fasta_download = [s for s in fasta_download if not re.search(fasta_to_remove, s)]
only_fasta = [re.sub(r'\.gz$', '', fasta) for fasta in fasta_download if re.search(r'\.fa\.gz$', fasta)]


# get anno files ###########
url = f"{BASE_URL}/{ASSEMBLYPATH}/{RELEASE}/gtf/{SPECIES}/"
with urllib.request.urlopen(url) as response:
    html_content = response.read()

soup = BeautifulSoup(html_content, "html.parser")
ALL_FILES_ANNO = list([link['href'] for link in soup.find_all('a') if link.get('href') and not link['href'].endswith('/')])

# patterns anno file ########
anno_patterns_include = config["genome"]["annotation"]["file_patterns_include"]
anno_patterns_exclude = config["genome"]["annotation"]["file_patterns_exclude"]

anno_to_keep = re.compile('|'.join(anno_patterns_include))
anno_to_remove = re.compile('|'.join(anno_patterns_exclude))

anno_download = [s for s in ALL_FILES_ANNO if re.search(anno_to_keep, s)]
anno_download = [s for s in anno_download if not re.search(anno_to_remove, s)]
anno_extensionless = [os.path.splitext(os.path.splitext(path)[0])[0] for path in anno_download]


print("\n\n" + "Running for " + ASSEMBLY + " " + RELEASE + " " + SPECIES + " " + SEQTYPE)


#################
#### Tagrets ####
#################


# Define paths for genome-related files
if config["genome"]["fasta"]["download"]:
    ALL_FA = expand(
        [os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, SEQTYPE, name) for name in only_fasta]
    )

if config["genome"]["annotation"]["download"]:
    ALL_ANNO = expand(
        [os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", name, ".gtf") for name in anno_extensionless]
    )

# Define paths for index files based on configuration
if config["indices"]["bowtie1"]["run"]:
    ALL_BOWTIE1 = [
        os.path.join(
            INDICES_BUILD_PATH, "indices", "bowtie1", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.{i}.ebwt"
        ) for i in range(1, 5)
    ] + [
        os.path.join(
            INDICES_BUILD_PATH, "indices", "bowtie1", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.rev.{i}.ebwt"
        ) for i in range(1, 3)
    ]

if config["indices"]["bowtie2"]["run"]:
    ALL_BOWTIE2 = [
        os.path.join(
            INDICES_BUILD_PATH, "indices", "bowtie2", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.{i}.bt2"
        ) for i in range(1, 5)
    ] + [
        os.path.join(
            INDICES_BUILD_PATH, "indices", "bowtie2", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.rev.{i}.bt2"
        ) for i in range(1, 3)
    ]

if config["indices"]["bwa"]["run"]:
    ALL_BWA = os.path.join(
        INDICES_BUILD_PATH, "indices", "bwa", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.sa"
    )

if config["indices"]["kallisto"]["run"]:
    ALL_KALLISTO = os.path.join(
        INDICES_BUILD_PATH, "indices", "kallisto", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.idx"
    )

if config["indices"]["star"]["run"]:
    ALL_STAR = os.path.join(
        INDICES_BUILD_PATH, "indices", "star", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "SA"
    )

if config["indices"]["cellranger"]["run"]:
    ALL_CELLRANGER = os.path.join(
        INDICES_BUILD_PATH, "indices", "cellranger", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "**"
    )

if config["indices"]["cellranger_vdj"]["run"]:
    ALL_CELLRANGER_VDJ = os.path.join(
        INDICES_BUILD_PATH, "indices", "cellranger_vdj", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "**"
    )

# Define paths for additional annotation files based on configuration
if config["annotation_files"]["gene_transcript_relation"]:
    GT = expand(
        os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}_gene-transcript.txt"),
        file=anno_extensionless
    )

if config["annotation_files"]["bed"]:
    BED = expand(
        os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}.bed"),
        file=anno_extensionless
    )

if config["annotation_files"]["refflat"]:
    REFFLAT = expand(
        os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}.refflat"),
        file=anno_extensionless
    )

if config["annotation_files"]["igv"]:
    IGV_GTF = expand(
        os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}.igv.gtf"),
        file=anno_extensionless
    )



if config["genome"]["fasta"]["download"]: ALL_FA = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "{file}"), file = only_fasta)
if config["genome"]["annotation"]["download"]: ALL_ANNO = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}.gtf"), file = anno_extensionless)

if config["indices"]["bowtie1"]["run"]: ALL_BOWTIE1 = [os.path.join(INDICES_BUILD_PATH, "indices", "bowtie1", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.{i}.ebwt") for i in range(1, 5)],
        [os.path.join(INDICES_BUILD_PATH, "indices", "bowtie1", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.rev.{i}.ebwt") for i in range(1, 3)]
if config["indices"]["bowtie2"]["run"]: ALL_BOWTIE2 = [os.path.join(INDICES_BUILD_PATH, "indices", "bowtie2", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.{i}.bt2") for i in range(1, 5)] +
        [os.path.join(INDICES_BUILD_PATH, "indices", "bowtie2", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.rev.{i}.bt2") for i in range(1, 3)]
if config["indices"]["bwa"]["run"]: ALL_BWA = os.path.join(INDICES_BUILD_PATH, "indices", "bwa", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.sa")
if config["indices"]["kallisto"]["run"]: ALL_KALLISTO = os.path.join(INDICES_BUILD_PATH, "indices", "kallisto", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.idx")
if config["indices"]["star"]["run"]: ALL_STAR = os.path.join(INDICES_BUILD_PATH, "indices", "star", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "SA")
if config["indices"]["cellranger"]["run"]: ALL_CELLRANGER = os.path.join(INDICES_BUILD_PATH, "indices", "cellranger", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "**")
if config["indices"]["cellranger_vdj"]["run"]: ALL_CELLRANGER_VDJ = os.path.join(INDICES_BUILD_PATH, "indices", "cellranger_vdj", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "**")

if config["annotation_files"]["gene_transcript_relation"]: GT = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}_gene-transcript.txt"), file = anno_extensionless)
if config["annotation_files"]["bed"]: BED = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}.bed"), file = anno_extensionless)
if config["annotation_files"]["refflat"]: REFFLAT = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}.refflat"), file = anno_extensionless)
if config["annotation_files"]["igv"]: IGV_GTF = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}.igv.gtf"), file = anno_extensionless)


#################
##### Rules #####
#################


rule all:
    input:
        ALL_FA,
        ALL_ANNO,
        ALL_BOWTIE1,
        ALL_BOWTIE2,
        ALL_BWA,
        ALL_KALLISTO,
        ALL_STAR,
        ALL_CELLRANGER,
        ALL_CELLRANGER_VDJ,
        GT,
        BED,
        REFFLAT,
        IGV_GTF


        #expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "{file}"), file=only_fasta), 
        #[os.path.join(INDICES_BUILD_PATH, "indices", "bowtie1", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.{i}.ebwt") for i in range(1, 5)] +
        #[os.path.join(INDICES_BUILD_PATH, "indices", "bowtie1", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.rev.{i}.ebwt") for i in range(1, 3)] +
        #[os.path.join(INDICES_BUILD_PATH, "indices", "bowtie2", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.{i}.bt2") for i in range(1, 5)] +
        #[os.path.join(INDICES_BUILD_PATH, "indices", "bowtie2", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.rev.{i}.bt2") for i in range(1, 3)],
        #os.path.join(INDICES_BUILD_PATH, "indices", "bwa", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.sa"),
        #os.path.join(INDICES_BUILD_PATH, "indices", "star", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "SA"),
        #os.path.join(INDICES_BUILD_PATH, "indices", "kallisto", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.idx"),
        #expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}.gtf"), file = anno_extensionless),
        #expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}_gene-transcript.txt"), file = anno_extensionless)

# downloads #######
rule download_genome:
    output:
        expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "{file}"), file = fasta_download)
    run:
        if config["genome"]["fasta"]["download"]:
            print("\n\nI will now try to download the fasta file(s) and if already mentioned, CHECKSUMS...\n\n")
            for file in fasta_download:
                print(file)
                url = f"{BASE_URL}/{ASSEMBLYPATH}/{RELEASE}/fasta/{SPECIES}/{SEQTYPE}/{file}"
                filename = os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, SEQTYPE, file)
                urllib.request.urlretrieve(url, filename)
        else:
            print("Skipping download_genome rule.")


rule download_annotation:
    output:
        expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}"), file = anno_download)

    run:
        if config["genome"]["annotation"]["download"]:
            print("\n\nI will now try to download the annotation file(s) and if already mentioned, CHECKSUMS...\n\n")
            for file in anno_download:
                url = f"{BASE_URL}/{ASSEMBLYPATH}/{RELEASE}/gtf/{SPECIES}/{file}"
                filename = os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", file)
                urllib.request.urlretrieve(url, filename)
        else:
            print("Skipping download_annotation rule.")


# Index building ###########
rule decompress:
    input:
        fasta_files_gz = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "{file}.gz"), file = only_fasta),
        gtf_file_gz = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}.gtf.gz"), file = anno_extensionless)

    output:
        fasta_files = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "{file}"), file = only_fasta),
        gtf_file = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}.gtf"), file = anno_extensionless)

    run:
        shell(
            """
            gzip -df {input}
            """
            )


# bowtie ##########
rule bowtie1_index:
    input:
        fasta = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "{file}"), file = only_fasta)
    
    output:
        [os.path.join(INDICES_BUILD_PATH, "indices", "bowtie1", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.{i}.ebwt") for i in range(1, 5)] +
        [os.path.join(INDICES_BUILD_PATH, "indices", "bowtie1", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.rev.{i}.ebwt") for i in range(1, 3)]

    params:
        fasta_files = lambda wildcards, input: ",".join(f for f in input),
        params = config["indices"]["bowtie1"]["tool_params"]

    run:
        if config["indices"]["bowtie1"]["run"]:
            print("Making the bowtie1 indices")
            shell(
                """
                module load apps/bowtie
                mkdir -p {INDICES_BUILD_PATH}/indices/bowtie1/{SPECIES}/{ASSEMBLY}/{RELEASE}/{SEQTYPE}/
                bowtie-build {params.params} {params.fasta_files} {INDICES_BUILD_PATH}/indices/bowtie1/{SPECIES}/{ASSEMBLY}/{RELEASE}/{SEQTYPE}/{SPECIES}
                module unload apps/bowtie
                """
                )


# bowtie2 ##########
rule bowtie2_index:
    input:
        fasta = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "{file}"), file = only_fasta)
    
    output:
        [os.path.join(INDICES_BUILD_PATH, "indices", "bowtie2", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.{i}.bt2") for i in range(1, 5)] +
        [os.path.join(INDICES_BUILD_PATH, "indices", "bowtie2", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.rev.{i}.bt2") for i in range(1, 3)]

    params:
        fasta_files = lambda wildcards, input: ",".join(f for f in input),
        params = config["indices"]["bowtie2"]["tool_params"]

    run:
        if config["indices"]["bowtie2"]["run"]:
            print("Making the bowtie2 indices")
            shell(
                """
                module load apps/bowtie2
                mkdir -p {INDICES_BUILD_PATH}/indices/bowtie2/{SPECIES}/{ASSEMBLY}/{RELEASE}//{SEQTYPE}/
                bowtie2-build {params.params} {params.fasta_files} {INDICES_BUILD_PATH}/indices/bowtie2/{SPECIES}/{ASSEMBLY}/{RELEASE}/{SEQTYPE}/{SPECIES}
                module unload apps/bowtie2
                """
                )


# bwa ##########
rule bwa_index:
    input:
        fasta = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "{file}"), file = only_fasta)
    
    output:
        os.path.join(INDICES_BUILD_PATH, "indices", "bwa", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.sa")

    params:
        params = config["indices"]["bwa"]["tool_params"]

    run:
        if config["indices"]["bwa"]["run"]:
            print("Making the bwa indices")
            shell(
                """
                module load apps/bwa
                mkdir -p {INDICES_BUILD_PATH}/indices/bwa/{SPECIES}/{ASSEMBLY}/{RELEASE}/
                bwa index {params.params} -p {INDICES_BUILD_PATH}/indices/bwa/{SPECIES}/{ASSEMBLY}/{RELEASE}/{SEQTYPE}/{SPECIES} <(cat {input})
                module unload apps/bwa
                """
                )


# kallisto ##########
rule kallisto_index:
    input:
        fasta = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "{file}"), file = only_fasta)
            
    output:
        os.path.join(INDICES_BUILD_PATH, "indices", "kallisto", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.idx")

    params:
        params = config["indices"]["kallisto"]["tool_params"]

    run:
        if config["indices"]["kallisto"]["run"]:
            print("Making the kallisto indices")
            shell(
                """
                module load apps/kallisto
                mkdir -p {INDICES_BUILD_PATH}/indices/kallisto/{SPECIES}/{ASSEMBLY}/{RELEASE}/{SEQTYPE}/
                kallisto index -i {INDICES_BUILD_PATH}/indices/kallisto/{SPECIES}/{ASSEMBLY}/{RELEASE}/{SEQTYPE}/{SPECIES}.idx {input.fasta}
                module unload apps/kallisto
                """
                )


# STAR ##########
rule star_index:
    input:
        fasta = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "{file}"), file = only_fasta),
        gtf_file = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}.gtf"), file = anno_extensionless)
            
    output:
        os.path.join(INDICES_BUILD_PATH, "indices", "star", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "SA")

    params:
        fasta_files = lambda wildcards, input: " ".join(f for f in input.fasta),
        params = config["indices"]["star"]["tool_params"]

    run:
        if config["indices"]["star"]["run"]:
            print("Making the STAR indices")
            shell(
                """
                module load apps/star
                mkdir -p {INDICES_BUILD_PATH}/indices/star/{SPECIES}/{ASSEMBLY}/{RELEASE}/{SEQTYPE}
                STAR --runMode genomeGenerate {params.params} --genomeDir {INDICES_BUILD_PATH}/indices/star/{SPECIES}/{ASSEMBLY}/{RELEASE}/{SEQTYPE} --sjdbGTFfile {input.gtf_file} --genomeFastaFiles {params.fasta_files}
                module unload apps/star
                """
                )


# derivatives of gtf ###########

# gene-transcript-relation file (kallisto)
rule create_gene_transcript_file:
    input:
        expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}.gtf"), file = anno_extensionless)
    
    output:
        expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}_gene-transcript.txt"), file = anno_extensionless)
        #os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{anno_extensionless}_gene-transcript.txt")

    shell:
        """
        awk -F'\t' '$3=="transcript" {{
            split($9, fields, ";");
            for (i in fields) {{
                split(fields[i], kv, " ");
                if (kv[1] == "gene_id") gi=substr(kv[2], 2, length(kv[2])-2);
                else if (kv[1] == "gene_name") gn=substr(kv[2], 2, length(kv[2])-2);
                else if (kv[1] == "transcript_id") ti=substr(kv[2], 2, length(kv[2])-2);
            }}
            print gi"\t"ti"\t"gn;
        }}' {input} > {output}
        """


# bed file
rule create_bed_file:
    input:
        os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{anno_extensionless}.gtf")
    
    #output:
    #    os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{anno_extensionless}.bed")

    shell:
        """
        
        """


# refflat file
rule create_refflat_file:
    input:
        os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{anno_extensionless}.gtf")
    
    #output:
    #    os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{anno_extensionless}.refflat")

    shell:
        """
        
        """


# some file
rule create_some_file:
    input:
        os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{anno_extensionless}.gtf")
    
    #output:
    #    os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{anno_extensionless}.some")

    shell:
        """
        
        """