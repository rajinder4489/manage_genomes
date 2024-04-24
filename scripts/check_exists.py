import os
import shutil
import urllib.request
import re
from bs4 import BeautifulSoup

def assembly_check(config):
    ALL_ASSEMBLIES = ['grch37', 'grch38']
    ASSEMBLY = config["genome"]["assembly"]
    
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

def release_check(config):
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

def species_check(config):
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

def seq_type_check(config):
    SEQTYPE = config["genome"]["seq_type"]
    if not SEQTYPE:
        raise WorkflowError("Please provide the seq type via 'seq_type' in the config.yaml!")

def fetch_files(config):
    url = f"{BASE_URL}/{ASSEMBLYPATH}/{RELEASE}/fasta/{SPECIES}/{SEQTYPE}/"
    with urllib.request.urlopen(url) as response:
        html_content = response.read()

    soup = BeautifulSoup(html_content, "html.parser")
    ALL_FILES = list([link['href'] for link in soup.find_all('a') if link.get('href') and not link['href'].endswith('/')])

    fasta_patterns_include = config["genome"]["fasta"]["file_patterns_include"]
    fasta_patterns_exclude = config["genome"]["fasta"]["file_patterns_exclude"]

    fasta_to_keep = re.compile('|'.join(fasta_patterns_include))
    fasta_to_remove = re.compile('|'.join(fasta_patterns_exclude))

    fasta_download = [s for s in ALL_FILES if re.search(fasta_to_keep, s)]
    fasta_download = [s for s in fasta_download if not re.search(fasta_to_remove, s)]
    only_fasta = [re.sub(r'\.gz$', '', fasta) for fasta in fasta_download if re.search(r'\.fa\.gz$', fasta)]
    
    return only_fasta

def fetch_annotation_files(config):
    url = f"{BASE_URL}/{ASSEMBLYPATH}/{RELEASE}/gtf/{SPECIES}/"
    with urllib.request.urlopen(url) as response:
        html_content = response.read()

    soup = BeautifulSoup(html_content, "html.parser")
    ALL_FILES_ANNO = list([link['href'] for link in soup.find_all('a') if link.get('href') and not link['href'].endswith('/')])

    anno_patterns_include = config["genome"]["annotation"]["file_patterns_include"]
    anno_patterns_exclude = config["genome"]["annotation"]["file_patterns_exclude"]

    anno_to_keep = re.compile('|'.join(anno_patterns_include))
    anno_to_remove = re.compile('|'.join(anno_patterns_exclude))

    anno_download = [s for s in ALL_FILES_ANNO if re.search(anno_to_keep, s)]
    anno_download = [s for s in anno_download if not re.search(anno_to_remove, s)]
    anno_extensionless = [os.path.splitext(os.path.splitext(path)[0])[0] for path in anno_download]
    
    return anno_extensionless

def print_run_info():
    print("\n\n" + "Running for " + ASSEMBLY + " " + RELEASE + " " + SPECIES + " " + SEQTYPE)
