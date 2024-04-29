import os
import shutil
import urllib.request
import re
import warnings
from itertools import product
from bs4 import BeautifulSoup


#### Colors for printing the available names
# ANSI color escape codes
GREEN = "\033[95m"      # assembly
BLUE = "\033[94m"       # release
RED = "\033[91m"        # species
YELLOW = "\033[93m"     # seq type
DARK_RED = "\033[31m"   # Fail message
DARK_GREEN = "\033[32m" # Success message
COLOR_END = "\033[0m"   # Reset color to default

colored_failed = f"{DARK_RED}Failed{COLOR_END}"
colored_success = f"{DARK_RED}Success{COLOR_END}"

class MissingValueCombination(UserWarning):
    """Class representing a missing combination."""

def combination_checks_ensembl(config):
    """Class to check for the available assemblies/releases/species/seq types"""

    all_assemblies = ['grch37', 'grch38']
    
    base_url = config["ensembl"]["base_url"]
    
    assembly = config["ensembl"]["assembly"]
    release = config["ensembl"]["release"]
    species = config["ensembl"]["species"]
    seqtype = config["ensembl"]["seq_type"]


    if not assembly:
        raise ValueError("Please provide the assembly via 'assembly' in the config.yaml!")

    if not release:
        raise ValueError("Please provide the release via 'release' in the config.yaml!")

    if not species:
        raise ValueError("Please provide the species via 'species' in the config.yaml!")

    if not seqtype:
        raise ValueError("Please provide the seq type via 'seq_type' in the config.yaml!")

    results = {}
    for idx, (a, r, s, t) in enumerate(product(assembly, release, species, seqtype)):

        print("\n" + "Running for " + a + " " + r + " " + s + " " + t)

        # assembly ############
        if a not in all_assemblies:
            print(colored_failed)
            terminal_width = shutil.get_terminal_size().columns
            assembly_per_line = terminal_width // max(len(assembly) for assembly in all_assemblies)
            print("\nAvailable assemblies:")
            for i in range(0, len(all_assemblies), assembly_per_line):
                assembly = "  ".join(all_assemblies[i:i+assembly_per_line])
                colored_assembly = f"{GREEN}{assembly}{COLOR_END}"
                print(colored_assembly)
            #raise ValueError(f"The provided assembly '{a}' is not available in Ensembl")
            warnings.warn(f"The provided assembly '{a}' is not available in Ensembl")
            continue

        if a == 'grch38':
            assembly_path = ''
        elif a == 'grch37':
            assembly_path = 'grch37'


        # release ############
        url = f"{base_url}/{assembly_path}/"
        with urllib.request.urlopen(url) as response:
            html_content = response.read()

        soup = BeautifulSoup(html_content, "html.parser")
        all_releases = [link.text.strip("/") for link in soup.find_all("a") if "release" in link.text]

        if r not in all_releases:
            print(colored_failed)
            terminal_width = shutil.get_terminal_size().columns
            release_per_line = terminal_width // max(len(release) for release in all_releases)
            print(f"\nAvailable release for assembly '{a}':")
            for i in range(0, len(all_releases), release_per_line):
                release = "  ".join(all_releases[i:i+release_per_line])
                colored_release = f"{BLUE}{release}{COLOR_END}"
                print(colored_release)
            warnings.warn(f"The provided release '{r}' is not available for assembly '{a}' in Ensembl.")
            continue


        # species ############
        url = f"{base_url}/{assembly_path}/{r}/fasta/"
        with urllib.request.urlopen(url) as response:
            html_content = response.read()

        soup = BeautifulSoup(html_content, "html.parser")
        all_species = [link.text.strip("/") for link in soup.find_all("a") if "/" in link.text and link.text.endswith("/")]

        if s not in all_species:
            print(colored_failed)
            terminal_width = shutil.get_terminal_size().columns
            species_per_line = terminal_width // max(len(species) for species in all_species)
            print(f"\nAvailable species for assembly '{a}' and release '{r}':")
            for i in range(0, len(all_species), species_per_line):
                species = "   ".join(all_species[i:i+species_per_line])
                colored_species = f"{RED}{species}{COLOR_END}"
                print(colored_species)
            warnings.warn(f"The provided species '{s}' is not available for assembly '{a}' and release '{r}' in Ensembl.")
            continue


        # seq type ###########
        url = f"{base_url}/{assembly_path}/{r}/fasta/{s}"
        with urllib.request.urlopen(url) as response:
            html_content = response.read()

        soup = BeautifulSoup(html_content, "html.parser")
        all_seqtype = [link.text.strip("/") for link in soup.find_all("a") if "/" in link.text and link.text.endswith("/")]

        if t not in all_seqtype:
            print(colored_failed)
            terminal_width = shutil.get_terminal_size().columns
            seqtype_per_line = terminal_width // max(len(seqtype) for seqtype in all_seqtype)
            print(f"\nAvailable types for assembly '{a}', release '{r}', and species '{s}':")
            for i in range(0, len(all_seqtype), seqtype_per_line):
                seqtype = "   ".join(all_seqtype[i:i+seqtype_per_line])
                colored_seqtype = f"{YELLOW}{seqtype}{COLOR_END}"
                print(colored_seqtype)
            warnings.warn(f"The provided type '{t}' is not available for assembly '{a}', release '{r}', and species '{s}' in Ensembl.")
            continue

        print(colored_success)
        
        results[idx] = {
        'assembly': a,
        'release': r,
        'species': s,
        'seqtype': t 
        }

    return(results)
