import os
import urllib.request
import re
import warnings
from bs4 import BeautifulSoup


#### Colors for printing the available names
# ANSI color escape codes
DARK_RED = "\033[31m"   # Fail message
DARK_GREEN = "\033[32m" # Success message
COLOR_END = "\033[0m"   # Reset color to default

colored_failed = f"{DARK_RED}Failed: No files left after filtering for the regex{COLOR_END}"
colored_success = f"{DARK_RED}Succeded: Got some files to fetch{COLOR_END}"

class MissingFiles(UserWarning):
    """Class representing missing file(s)."""


def downloadable_ensembl(config, combn_dict):
    """Class to check for the available files"""

    base_url = config["ensembl"]["base_url"]

    for key in combn_dict:

        a = combn_dict[key]['assembly']
        r = combn_dict[key]['release']
        s = combn_dict[key]['species']
        t = combn_dict[key]['seqtype']

        if a == 'grch38':
            assembly_path = ''
        elif a == 'grch37':
            assembly_path = 'grch37'
    
    
        # get fasta files ###########
        url = f"{base_url}/{assembly_path}/{r}/fasta/{s}/{t}/"
        with urllib.request.urlopen(url) as response:
            html_content = response.read()

        soup = BeautifulSoup(html_content, "html.parser")
        all_files = list([link['href'] for link in soup.find_all('a') if link.get('href') and not link['href'].endswith('/')])

        # patterns fasta file ########
        fasta_patterns_include = config["ensembl"]["fasta"]["file_patterns_include"]
        fasta_patterns_exclude = config["ensembl"]["fasta"]["file_patterns_exclude"]

        fasta_to_keep = re.compile('|'.join(fasta_patterns_include))
        fasta_to_remove = re.compile('|'.join(fasta_patterns_exclude))

        fasta_download = [s for s in all_files if re.search(fasta_to_keep, s)]
        fasta_download = [s for s in fasta_download if not re.search(fasta_to_remove, s)]
        # only_fasta = [re.sub(r'\.gz$', '', fasta) for fasta in fasta_download if re.search(r'\.fa\.gz$', fasta)]

        if not fasta_download:
            warnings.warn(colored_failed, MissingFiles)
            combn_dict[key]["fasta"] = ""
        else:
            print(colored_success)
            combn_dict[key]["fasta"] = fasta_download


        # get anno files ###########
        url = f"{base_url}/{assembly_path}/{r}/gtf/{s}/"
        with urllib.request.urlopen(url) as response:
            html_content = response.read()

        soup = BeautifulSoup(html_content, "html.parser")
        all_files_anno = list([link['href'] for link in soup.find_all('a') if link.get('href') and not link['href'].endswith('/')])

        # patterns anno file ########
        anno_patterns_include = config["ensembl"]["annotation"]["file_patterns_include"]
        anno_patterns_exclude = config["ensembl"]["annotation"]["file_patterns_exclude"]

        anno_to_keep = re.compile('|'.join(anno_patterns_include))
        anno_to_remove = re.compile('|'.join(anno_patterns_exclude))

        anno_download = [s for s in all_files_anno if re.search(anno_to_keep, s)]
        anno_download = [s for s in anno_download if not re.search(anno_to_remove, s)]
        # anno_extensionless = [os.path.splitext(os.path.splitext(path)[0])[0] for path in anno_download]

        if not anno_download:
            warnings.warn(colored_failed, MissingFiles)
            combn_dict[key]["annotation"] = ""
        else:
            print(colored_success)
            combn_dict[key]["annotation"] = anno_download
    
    return(combn_dict)