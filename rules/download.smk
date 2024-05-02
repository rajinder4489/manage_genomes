rule download_annotation_ensembl:
    output:
        dynamic(os.path.join(annotation_download_path, "{species}", "{assembly}", "{release}", "annotation", "{a_file}")) if config["ensembl"]["annotation"]["download"] else "."

    wildcard_constraints:
        a_file = ".*.gtf.gz|CHECKSUMS"
    
    run:
        if config["ensembl"]["annotation"]["download"]:
            for key, values in downloadable_files.items():
                assembly = values["assembly"]
                release = values["release"]
                species = values["species"]
                anno_files = values["annotation"]

                if assembly == 'grch38':
                    assembly_path = ''
                elif assembly == 'grch37':
                    assembly_path = 'grch37'
                
                for a_file in anno_files:
                    url = f"{base_url}/{assembly_path}/{release}/gtf/{species}/{a_file}"
                    print(f"{a_file} : {url}")
                    filename = os.path.join(annotation_download_path, species, assembly, release, "annotation", a_file)
                    urllib.request.urlretrieve(url, filename)


rule download_genome_ensembl:
    output:
        dynamic(os.path.join(fasta_download_path, "{species}", "{assembly}", "{release}", "{seqtype}", "{g_file}")) if config["ensembl"]["fasta"]["download"] else "."

    wildcard_constraints:
        g_file = ".*.fa.gz|CHECKSUMS",
        seq_type = "dna|cdna|cds|ncrna|pep|dna_index"
    
    run:
        if config["ensembl"]["fasta"]["download"]:
            for key, values in downloadable_files.items():
                assembly = values["assembly"]
                release = values["release"]
                species = values["species"]
                seqtype = values["seqtype"]
                genome_files = values["fasta"]


                if assembly == 'grch38':
                    assembly_path = ''
                elif assembly == 'grch37':
                    assembly_path = 'grch37'

                for g_file in genome_files:
                    url = f"{base_url}/{assembly_path}/{release}/fasta/{species}/{seqtype}/{g_file}"
                    filename = os.path.join(fasta_download_path, species, assembly, release, seqtype, g_file)
                    urllib.request.urlretrieve(url, filename)


rule download_repeats_ensembl:
    output:
        dynamic(os.path.join(ensembl_repeats_download_path, "{species}", "{assembly}", "{release}", "{repeats}", "{g_file}")) if config["ensembl"]["fasta"]["download"] else "."

    wildcard_constraints:
        g_file = ".*.fa.gz|CHECKSUMS",
        seq_type = "dna|cdna|cds|ncrna|pep|dna_index"
    
    run:
        if config["ensembl"]["fasta"]["download"]:
            for key, values in downloadable_files.items():
                assembly = values["assembly"]
                release = values["release"]
                species = values["species"]


                if assembly == 'grch38':
                    assembly_path = ''
                elif assembly == 'grch37':
                    assembly_path = 'grch37'

                for g_file in genome_files:
                    url = f"{base_url}/{assembly_path}/{release}/fasta/{species}/{seqtype}/{g_file}"
                    filename = os.path.join(fasta_download_path, species, assembly, release, seqtype, g_file)
                    urllib.request.urlretrieve(url, filename)