import urllib.request

rule download_annotation_ensembl:
    output:
        #os.path.join(annotation_download_path, "{s}", "{a}_{r}_annotation", "gtf.gtf") if config["ensembl"]["annotation"]["download"] else "."
        os.path.join(annotation_download_path, "{s}", "{a}_{r}_annotation", "gtf.gtf")
    
    #log:
    #    "logs/{s}_{a}_{r}_downloads_ensembl/annotation.log"
    
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
                    filename = os.path.join(annotation_download_path, species, f"{assembly}_{release}_annotation", a_file)
                    urllib.request.urlretrieve(url, filename)
                
                shell(
                    """
                    (
                        ls -lha {fasta_download_path}/{species}/{assembly}_{release}_annotation/*.gtf.gz > {fasta_download_path}/{species}/{assembly}_{release}_annotation/readme.txt
                        zcat {fasta_download_path}/{species}/{assembly}_{release}_annotation/*.gtf.gz > {fasta_download_path}/{species}/{assembly}_{release}_annotation/gtf.gtf
                        rm -rf {fasta_download_path}/{species}/{assembly}_{release}_annotation/*.gtf.gz
                    )
                    """
                    )


rule download_genome_ensembl:
    output:
        os.path.join(fasta_download_path, "{s}", "{a}_{r}_{t}", "genome.fa")

    wildcard_constraints:
        t = "dna|cdna|cds|ncrna|pep|dna_index"
    
    #log:
    #    "logs/{s}_{a}_{r}_download_ensembl/{t}.log"

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
                    filename = os.path.join(fasta_download_path, species, assembly + "_" + release + "_" + seqtype, g_file)
                    urllib.request.urlretrieve(url, filename)
                
                shell(
                    """
                    (
                        ls -lha {fasta_download_path}/{species}/{assembly}_{release}_{seqtype}/*.fa.gz > {fasta_download_path}/{species}/{assembly}_{release}_{seqtype}/readme.txt
                        zcat {fasta_download_path}/{species}/{assembly}_{release}_{seqtype}/*.fa.gz > {fasta_download_path}/{species}/{assembly}_{release}_{seqtype}/genome.fa
                        rm -rf {fasta_download_path}/{species}/{assembly}_{release}_{seqtype}/*.fa.gz
                    )
                    """
                    )


rule download_repeats_ensembl:
    output:
        os.path.join(ensembl_repeats_download_path, "{s}", "{a}_{r}_repeats", "{r_file}")

    wildcard_constraints:
        r_file = ".*.txt.gz",

    #log:
    #    "logs/{s}_{a}_{r}_downloads_ensembl/repeats.log"
    
    run:
        if config["ensembl"]["fasta"]["download"]:
            for key, values in downloadable_files.items():
                assembly = values["assembly"]
                release = values["release"]
                species = values["species"]
                repeats_files = values["repeats"]


                if assembly == 'grch38':
                    assembly_path = ''
                elif assembly == 'grch37':
                    assembly_path = 'grch37'

                for r_file in repeats_files:
                    url = f"{base_url}/{assembly_path}/{release}/mysql/{species}_core_{release}_*/{seqtype}/{g_file}"
                    filename = os.path.join(fasta_download_path, species, assembly + "_" + release + "_" + repeats, r_file)
                    urllib.request.urlretrieve(url, filename)

#                if(config[source]["ensembl_repeats"]["unzip"]):
#                        subprocess.run(["gzip", "-df", filename])