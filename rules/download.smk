rule download_ensembl
    output:
        dynamic(os.path.join(genome_download_path, "{species}", "{assembly}", "{release}", "{seqtype}", "{file}"))
    run:
        if config["ensembl"]["fasta"]["download"]:
            for key, values in my_dict.items():
                assembly = values["assembly"]
                release = values["release"]
                species = values["species"]
                seqtype = values["seqtype"]
                genome_files = values["fasta"]
                anno_files = values["annotation"]

                for file in genome_files:
                    print(file)
                    url = f"{base_url}/{assembly}/{release}/fasta/{species}/{seqtype}/{file}"
                    filename = os.path.join(genome_download_path, species, assembly, release, seqtype, file)
                    urllib.request.urlretrieve(url, filename)
                
                for file in anno_files:
                    print(file)
                    url = f"{base_url}/{assembly}/{release}/fasta/{species}/{seqtype}/{file}"
                    filename = os.path.join(genome_download_path, species, assembly, release, seqtype, file)
                    urllib.request.urlretrieve(url, filename)
        else:
            print("Skipping download_genome rule.")
