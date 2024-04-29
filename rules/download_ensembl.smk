
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

