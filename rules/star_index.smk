# STAR ##########
rule star_index:
    input:
        fasta_files = glob.glob(os.path.join(fasta_download_path, f"{species}", f"{assembly}_{release}_{seqtype}", "*fa")),
        gtf_file = glob.glob(os.path.join(annotation_download_path, f"{species}", f"{assembly}_{release}_{seqtype}", "*gtf"))
            
    output:
        indices = (os.path.join(star_indices_path, f"{species}", f"{assembly}_{release}_{seqtype}", "SA")) if config["build_indices"]["star"]["run"] else "."

    params:
        fasta_files = lambda wildcards, input: " ".join(f for f in input.fasta_files),
        params = config["build_indices"]["star"]["tool_params"]

    run:
        if config["build_indices"]["star"]["run"]:
            print("Making the STAR indices")
            shell(
                """
                module load apps/star
                mkdir -p {star_indices_path}/{species}/{assembly}_{release}_{seqtype}
                STAR --runMode genomeGenerate {params.params} --genomeDir {star_indices_path}/{species}/{assembly}_{release}_{seqtype} --sjdbGTFfile {input.gtf_file} --genomeFastaFiles {params.fasta_files}
                module unload apps/star
                """
                )
