# bwa ##########
rule bwa_index:
    input:
        fasta_files = glob.glob(os.path.join(fasta_download_path, f"{species}", f"{assembly}_{release}_{seqtype}", "*fa"))
    
    output:
        indices = (os.path.join(bwa_indices_path, f"{species}", f"{assembly}_{release}_{seqtype}", f"{species}.{i}") for i in ["amb", "ann", "bwt", "pac", "sa"]) if config["build_indices"]["bwa"]["run"] else "."

    params:
        params = config["build_indices"]["bwa"]["tool_params"]

    run:
        if config["build_indices"]["bwa"]["run"]:
            print("Making the bwa indices")
            shell(
                """
                module load apps/bwa
                mkdir -p {bwa_indices_path}/{species}/{assembly}_{release}_{seqtype}
                bwa index {params.params} -p {bwa_indices_path}/{species}/{assembly}_{release}_{seqtype}/{species} <(cat {input})
                module unload apps/bwa
                """
                )
