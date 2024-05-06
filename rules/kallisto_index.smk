# kallisto ##########
rule kallisto_index:
    input:
        fasta_files = glob.glob(os.path.join(fasta_download_path, f"{species}", f"{assembly}_{release}_{seqtype}", "*fa")),
            
    output:
        indices = (os.path.join(kallisto_indices_path, f"{species}", f"{assembly}_{release}_{seqtype}", f"{species}.idx")) if config["build_indices"]["kallisto"]["run"] else "."

    wildcard_constraints:
        seqtype = "dna|cdna|cds|ncrna"
    
    params:
        params = config["build_indices"]["kallisto"]["tool_params"]

    run:
        if config["build_indices"]["kallisto"]["run"]:
            print("Making the kallisto indices")
            shell(
                """
                module load apps/kallisto
                mkdir -p {kallisto_indices_path}/{species}/{assembly}_{release}_{seqtype}
                kallisto index -i {kallisto_indices_path}/{species}/{assembly}_{release}_{seqtype}/{species}.idx {input.fasta_files}
                module unload apps/kallisto
                """
                )
