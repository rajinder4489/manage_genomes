# xengsort ##########
rule xengsort_index:
    input:
        fasta_files = os.path.join(fasta_download_path, "{s}", "{a}_{r}_{t}", "genome.fa")
            
    output:
        indices = os.path.join(xengsort_indices_path, "{s}", "{a}_{r}_{t}", "{s}.idx") if config["build_indices"]["xengsort"]["run"] else "."

    wildcard_constraints:
        t = "dna|cdna|cds|ncrna"
    
    params:
        params = config["build_indices"]["xengsort"]["tool_params"]

    log:
        "logs/xengsort_index/{s}/{a}_{r}_{t}/log.log"
    
    run:
        if config["build_indices"]["xengsort"]["run"]:
            shell(
                """
                (
                    echo "Making the Xengsort indices"
                    module load apps/xengsort
                    mkdir -p {xengsort_indices_path}/{all_species}/{all_assembly}_{all_release}_{seqtype}
                    ...
                    module unload apps/xengsort
                ) &> {log}
                """
                )
