# xengsort ##########
rule xengsort_index:
    input:
        fasta_files_host = os.path.join(fasta_download_path, "{host_s}", "{host_a}_{host_r}_{host_t}", "genome.fa"),
        fasta_files_graft = os.path.join(fasta_download_path, "{graft_s}", "{graft_a}_{graft_r}_{graft_t}", "genome.fa")

    output:
        indices = os.path.join(xengsort_indices_path, "{host_s}_{host_a}_{graft_s}_{graft_a}", "xengsort")
    
    wildcard_constraints:
        host_t = "dna|cdna|cds|ncrna",
        graft_t = "dna|cdna|cds|ncrna"
    
    params:
        params = config["build_indices"]["xengsort"]["tool_params"]

    #log:
    #    "logs/xengsort_index/{s}/{a}_{r}_{t}/log.log"
    
    run:
        if config["build_indices"]["xengsort"]["run"]:
            shell(
                """
                    echo "Making the Xengsort indices"
                    module load apps/xengsort
                    current_dir=$(pwd)
                    mkdir -p {xengsort_indices_path}/{wildcards.host_s}_{wildcards.host_a}_{wildcards.graft_s}_{wildcards.graft_a}
                    cd {xengsort_indices_path}/{wildcards.s}/{wildcards.a}_{wildcards.r}_{wildcards.s}
                    xengsort index --index xengsort_index -H input.fasta_files.host -G input.fasta_files.graft {params.params}
                    cd "$current_dir"
                    module unload apps/xengsort
                """
                )
