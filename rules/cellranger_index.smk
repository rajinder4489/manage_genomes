# cellranger ##########
rule cellranger_index:
    input:
        fasta_files = os.path.join(fasta_download_path, "{s}", "{a}_{r}_{t}", "genome.fa"),
        gtf_file = os.path.join(fasta_download_path, "{s}", "{a}_{r}_annotation", "gtf.gtf")
        
    output:
        indices = os.path.join(cellranger_indices_path, "{s}", "{a}_{r}_{t}", "star", "SA")

    wildcard_constraints:
        t = "dna|cdna|cds|ncrna"
    
    params:
        params = config["build_indices"]["cellranger"]["tool_params"]

#    log:
#        "logs/cellranger_index/{s}/{a}_{r}_{t}/log.log"
    
    run:
        if config["build_indices"]["cellranger"]["run"]:
            shell(
                """
                    echo "Making the Cellranger indices"
                    module load apps/cellranger
                    cd {cellranger_indices_path}/{wildcards.s}/{wildcards.a}_{wildcards.r}_{wildcards.t}
                    cellranger mkref --genome {wildcards.s} --fasta {input.fasta_files} --genes {input.gtf_file} --output-dir {cellranger_indices_path}/{wildcards.s}/{wildcards.a}_{wildcards.r}_{wildcards.t}
                    module unload apps/cellranger
                """
                )