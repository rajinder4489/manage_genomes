# kallisto ##########
rule kallisto_index:
    input:
        fasta_files = os.path.join(fasta_download_path, "{s}", "{a}_{r}_{t}", "genome.fa")
            
    output:
        indices = os.path.join(kallisto_indices_path, "{s}", "{a}_{r}_{t}", "{s}.idx")

    wildcard_constraints:
        t = "dna|cdna|cds|ncrna"
    
    params:
        params = config["build_indices"]["kallisto"]["tool_params"]

#    log:
#        "logs/{wildcards.s}_{wildcards.a}_{wildcards.r}_{wildcards.t}/kallisto_index.log"

    run:
        if config["build_indices"]["kallisto"]["run"]:
            shell(
                """
                (
                    echo "Making the kallisto indices"
                    module load apps/kallisto
                    mkdir -p {kallisto_indices_path}/{wildcards.s}/{wildcards.a}_{wildcards.r}_{wildcards.t}
                    kallisto index -i {kallisto_indices_path}/{wildcards.s}/{wildcards.a}_{wildcards.r}_{wildcards.t}/{wildcards.s}.idx {input.fasta_files}
                    module unload apps/kallisto
                )
                """
                )
# &> {log}