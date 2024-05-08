# bwa ##########
rule bwa_index:
    input:
        fasta_files = os.path.join(fasta_download_path, "{s}", "{a}_{r}_{t}", "genome.fa")
    
    output:
        indices = os.path.join(bwa_indices_path, "{s}", "{a}_{r}_{t}", "{s}.sa")

    wildcard_constraints:
        t = "dna|cdna|cds|ncrna"
    
    params:
        params = config["build_indices"]["bwa"]["tool_params"]

#    log:
#        "logs/{s}_{a}_{r}_{t}/bwa_index.log"
    
    run:
        if config["build_indices"]["bwa"]["run"]:
            shell(
                """
                (
                    echo "Making the bwa indices"
                    module load apps/bwa
                    mkdir -p {bwa_indices_path}/{wildcards.s}/{wildcards.a}_{wildcards.r}_{wildcards.t}
                    bwa index {params.params} -p {bwa_indices_path}/{wildcards.s}/{wildcards.a}_{wildcards.r}_{wildcards.t}/{wildcards.s} <(cat {input.fasta_files})
                    module unload apps/bwa
                )
                """
                )
# &> {log}