import os

# bowtie2 ##########
rule bowtie2_index:
    input:
        fasta_files = os.path.join(fasta_download_path, "{s}", "{a}_{r}_{t}", "genome.fa")

    output:
        indices = [os.path.join(bowtie2_indices_path, "{s}", "{a}_{r}_{t}", "{s}.4.bt2") if config["build_indices"]["bowtie2"]["run"] else "."] +
                    [os.path.join(bowtie2_indices_path, "{s}", "{a}_{r}_{t}", "{s}.rev.2.bt2") if config["build_indices"]["bowtie2"]["run"] else "."]

    wildcard_constraints:
        t = "dna|cdna|cds|ncrna"
    
    params:
        params = config["build_indices"]["bowtie2"]["tool_params"]

    log:
        "logs/bowtie2_index/{s}/{a}_{r}_{t}/log.log"

    run:
        if config["build_indices"]["bowtie2"]["run"]:
            shell(
                """
                (
                    echo "Making the bowtie2 indices"
                    module load apps/bowtie2
                    mkdir -p {bowtie2_indices_path}/{wildcards.s}/{wildcards.a}_{wildcards.r}_{wildcards.t}
                    bowtie2-build {params.params} {input.fasta_files} {bowtie2_indices_path}/{wildcards.s}/{wildcards.a}_{wildcards.r}_{wildcards.t}/{wildcards.s}
                    module unload apps/bowtie2
                ) &> {log}
                """
                )
