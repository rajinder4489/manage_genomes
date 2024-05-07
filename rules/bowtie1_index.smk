import os
import glob

# bowtie ##########
rule bowtie1_index:
    input:
        fasta_files = os.path.join(fasta_download_path, "{s}", "{a}_{r}_{t}", "genome.fa")

    output:
        indices = [os.path.join(bowtie1_indices_path, "{s}", "{a}_{r}_{t}", "{s}.4.ebwt") if config["build_indices"]["bowtie1"]["run"] else "."] +
            [os.path.join(bowtie1_indices_path, "{s}", "{a}_{r}_{t}", "{s}.rev.2.ebwt") if config["build_indices"]["bowtie1"]["run"] else "."]

    wildcard_constraints:
        seqtype = "dna|cdna|cds|ncrna"
    
    params:
        params = config["build_indices"]["bowtie1"]["tool_params"]

    log:
        "logs/bowtie1_indices/{s}/{a}_{r}_{t}/log.log"

    run:
        if config["build_indices"]["bowtie1"]["run"]:
            shell(
                """
                (
                    echo "Making the bowtie1 indices"
                    module load apps/bowtie
                    mkdir -p {bowtie1_indices_path}/{wildcards.s}/{wildcards.a}_{wildcards.r}_{wildcards.t}
                    bowtie-build {params.params} {input.fasta_files} {bowtie1_indices_path}/{wildcards.s}/{wildcards.a}_{wildcards.r}_{wildcards.t}/{wildcards.s}
                    module unload apps/bowtie
                ) &> {log}
                """
                )