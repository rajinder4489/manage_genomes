import os
import glob

# bowtie ##########
rule bowtie1_index:
    input:
        fasta_files = glob.glob(os.path.join(fasta_download_path, f"{species}", f"{assembly}_{release}_{seqtype}", "*fa"))

    output:
        indices = [(os.path.join(bowtie1_indices_path, f"{species}", f"{assembly}_{release}_{seqtype}", f"{species}.{i}.ebwt") for i in range(1, 5)) if config["build_indices"]["bowtie1"]["run"] else "."] +
                    [(os.path.join(bowtie1_indices_path, f"{species}", f"{assembly}_{release}_{seqtype}", f"{species}.rev.{i}.ebwt") for i in range(1, 3)) if config["build_indices"]["bowtie1"]["run"] else "."]

    wildcard_constraints:
        seqtype = "dna|cdna|cds|ncrna"
    
    params:
        fasta_files = lambda wildcards, input: ",".join(f for f in input),
        params = config["build_indices"]["bowtie1"]["tool_params"]

    run:
        if config["build_indices"]["bowtie1"]["run"]:
            print("Making the bowtie1 indices")
            print(os.path.join(fasta_download_path, f"{species}", f"{assembly}_{release}_{seqtype}", "*fa"))

            print(f"bowtie-build {params.params} {params.fasta_files} {bowtie1_indices_path}/{species}/{assembly}_{release}_{seqtype}/{species}")
            shell(
                """
                module load apps/bowtie
                mkdir -p {bowtie1_indices_path}/{species}/{assembly}_{release}_{seqtype}
                bowtie-build {params.params} {params.fasta_files} {bowtie1_indices_path}/{species}/{assembly}_{release}_{seqtype}/{species}
                module unload apps/bowtie
                """
                )