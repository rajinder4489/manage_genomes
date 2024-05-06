import os
import glob

# bowtie2 ##########
rule bowtie2_index:
    input:
        fasta_files = glob.glob(os.path.join(fasta_download_path, f"{species}", f"{assembly}_{release}_{seqtype}", "*fa"))
    
    output:
        indices = [(os.path.join(bowtie2_indices_path, f"{species}", f"{assembly}_{release}_{seqtype}", f"{species}.{i}.bt2") for i in range(1, 5)) if config["build_indices"]["bowtie2"]["run"] else "."] +
                    [(os.path.join(bowtie2_indices_path, f"{species}", f"{assembly}_{release}_{seqtype}", f"{species}.rev.{i}.bt2") for i in range(1, 3)) if config["build_indices"]["bowtie2"]["run"] else "."]

    wildcard_constraints:
        seqtype = "dna|cdna|cds|ncrna"
    
    params:
        fasta_files = lambda wildcards, input: ",".join(f for f in input),
        params = config["build_indices"]["bowtie2"]["tool_params"]

    run:
        if config["build_indices"]["bowtie2"]["run"]:
            print("Making the bowtie2 indices")
            shell(
                """
                module load apps/bowtie2
                mkdir -p {bowtie2_indices_path}/{species}/{assembly}_{release}_{seqtype}
                bowtie2-build {params.params} {params.fasta_files} {bowtie2_indices_path}/{species}/{assembly}_{release}_{seqtype}/{species}
                module unload apps/bowtie2
                """
                )
