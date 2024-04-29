# bowtie2 ##########
rule bowtie2_index:
    input:
        fasta = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "{file}"), file = only_fasta)
    
    output:
        [os.path.join(INDICES_BUILD_PATH, "indices", "bowtie2", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.{i}.bt2") for i in range(1, 5)] +
        [os.path.join(INDICES_BUILD_PATH, "indices", "bowtie2", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.rev.{i}.bt2") for i in range(1, 3)]

    params:
        fasta_files = lambda wildcards, input: ",".join(f for f in input),
        params = config["indices"]["bowtie2"]["tool_params"]

    run:
        if config["indices"]["bowtie2"]["run"]:
            print("Making the bowtie2 indices")
            shell(
                """
                module load apps/bowtie2
                mkdir -p {INDICES_BUILD_PATH}/indices/bowtie2/{SPECIES}/{ASSEMBLY}/{RELEASE}//{SEQTYPE}/
                bowtie2-build {params.params} {params.fasta_files} {INDICES_BUILD_PATH}/indices/bowtie2/{SPECIES}/{ASSEMBLY}/{RELEASE}/{SEQTYPE}/{SPECIES}
                module unload apps/bowtie2
                """
                )
