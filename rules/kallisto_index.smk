# kallisto ##########
rule kallisto_index:
    input:
        fasta = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "{file}"), file = only_fasta)
            
    output:
        os.path.join(INDICES_BUILD_PATH, "indices", "kallisto", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.idx")

    params:
        params = config["indices"]["kallisto"]["tool_params"]

    run:
        if config["indices"]["kallisto"]["run"]:
            print("Making the kallisto indices")
            shell(
                """
                module load apps/kallisto
                mkdir -p {INDICES_BUILD_PATH}/indices/kallisto/{SPECIES}/{ASSEMBLY}/{RELEASE}/{SEQTYPE}/
                kallisto index -i {INDICES_BUILD_PATH}/indices/kallisto/{SPECIES}/{ASSEMBLY}/{RELEASE}/{SEQTYPE}/{SPECIES}.idx {input.fasta}
                module unload apps/kallisto
                """
                )
