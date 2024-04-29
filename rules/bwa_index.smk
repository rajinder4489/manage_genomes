# bwa ##########
rule bwa_index:
    input:
        fasta = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "{file}"), file = only_fasta)
    
    output:
        os.path.join(INDICES_BUILD_PATH, "indices", "bwa", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, f"{SPECIES}.sa")

    params:
        params = config["indices"]["bwa"]["tool_params"]

    run:
        if config["indices"]["bwa"]["run"]:
            print("Making the bwa indices")
            shell(
                """
                module load apps/bwa
                mkdir -p {INDICES_BUILD_PATH}/indices/bwa/{SPECIES}/{ASSEMBLY}/{RELEASE}/
                bwa index {params.params} -p {INDICES_BUILD_PATH}/indices/bwa/{SPECIES}/{ASSEMBLY}/{RELEASE}/{SEQTYPE}/{SPECIES} <(cat {input})
                module unload apps/bwa
                """
                )
