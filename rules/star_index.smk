# STAR ##########
rule star_index:
    input:
        fasta = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "{file}"), file = only_fasta),
        gtf_file = expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}.gtf"), file = anno_extensionless)
            
    output:
        os.path.join(INDICES_BUILD_PATH, "indices", "star", SPECIES, ASSEMBLY, RELEASE, SEQTYPE, "SA")

    params:
        fasta_files = lambda wildcards, input: " ".join(f for f in input.fasta),
        params = config["indices"]["star"]["tool_params"]

    run:
        if config["indices"]["star"]["run"]:
            print("Making the STAR indices")
            shell(
                """
                module load apps/star
                mkdir -p {INDICES_BUILD_PATH}/indices/star/{SPECIES}/{ASSEMBLY}/{RELEASE}/{SEQTYPE}
                STAR --runMode genomeGenerate {params.params} --genomeDir {INDICES_BUILD_PATH}/indices/star/{SPECIES}/{ASSEMBLY}/{RELEASE}/{SEQTYPE} --sjdbGTFfile {input.gtf_file} --genomeFastaFiles {params.fasta_files}
                module unload apps/star
                """
                )
