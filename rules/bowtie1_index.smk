# bowtie ##########
rule bowtie1_index:
    input:
#        fasta = os.path.join(fasta_download_path, "{species}", '_'.join(["{assembly}", "{release}", "{seqtype}"]), '.'join("{g_file}", "fa"))
        fasta = os.path.join(annotation_download_path, "{species}", '_'.join([{assembly}, {release}, "annotation"]) "{g_file}")
    
    output:
        [os.path.join(bowtie1_indices_path, "{species}", '_'.join([{assembly}, {release}, {seqtype}]), f"{species}.{i}.ebwt") for i in range(1, 5)] +
        [os.path.join(bowtie1_indices_path, "{species}", '_'.join([{assembly}, {release}, {seqtype}]), f"{species}.rev.{i}.ebwt") for i in range(1, 3)]

    params:
        fasta_files = lambda wildcards, input: ",".join(f for f in input),
        params = config["build_indices"]["bowtie1"]["tool_params"]

    run:
        if config["build_indices"]["bowtie1"]["run"]:
            print("Making the bowtie1 indices")
            shell(
                """
                module load apps/bowtie
                mkdir -p {bowtie1_indices_path}{species}/{assembly}/{release}/{seqtype}/
                bowtie-build {params.params} {params.fasta_files} {bowtie1_indices_path}{species}/{assembly}_{release}_{seqtype}/{species}
                module unload apps/bowtie
                """
                )
