# STAR ##########
rule star_index:
    input:
        fasta_files = os.path.join(fasta_download_path, "{s}", "{a}_{r}_{t}", "genome.fa"),
        gtf_file = os.path.join(annotation_download_path, "{s}", "{a}_{r}_annotation", "gtf.gtf")
        
    output:
        indices = (os.path.join(star_indices_path, "{s}", "{a}_{r}_{t}", "SA")) if config["build_indices"]["star"]["run"] else "."

    wildcard_constraints:
        t = "dna|cdna|cds|ncrna"
    
    params:
        params = config["build_indices"]["star"]["tool_params"]

    log:
        "logs/star_indices/{s}/{a}_{r}_{t}/log.log"

    run:
        if config["build_indices"]["star"]["run"]:
            # print("Making the STAR indices")
            shell(
                """
                (
                    echo "Making the STAR indices"
                    module load apps/star
                    mkdir -p {star_indices_path}/{wildcards.s}/{wildcards.a}_{wildcards.r}_{wildcards.t}
                    STAR --runMode genomeGenerate --genomeDir {star_indices_path}/{wildcards.s}/{wildcards.a}_{wildcards.r}_{wildcards.t} --sjdbGTFfile {input.gtf_file} --genomeFastaFiles {input.fasta_files} {params.params}
                    module unload apps/star
                ) &> {log}
                """
                )
