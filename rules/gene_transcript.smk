# gene-transcript-relation file (kallisto)
rule create_gene_transcript_file:
    input:
        expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}.gtf"), file = anno_extensionless)
    
    output:
        expand(os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{file}_gene-transcript.txt"), file = anno_extensionless)
        #os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{anno_extensionless}_gene-transcript.txt")

    shell:
        """
        awk -F'\t' '$3=="transcript" {{
            split($9, fields, ";");
            for (i in fields) {{
                split(fields[i], kv, " ");
                if (kv[1] == "gene_id") gi=substr(kv[2], 2, length(kv[2])-2);
                else if (kv[1] == "gene_name") gn=substr(kv[2], 2, length(kv[2])-2);
                else if (kv[1] == "transcript_id") ti=substr(kv[2], 2, length(kv[2])-2);
            }}
            print gi"\t"ti"\t"gn;
        }}' {input} > {output}
        """