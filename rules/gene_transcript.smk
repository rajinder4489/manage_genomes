# gene-transcript-relation file (kallisto)
rule create_gene_transcript_file:
    input:
        gtf_file = os.path.join(annotation_download_path, "{s}", "{a}_{r}_annotation", "annotation.gtf")

    output:
        gt_file = os.path.join(annotation_download_path, "{s}", "{a}_{r}_annotation", "gene-transcript.txt")

    #log:
    #    "logs/{s}_{a}_{wildcar.r}_annotation/gene_transcript.log"
    
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
            }}' {input.gtf_file} > {output.gt_file}
        """