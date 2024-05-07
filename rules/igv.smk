# igv file
rule create_igv:
    input:
        gtf_file = os.path.join(annotation_download_path, "{s}", "{a}_{r}_annotation", "gtf.gtf")
    
    output:
        igv_file = os.path.join(annotation_download_path, "{s}", "{a}_{r}_annotation", "annotation.igv.gtf")
    
    log:
        f"logs/igv.gtf/{species}/{assembly}_{release}_annotation/log.log"

    shell:
        """
        (
            igvtools sort {input.gtf} {output.gtf}
		    igvtools index {output.gtf}) >& {log}
        ) &> {log}
        """