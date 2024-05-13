# refflat file
rule create_refflat_file:
    input:
        gtf_file = os.path.join(annotation_download_path, "{s}", "{a}_{r}_annotation", "annotation.gtf")
    
    output:
        rf_file = os.path.join(annotation_download_path, "{s}", "{a}_{r}_annotation", "annotation.refflat")

    shell:
        """
        /share/apps/ucsctools/416/gtfToGenePred -genePredExt {input.gtf_file} stdout | awk 'BEGIN{{FS=OFS="\\t"}}{{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' > {output.rf_file}
        """
