# bed file
rule create_bed_file:
    input:
        gtf_file = os.path.join(annotation_download_path, "{s}", "{a}_{r}_annotation", "annotation.gtf")
    
    output:
        bed_file = os.path.join(annotation_download_path, "{s}", "{a}_{r}_annotation", "annotation.bed")

    shell:
        """
        /share/apps/ucsctools/416/gtfToGenePred {input.gtf_file} stdout | /share/apps/ucsctools/416/genePredToBed stdin stdout > {output.bed_file}
        """
