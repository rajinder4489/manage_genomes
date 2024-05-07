# bed file
rule create_bed_file:
    input:
        gtf_file = glob.glob(os.path.join(annotation_download_path, f"{species}", f"{assembly}_{release}_{seqtype}", "*gtf"))
    
    output:
        bed_file = os.path.join(annotation_download_path, f"{species}", f"{assembly}_{release}_{seqtype}", "anno.bed")

    shell:
        """
        
        """
