# igv file
rule create_igv:
    input:
        gtf_file = os.path.join(annotation_download_path, "{s}", "{a}_{r}_annotation", "gtf.gtf")
    
    output:
        igv_file = os.path.join(annotation_download_path, "{s}", "{a}_{r}_annotation", "annotation.igv.gtf")
    
    #log:
    #    f"logs/{wildcards.s}_{wildcards.a}_{wildcards.r}_annotation/igv.gtf.log"

    shell:
        """
            module load apps/java/20.0.1
            /home/sequencing/ragu397g/IGV_2.17.4/igvtools sort {input.gtf_file} {output.igv_file}
		    /home/sequencing/ragu397g/IGV_2.17.4/igvtools index {output.igv_file}
            module unload apps/java/20.0.1
        """