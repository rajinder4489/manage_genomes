# bed file
rule create_bed_file:
    input:
        os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{anno_extensionless}.gtf")
    
    #output:
    #    os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{anno_extensionless}.bed")

    shell:
        """
        
        """
