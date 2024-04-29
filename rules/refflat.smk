# refflat file
rule create_refflat_file:
    input:
        os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{anno_extensionless}.gtf")
    
    #output:
    #    os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{anno_extensionless}.refflat")

    shell:
        """
        
        """
