
# igv file
rule create_igv:
    input:
        os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{anno_extensionless}.gtf")
    
    #output:
    #    os.path.join(GENOME_DOWNLOAD_PATH, SPECIES, ASSEMBLY, RELEASE, "annotation", "{anno_extensionless}.some")

    shell:
        """
        
        """