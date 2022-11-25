## qualimap on aligned reads
rule qualimapQC:
    input:
        SORTEDBAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam'
    output:
        qualimapDir = directory('{OUTDIR}/{sample}/qualimap_out'),
        fastqcReport = '{OUTDIR}/{sample}/qualimap_out/qualimapReport.html'
    params:
        GENES_GTF = config['GENES_GTF']
    threads:
        config['CORES']
    shell:
        """
        mkdir -p {output.qualimapDir}
        cd {output.qualimapDir}

        qualimap rnaseq \
        -bam {input.SORTEDBAM} \
        -gtf {params.GENES_GTF} \
        --sequencing-protocol strand-specific-forward \
        --sorted \
        --java-mem-size=8G \
        -outdir {output.qualimapDir} \
        -outformat html
        """
        # -nt {threads} \
