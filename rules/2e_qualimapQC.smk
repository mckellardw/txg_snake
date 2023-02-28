## qualimap on aligned reads
rule qualimapQC:
    input:
        SORTEDBAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam' #TODO - change to deduped bam
    output:
        qualimapDir = directory('{OUTDIR}/{sample}/qualimap'),
        qualimapReport = '{OUTDIR}/{sample}/qualimap/qualimapReport.html'
    params:
        MEM = "24G"
    threads:
        config['CORES']
    run:
        GENES_GTF = GTF_DICT[wildcards.sample]
        shell(
            f"""
            mkdir -p {output.qualimapDir}
            cd {output.qualimapDir}

            {QUALIMAP_EXEC} rnaseq \
            -bam {input.SORTEDBAM} \
            -gtf {GENES_GTF} \
            --sequencing-protocol strand-specific-forward \
            --sorted \
            --java-mem-size={params.MEM} \
            -outdir {output.qualimapDir} \
            -outformat html
            """
        )
        # -nt {threads} \
