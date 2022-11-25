
rule bamToSplitBigWig:
    input:
        BAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.bam',
        BAI = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.bam.bai'
    output:
        POS_BW = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out_plus.bw',
        MERGED_BW = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out_merged.bw'
    params:
        BAM2SPLITBW=config['BAM2SPLITBW'],
        STAR_REF = config['STAR_REF'],
        OUTPUT_DIR = '{OUTDIR}/{sample}'
    threads:
        config['CORES']
    shell:
        """
        {params.BAM2SPLITBW} {input.BAM} {threads} {params.OUTPUT_DIR} {STAR_REF}/chrNameLength.txt
        """
