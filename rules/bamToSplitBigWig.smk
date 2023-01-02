
rule bamToSplitBigWig:
    input:
        BAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.bam',
        BAI = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.bam.bai'
    output:
        POS_BW = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out_plus.bw',
        MERGED_BW = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out_merged.bw'
    params:
        BAM2SPLITBW=config['BAM2SPLITBW'],
        OUTPUT_DIR = '{OUTDIR}/{sample}'
    threads:
        config['CORES']
    run:
        tmp_chemistry = CHEM_DICT[wildcards.sample]
        STAR_REF = REF_DICT[wildcards.sample]

        shell(
            f"""
            {params.BAM2SPLITBW} {input.BAM} {threads} {params.OUTPUT_DIR} {STAR_REF}/chrNameLength.txt
            """
        )
