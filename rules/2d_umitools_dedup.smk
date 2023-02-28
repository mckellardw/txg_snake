# Remove reads that don't have a corrected spot/cell barcode with samtools, then remove duplicates w/ **umi-tools**
## High mem usage? Check here! https://umi-tools.readthedocs.io/en/latest/faq.html
rule umitools_dedupBAM:
    input:
        SORTEDBAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.out.bam',
        SORTEDBAMINDEX = '{OUTDIR}/{sample}/Aligned.sortedByCoord.out.bam.bai'
    output:
        DEDUPBAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.bam'
    params:
        OUTPUT_PREFIX='{OUTDIR}/{sample}/umitools_dedup/{sample}'
    threads:
        config['CORES']
    log:
        '{OUTDIR}/{sample}/dedup.log'
    run:
        whitelist = CHEMISTRY_SHEET["whitelist"][CHEM_DICT[wildcards.sample]]

        # shell(f"""
        #     samtools view -1 -b \
        #     -@ {threads} \
        #     --tag-file CB:{CB_WHITELIST} \
        #     {input.SORTEDBAM} \
        #     > {output.TMPBAM}

        #     samtools index \
        #     -@ {threads} \
        #     {output.TMPBAM}

        #     umi_tools dedup \
        #     -I {output.TMPBAM} \
        #     --extract-umi-method=tag \
        #     --umi-tag=UB \
        #     --cell-tag=CB \
        #     --method=unique \
        #     --per-cell \
        #     --unmapped-reads=discard \
        #     --output-stats={params.OUTPUT_PREFIX} \
        #     --log {log} \
        #     -S {output.DEDUPBAM}
        # """
        # )
        shell(
            f"""
            bash scripts/split_dedup.sh {input.SORTEDBAM} {whitelist} {threads} {output.DEDUPBAM} {OUTDIR}/{wildcards.sample}/tmp/dedup | tee {log}
            """
        )

# Index the deduped .bam file
rule umitools_indexDedupBAM:
    input:
        SORTEDBAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.bam'
    output:
        BAI = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.bam.bai'
    threads:
        config['CORES']
    shell:
        """
        samtools index -@ {threads} {input.SORTEDBAM}
        """

# Split .bam file by strand for IGV browsing
rule strand_split_dedup_bam:
    input:
        DEDUPBAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.bam'
    output:
        FWDBAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.fwd.bam',
        REVBAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.rev.bam'
    threads:
        1
    run:
        shell(
            f"""
            {SAMTOOLS_EXEC} view -b -F 0x10 {input.DEDUPBAM} > {output.FWDBAM}
            {SAMTOOLS_EXEC} view -b -f 0x10 {input.DEDUPBAM} > {output.REVBAM}
            """
        )

# Index the split/deduped bam files
rule indexSplitBAMs:
    input:
        FWDBAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.fwd.bam',
        REVBAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.rev.bam'
    output:
        FWDBAI = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.fwd.bam.bai',
        REVBAI = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.rev.bam.bai'
    threads:
        config['CORES']
    shell:
        """
        {SAMTOOLS_EXEC} index -@ {threads} {input.FWDBAM}
        {SAMTOOLS_EXEC} index -@ {threads} {input.REVBAM}
        """