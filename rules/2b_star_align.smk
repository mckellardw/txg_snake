# Make output directory, align fastqs, and generate raw/filtered feature/cell-barcode matrices
#   Info for STARsolo command line paramaters: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md

rule STARsolo_align:
    input:
        R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_final_filtered.fq.gz',
        R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final_filtered.fq.gz',
    output:
        SORTEDBAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam', #TODO: add temp()
        UNMAPPED1 = '{OUTDIR}/{sample}/STARsolo/Unmapped.out.mate1',
        UNMAPPED2 = '{OUTDIR}/{sample}/STARsolo/Unmapped.out.mate2',
        GENE = directory('{OUTDIR}/{sample}/STARsolo/Solo.out/Gene'),
        GENEFULL = directory('{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull'),
        SJ = directory('{OUTDIR}/{sample}/STARsolo/Solo.out/SJ'),
        VEL = directory('{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto'),
        GENEMAT = '{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw/matrix.mtx',
        GENEFULLMAT = '{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/matrix.mtx',
        SJMAT = '{OUTDIR}/{sample}/STARsolo/Solo.out/SJ/raw/matrix.mtx',
        VELMAT = '{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw/spliced.mtx'
    params:
        STAR_EXEC = config['STAR_EXEC'],
        MEMLIMIT = config['MEMLIMIT']
    threads:
        config['CORES']
    run: 
        tmp_chemistry = CHEM_DICT[wildcards.sample]
        STAR_REF = REF_DICT[wildcards.sample]

        UMIlen = CHEMISTRY_SHEET["STAR.UMIlen"][tmp_chemistry]
        SOLOtype = CHEMISTRY_SHEET["STAR.soloType"][tmp_chemistry]
        CB_WHITELIST = CHEMISTRY_SHEET["whitelist"][tmp_chemistry]
        extraSTAR = CHEMISTRY_SHEET["STAR.extra"][tmp_chemistry]

        #TODO- add whitelist check/gunzip 

        shell(
            f"""
            mkdir -p {OUTDIR}/{wildcards.sample}

            {params.STAR_EXEC} \
            --runThreadN {threads} \
            --outFileNamePrefix {OUTDIR}/{wildcards.sample}/STARsolo/ \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
            --readFilesCommand zcat \
            --soloUMIlen {UMIlen} \
            --genomeDir {STAR_REF} \
            --genomeLoad LoadAndRemove \
            --limitBAMsortRAM={params.MEMLIMIT} \
            --readFilesIn {input.R2_FQ} {input.R1_FQ} \
            --clipAdapterType CellRanger4 \
            --outReadsUnmapped Fastx \
            --soloType {SOLOtype}  {extraSTAR} \
            --soloBarcodeReadLength 0 \
            --soloCBwhitelist {CB_WHITELIST} \
            --soloCellFilter EmptyDrops_CR \
            --soloFeatures Gene GeneFull SJ Velocyto \
            --soloMultiMappers EM
            """
        )
            #STRS parameters:
            # --outFilterMismatchNoverLmax 0.05 \
            # --outFilterMatchNmin 12 \
            # --outFilterScoreMinOverLread 0 \
            # --outFilterMatchNminOverLread 0 \

# compress outputs from STAR (count matrices, cell barcodes, and gene lists)
rule compress_STAR_outs:
    input:
        VELMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw/spliced.mtx",
        SJMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/SJ/raw/matrix.mtx",
        GENEMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw/matrix.mtx",
        GENEFULLMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/matrix.mtx"
    output:
        VELMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw/spliced.mtx.gz",
        SJMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/SJ/raw/matrix.mtx.gz",
        GENEMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw/matrix.mtx.gz",
        GENEFULLMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/matrix.mtx.gz"
    params:
        # SOLODIR = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/"
        VELDIR = directory("{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto"),
        SJDIR = directory("{OUTDIR}/{sample}/STARsolo/Solo.out/SJ"),
        GENEDIR = directory("{OUTDIR}/{sample}/STARsolo/Solo.out/Gene"),
        GENEFULLDIR = directory("{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull")
    threads:
        config['CORES']        
    run:
        shell(
            f"""
            pigz -p{threads} {OUTDIR}/{wildcards.sample}/STARsolo/*/*/*.tsv {OUTDIR}/{wildcards.sample}/STARsolo/*/*/*.mtx
            """
        )
# {params.VELDIR}/*/*.tsv {params.VELDIR}/*/*.mtx \
#              {params.GENEDIR}/*/*.tsv {params.GENEDIR}/*/*.mtx \
#              {params.GENEFULLDIR}/*/*.tsv {params.GENEFULLDIR}/*/*.mtx \
#              {params.SJDIR}/*/*.tsv {params.SJDIR}/*/*.mtx

rule indexSortedBAM:
    input:
        SORTEDBAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam'
    output:
        BAI = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam.bai'
    threads:
        config['CORES']
    run:
        shell(
            f"""
            samtools index -@ {threads} {input.SORTEDBAM}
            """
        )
