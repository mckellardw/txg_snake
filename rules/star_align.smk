# Make output directory, align fastqs, and generate raw/filtered feature/cell-barcode matrices
#   Info for STARsolo command line paramaters: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md

rule STARsolo_align:
    input:
        CB_WHITELIST = config['CB_WHITELIST'],
        FINAL_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz',
        FINAL_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz'
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
        OUTDIR = config['OUTDIR'],
        STAR_EXEC = config['STAR_EXEC'],
        STAR_REF = config['STAR_REF'],
        UMIlen = config['UMIlen'],
        MEMLIMIT = config['MEMLIMIT']
    threads:
        config['CORES']
    shell: #*NOTE* ``--soloBarcodeReadLength 0` should be removed to ensure read-length checking on BC read (R1) is performed
        """
        mkdir -p {params.OUTDIR}/{wildcards.sample}

        {params.STAR_EXEC} \
        --runThreadN {threads} \
        --outFileNamePrefix {params.OUTDIR}/{wildcards.sample}/STARsolo/ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
        --readFilesCommand zcat \
        --soloUMIlen {params.UMIlen} \
        --genomeDir {params.STAR_REF} \
        --genomeLoad LoadAndRemove \
        --limitBAMsortRAM={params.MEMLIMIT} \
        --readFilesIn {input.FINAL_R2_FQ} {input.FINAL_R1_FQ} \
        --clipAdapterType CellRanger4 \
        --outReadsUnmapped Fastx \
        --outFilterMismatchNoverLmax 0.05 \
        --outFilterMatchNmin 12 \
        --outFilterScoreMinOverLread 0 \
        --outFilterMatchNminOverLread 0 \
        --outFilterMultimapNmax 50 \
        --soloType CB_UMI_Simple \
        --soloBarcodeReadLength 0 \
        --soloCBwhitelist {input.CB_WHITELIST} \
        --soloCellFilter EmptyDrops_CR \
        --soloFeatures Gene GeneFull SJ Velocyto \
        --soloMultiMappers EM
        """

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
        1
        # config["CORES_LO"]
    run:
        shell(
            f"""
            gzip -qf {params.VELDIR}/*/*.tsv {params.VELDIR}/*/*.mtx
            gzip -qf {params.GENEDIR}/*/*.tsv {params.GENEDIR}/*/*.mtx
            gzip -qf {params.GENEFULLDIR}/*/*.tsv {params.GENEFULLDIR}/*/*.mtx
            gzip -qf {params.SJDIR}/*/*.tsv {params.SJDIR}/*/*.mtx
            """
        )

rule indexSortedBAM:
    input:
        SORTEDBAM = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam'
    output:
        BAI = '{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam.bai'
    threads:
        config['CORES']
    shell:
        """
        samtools index -@ {threads} {input.SORTEDBAM}
        """
