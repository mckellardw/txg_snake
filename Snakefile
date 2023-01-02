########################################################################################################
# txg_snake
#   Snakemake workflow to align and quantify 10x Chromium/Visium/etc-ium datasets
#   v1.0
#   Written by David McKellar
########################################################################################################

import pandas as pd
import scipy.io
import scipy.sparse

########################################################################################################
# Config file
########################################################################################################
configfile:'config.yaml'
CHEMISTRY_SHEET = pd.read_csv(config["CHEMISTRY_SHEET"], na_filter=False,index_col=0) #"resources/chemistry_sheet.csv"
########################################################################################################
# Directories and locations
########################################################################################################
TMPDIR = config['TMPDIR']
OUTDIR = config['OUTDIR']

shell("mkdir -p {TMPDIR}")
shell("mkdir -p {OUTDIR}")

########################################################################################################
# Variables and references
########################################################################################################
SAMPLE_SHEET = pd.read_csv(config["SAMPLE_SHEET_PATH"], na_filter=False)
SAMPLES = [i for i in list(SAMPLE_SHEET['sampleID']) if i] #grab non-empty entries in SAMPLE_SHEET$sampleID

R1_FQS = dict(zip(SAMPLES, list(SAMPLE_SHEET['fastq_R1'])))
R2_FQS = dict(zip(SAMPLES, list(SAMPLE_SHEET['fastq_R2'])))

########################################################################################################
# Executables
########################################################################################################
STAR_EXEC = config['STAR_EXEC']
KB_EXEC = config['KB_EXEC']
KALLISTO_EXEC = config['KALLISTO_EXEC']
BUST_EXEC = config['BUST_EXEC']
FASTQC_EXEC = config["FASTQC_EXEC"]
CUTADAPT_EXEC = config["CUTADAPT_EXEC"]
SAMTOOLS_EXEC = config["SAMTOOLS_EXEC"]
UMITOOLS_EXEC = config["UMITOOLS_EXEC"]
QUALIMAP_EXEC = config["QUALIMAP_EXEC"]
# MULTIQC_EXEC = config["MULTIQC_EXEC"]
BAM2SPLITBW = config["BAM2SPLITBW"]
FASTX_COLLAPSER = config["FASTX_COLLAPSER"]
BLASTDB = config["BLASTDB"]
PICARD_EXEC = config["PICARD_EXEC"]

########################################################################################################
# Pre-run setup
########################################################################################################
CHEM_DICT = {}
REF_DICT = {}
GTF_DICT = {}
for i in range(0,SAMPLE_SHEET.shape[0]):
    tmp_id = list(SAMPLE_SHEET["sampleID"])[i]
    CHEM_DICT[tmp_id] = list(SAMPLE_SHEET["chemistry"])[i]
    REF_DICT[tmp_id] = list(SAMPLE_SHEET["STAR_ref"])[i]
    GTF_DICT[tmp_id] = list(SAMPLE_SHEET["genes_gtf"])[i]

########################################################################################################
rule all:
    input:
        expand('{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw/matrix.mtx.gz', OUTDIR=config['OUTDIR'], sample=SAMPLES),
        # expand('{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw_feature_bc_matrix_h5.h5', OUTDIR=config['OUTDIR'], sample=SAMPLES),
        # expand('{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam.bai', OUTDIR=config['OUTDIR'], sample=SAMPLES), # umi_tools deduplicated .bam
        # expand('{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out_plus.bw', OUTDIR=config['OUTDIR'], sample=SAMPLES), # strand-split bigWigs
        # expand('{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out_merged.bw', OUTDIR=config['OUTDIR'], sample=SAMPLES), #
        expand('{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam.bai', OUTDIR=config['OUTDIR'], sample=SAMPLES), #non-deduplicated .bam; used for saturation estimation
        # expand('{OUTDIR}/{sample}/kb_wrapper/counts_unfiltered/adata.h5ad', OUTDIR=config['OUTDIR'], sample=SAMPLES),
        # expand('{OUTDIR}/{sample}/kb/counts_unfiltered/output.mtx', OUTDIR=config['OUTDIR'], sample=SAMPLES),
        expand('{OUTDIR}/{sample}/preTrim_fastqc_R2_out', OUTDIR=config['OUTDIR'], sample=SAMPLES), # raw R2 fastQC results
        expand('{OUTDIR}/{sample}/postTrim_fastqc_R2_out', OUTDIR=config['OUTDIR'], sample=SAMPLES), # adapter/polyA/ployG-trimmed R2 fastQC results
        expand('{OUTDIR}/{sample}/cutadapt.log', OUTDIR=config['OUTDIR'], sample=SAMPLES),
        expand('{OUTDIR}/{sample}/qualimap_out/qualimapReport.html', OUTDIR=config['OUTDIR'], sample=SAMPLES), # alignment QC qith qualimap plotgardener)
        expand('{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw/spliced.mtx.gz', OUTDIR=config['OUTDIR'], sample=SAMPLES), # count mats
        # expand('{OUTDIR}/{sample}/mirbase/Solo.out/Gene/raw/matrix.mtx', OUTDIR=config['OUTDIR'], sample=SAMPLES), # mirbase count mat
        expand('{OUTDIR}/{sample}/Unmapped_fastqc_out', OUTDIR=config['OUTDIR'], sample=SAMPLES), #fastQC results for unmapped reads
        expand('{OUTDIR}/{sample}/Unmapped.out.mate2_blastResults.txt', OUTDIR=config['OUTDIR'], sample=SAMPLES), # blastn results for unmapped R1 reads non-strand-split bigWigs (for


# Aggregating and QCing raw read data
include: "rules/merge_fqs.smk"
include: "rules/trimQC.smk"

# Alignment, post-alignment clean-up, and QC with STARsolo
# include: "rules/build_ref.smk"
include: "rules/star_align.smk"
include: "rules/umitools_dedup.smk"
include: "rules/star_unmapped.smk"
include: "rules/qualimapQC.smk"

# kallisto
include: "rules/kb_standard.smk"
# include: "rules/kallisto_quant.smk"

include: "rules/bamToSplitBigWig.smk"

#############################################
## Alignment
#############################################


# convert .mtx format to .h5
# Shout out Alex Wolf- https://falexwolf.me/2017/sparse-matrices-with-h5py/
rule STAR_mtx2h5:
    input:
        # VELMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw/spliced.mtx.gz",
        # SJMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/SJ/raw/matrix.mtx.gz",
        GENEMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw/matrix.mtx.gz",
        GENEFULLMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw/matrix.mtx.gz"
    output:
        # VELMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto/raw_feature_bc_matrix_h5.h5",
        # SJMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/SJ/raw_feature_bc_matrix_h5.h5",
        GENEMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw_feature_bc_matrix_h5.h5",
        GENEFULLMAT = "{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull/raw_feature_bc_matrix_h5.h5"
    params:
        # VELDIR = directory("{OUTDIR}/{sample}/STARsolo/Solo.out/Velocyto"),
        # SJDIR = directory("{OUTDIR}/{sample}/STARsolo/Solo.out/SJ"),
        GENEDIR = directory("{OUTDIR}/{sample}/STARsolo/Solo.out/Gene"),
        GENEFULLDIR = directory("{OUTDIR}/{sample}/STARsolo/Solo.out/GeneFull")
    threads:
        1
    run:

        X = scipy.io.mmread(input.GENEFULLMAT)
        bcs = pd.read_csv(
            "features.tsv.gz",
            sep="\t",
            header=None,
            usecols=[1]
        )
        feats = pd.read_csv(
            "features.tsv.gz",
            sep="\t",
            header=None,
            usecols=[1]
        )
        f = h5py.File(output.GENEFULLMAT)
        f.create_dataset('X', data=X)

        f.close()




# # Mark/remove duplicates - **STAR**
# rule star_dedupBAM:
#     input:
#         CB_WHITELIST = config['CB_WHITELIST'],
#         SORTEDBAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.out.bam',
#     output:
#         DEDUPBAM = '{OUTDIR}/{sample}/Processed.out.bam'
#     params:
#         PICARD_EXEC = config['PICARD_EXEC'],
#         OUTDIR = config['OUTDIR'],
#         TMPDIR = '{OUTDIR}/{sample}/_STARtmp',
#         STAR_EXEC = config['STAR_EXEC']
#         # TMPBAM = '{OUTDIR}/{sample}/tmp.bam'
#     threads:
#         config['CORES']*2
#     shell:
#         """
#         {params.STAR_EXEC} \
#         --runMode inputAlignmentsFromBAM \
#         --inputBAMfile {input.SORTEDBAM} \
#         --readFilesType SAM SE \
#         --outFileNamePrefix {params.OUTDIR}/{wildcards.sample}/ \
#         --limitBAMsortRAM 32000000000 \
#         --outTmpDir {params.TMPDIR} \
#         --bamRemoveDuplicatesType UniqueIdentical
#         """
#
# rule star_indexDedupBAM:
#     input:
#         SORTEDBAM = '{OUTDIR}/{sample}/Processed.out.bam'
#     output:
#         BAI = '{OUTDIR}/{sample}/Processed.out.bam.bai'
#     threads:
#         config['CORES']
#     shell:
#         """
#         samtools index -@ {threads} {input.SORTEDBAM}
#         """



# Mark/remove duplicates - **picard**
## https://broadinstitute.github.io/picard/command-line-overview.html#UmiAwareMarkDuplicatesWithMateCigar
# rule dedupBAM:
#     input:
#         CB_WHITELIST = config['CB_WHITELIST'],
#         SORTEDBAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.out.bam'
#     output:
#         DEDUPBAM = '{OUTDIR}/{sample}/Aligned.sortedByCoord.dedup.out.bam',
#         UMI_METRICS = '{OUTDIR}/{sample}/picard/umi_metrics.txt',
#         DUPLICATE_METRICS = '{OUTDIR}/{sample}/picard/duplicate_metrics.txt'
#     params:
#         PICARD_EXEC = config['PICARD_EXEC'],
#         TMPBAM = '{OUTDIR}/{sample}/tmp.bam'
#     log:
#         '{OUTDIR}/{sample}/picard/UmiAwareMarkDuplicatesWithMateCigar.log'
#     threads:
#         config['CORES']
#         #1
#     shell:
#         """
#         samtools view -1 \
#         -@ {threads} \
#         --tag-file CB:{input.CB_WHITELIST} \
#         {input.SORTEDBAM} \
#         > {params.TMPBAM}
#
#         samtools index \
#         -@ {threads} \
#         {params.TMPBAM}
#
#         java -jar {params.PICARD_EXEC} \
#         UmiAwareMarkDuplicatesWithMateCigar \
#         UMI_TAG_NAME=UB \
#         BARCODE_TAG=CB \
#         REMOVE_DUPLICATES=true \
#         ASSUME_SORT_ORDER=coordinate \
#         I={params.TMPBAM} \
#         O={output.DEDUPBAM} \
#         M={output.DUPLICATE_METRICS} \
#         UMI_METRICS={output.UMI_METRICS} \
#         2> {log}
#         """





# cat {input.UNMAPPED1_FQ} | awk '{{if(NR%4==1) {{printf(">%s\n",substr($0,2));}} else if(NR%4==2) print;}}' > {params.TMP_FA}


#############################################
## mirbase alignment
#############################################
# Filter reads by R2 length, so that only short reads are included in miRNA quantification
rule length_filter_R2:
    input:
        FINAL_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz',
        FINAL_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz'
    output:
        FINAL_R1_TRIMMED = '{OUTDIR}/{sample}/tmp/{sample}_R1_final_filtered.fq.gz',
        FINAL_R2_TRIMMED = '{OUTDIR}/{sample}/tmp/{sample}_R2_final_filtered.fq.gz'
    params:
        MIN_LENGTH = 14,
        MAX_LENGTH = 30,
        FIVE_PRIME_R2 = "AAGCAGTGGTATCAACGCAGAGTACATGGG"
    log:
        '{OUTDIR}/{sample}/mirbase/'
    threads:
        config['CORES']
    shell:
        """
        cutadapt \
        --minimum-length {params.MIN_LENGTH} \
        --maximum-length {params.MAX_LENGTH} \
        --pair-filter=first \
 		-a {params.FIVE_PRIME_R2} \
 		-g X{params.FIVE_PRIME_R2} \
 		-o {output.FINAL_R2_TRIMMED} \
        -p {output.FINAL_R1_TRIMMED} \
        --cores {threads} \
        {input.FINAL_R2_FQ} {input.FINAL_R1_FQ} 1> {log}
        """

# Strict STAR alignment with only 1 mismatch allowed and no soft-clipping
# rule mirbase_align:
#     input:
#         CB_WHITELIST = config['CB_WHITELIST'],
#         FINAL_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_final_filtered.fq.gz',
#         FINAL_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final_filtered.fq.gz'
#     output:
#         SORTEDBAM = '{OUTDIR}/{sample}/mirbase/Aligned.sortedByCoord.out.bam', #TODO: add temp()?
#         UNMAPPED1 = '{OUTDIR}/{sample}/mirbase/Unmapped.out.mate1',
#         UNMAPPED2 = '{OUTDIR}/{sample}/mirbase/Unmapped.out.mate2',
#         GENE = directory('{OUTDIR}/{sample}/mirbase/Solo.out/Gene'),
#         GENEMAT = '{OUTDIR}/{sample}/mirbase/Solo.out/Gene/raw/matrix.mtx'
#     params:
#         OUTDIR = config['OUTDIR'],
#         STAR_EXEC = config['STAR_EXEC'],
#         STAR_REF = config['STAR_MIRBASE_REF'],
#         UMIlen = config['UMIlen'],
#         MEMLIMIT = config['MEMLIMIT']
#     threads:
#         config['CORES']
#     shell:
#         """
#         mkdir -p {params.OUTDIR}/{wildcards.sample}/mirbase/
#
#         {params.STAR_EXEC} \
#         --runThreadN {threads} \
#         --outFileNamePrefix {params.OUTDIR}/{wildcards.sample}/mirbase/ \
#         --outSAMtype BAM SortedByCoordinate \
#         --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
#         --readFilesCommand zcat \
#         --soloUMIlen {params.UMIlen} \
#         --genomeDir {params.STAR_REF} \
#         --genomeLoad LoadAndKeep \
#         --limitBAMsortRAM={params.MEMLIMIT} \
#         --readFilesIn {input.FINAL_R2_FQ} {input.FINAL_R1_FQ} \
#         --clipAdapterType CellRanger4 \
#         --outReadsUnmapped Fastx \
#         --soloType CB_UMI_Simple \
#         --soloBarcodeReadLength 0 \
#         --soloCBwhitelist {input.CB_WHITELIST} \
#         --soloCellFilter EmptyDrops_CR \
#         --soloFeatures Gene GeneFull \
#         --soloMultiMappers EM \
#         --alignEndsType EndToEnd \
#         --outFilterMismatchNoverReadLmax 0.2 \
#         --outFilterMultimapNmax 10 \
#         --outFilterScoreMinOverLread 0 \
#         --outFilterMatchNminOverLread 0 \
#         --alignSJDBoverhangMin 1000 \
#         --alignIntronMax 1
#         """
        # --outFilterMultimapScoreRange 0 \
        # --outFilterMatchNmin 16 \

        # gzip -qf {output.GENE}/raw/*
        # gzip -qf {output.GENE}/filtered/*

#############################################
## miRge3.0 analysis
#############################################
#TODO- conda env mismatch issue
# rule miRge3:
#     input:
#         FINAL_R1_FQ = '{OUTDIR}/{sample}/tmp/merged_R1_trimmed.fastq.gz'
#     output:
#         # MIRGE_HTML = '{OUTDIR}/{SAMPLE}/miRge/annotation.report.html'
#         MIRGE_CHECK= '{OUTDIR}/{sample}/miRge_check.txt'
#     params:
#         OUTDIR = config['OUTDIR'],
#         MIRGE_EXEC = config['MIRGE_EXEC'],
#         MIRGE_LIB = config['MIRGE_LIB'],
#         SPECIES = config['SPECIES'],
#         # UMIlen = config['UMIlen'],
#         MEMLIMIT = config['MEMLIMIT']
#     threads:
#         config['CORES']
#     shell:
#         """
#         {params.MIRGE_EXEC} \
#         -s {input.FINAL_R1_FQ} \
#         -lib {params.MIRGE_LIB} \
#         -on {params.SPECIES} \
#         -db mirbase \
#         -o {params.OUTDIR}/{wildcards.SAMPLE} \
#         --threads {threads} \
#         -gff -nmir -ai && touch {output.MIRGE_CHECK}
#         """
# mkdir -p {params.OUTDIR}/{wildcards.SAMPLE}/miRge
# -a illumina \
# -trf


#############################################
## Additional files for visualization
#############################################
