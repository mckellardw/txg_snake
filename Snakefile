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
# Build dictionaries of chemistries & species to use for alignment
CHEM_DICT = {} # Dictionary of chemistry recipe to use for each sample
REF_DICT = {} # Dictionary of reference genomes to use
GTF_DICT = {} # Dictionary of gene annotations (.gtf format)
IDX_DICT = {} # Dictionary of kallisto indices
T2G_DICT = {} # Dictionary of kallisto transcript-to-gene maps
for i in range(0,SAMPLE_SHEET.shape[0]):
    tmp_id = list(SAMPLE_SHEET["sampleID"])[i]
    CHEM_DICT[tmp_id] = list(SAMPLE_SHEET["chemistry"])[i]
    REF_DICT[tmp_id] = list(SAMPLE_SHEET["STAR_ref"])[i]
    GTF_DICT[tmp_id] = list(SAMPLE_SHEET["genes_gtf"])[i]
    IDX_DICT[tmp_id] = list(SAMPLE_SHEET["kb_idx"])[i]
    T2G_DICT[tmp_id] = list(SAMPLE_SHEET["kb_t2g"])[i]

########################################################################################################
rule all:
    input:
        # expand('{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw_feature_bc_matrix_h5.h5', OUTDIR=config['OUTDIR'], sample=SAMPLES),
        # expand('{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam.bai', OUTDIR=config['OUTDIR'], sample=SAMPLES), # umi_tools deduplicated .bam
        # expand('{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out_plus.bw', OUTDIR=config['OUTDIR'], sample=SAMPLES), # strand-split bigWigs
        # expand('{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out_merged.bw', OUTDIR=config['OUTDIR'], sample=SAMPLES), 
        # expand('{OUTDIR}/{sample}/kb_wrapper/counts_unfiltered/adata.h5ad', OUTDIR=config['OUTDIR'], sample=SAMPLES),
        expand('{OUTDIR}/{sample}/kb/counts_unfiltered/output.mtx.gz', OUTDIR=config['OUTDIR'], sample=SAMPLES),
        expand('{OUTDIR}/{sample}/STARsolo/Aligned.sortedByCoord.out.bam.bai', OUTDIR=config['OUTDIR'], sample=SAMPLES), #non-deduplicated .bam; used for saturation estimation
        expand('{OUTDIR}/{sample}/qualimap_out/qualimapReport.html', OUTDIR=config['OUTDIR'], sample=SAMPLES), # alignment QC qith qualimap 
        expand('{OUTDIR}/{sample}/STARsolo/Solo.out/Gene/raw/matrix.mtx.gz', OUTDIR=config['OUTDIR'], sample=SAMPLES),
        expand('{OUTDIR}/{sample}/Unmapped_fastqc_out', OUTDIR=config['OUTDIR'], sample=SAMPLES), #fastQC results for unmapped reads
        expand('{OUTDIR}/{sample}/Unmapped.out.mate2_blastResults.txt', OUTDIR=config['OUTDIR'], sample=SAMPLES), # blastn results for unmapped R1 reads non-strand-split bigWigs (for
        expand('{OUTDIR}/{sample}/preTrim_fastqc_R1_out', OUTDIR=config['OUTDIR'], sample=SAMPLES), # raw R1 fastQC results
        expand('{OUTDIR}/{sample}/preTrim_fastqc_R2_out', OUTDIR=config['OUTDIR'], sample=SAMPLES), # raw R2 fastQC results
        expand('{OUTDIR}/{sample}/postTrim_fastqc_R2_out', OUTDIR=config['OUTDIR'], sample=SAMPLES), # adapter/polyA/ployG-trimmed R2 fastQC results
        expand('{OUTDIR}/{sample}/cutadapt.log', OUTDIR=config['OUTDIR'], sample=SAMPLES)


# Aggregating and QCing raw read data
include: "rules/1_merge_fqs.smk"
include: "rules/1_trimQC.smk"

# Alignment, post-alignment clean-up, and QC with STARsolo
# include: "rules/build_ref.smk"
include: "rules/2_star_align.smk"
include: "rules/2_umitools_dedup.smk"
include: "rules/2_star_unmapped.smk"
include: "rules/2_qualimapQC.smk"

# kallisto
include: "rules/3_kallisto_align.smk"
# include: "rules/3_kallisto_quant.smk"

include: "rules/bamToSplitBigWig.smk"
