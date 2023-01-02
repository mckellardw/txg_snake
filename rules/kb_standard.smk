#############################################
## Pseudoalignment & counting
#############################################

# Make output directory, align fastqs, and generate raw/filtered feature/cell-barcode matrices
#   Info on kallisto/bustools: https://www.kallistobus.tools/
#TODO gzip outputs w/ pigz
rule kb_wrapper:
    input:
        FINAL_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz',
        FINAL_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz'
    output:
        COUNTMAT = '{OUTDIR}/{sample}/kb_wrapper/counts_unfiltered/adata.h5ad'
    params:
        OUTDIR = config['OUTDIR'],
        # KB_EXEC = config['KB_EXEC'], #TODO
        # KB_IDX = config['KB_IDX'],
        # KB_T2G = config['KB_T2G'],
        MEMLIMIT = config['MEMLIMIT']
    log:
        '{OUTDIR}/{sample}/kb_wrapper/kb_wrapper.log'
    threads:
        config['CORES']
    run:
        tmp_chemistry = CHEM_DICT[wildcards.sample]
        STAR_REF = REF_DICT[wildcards.sample]

        UMIlen = CHEMISTRY_SHEET["STAR.UMIlen"][tmp_chemistry]
        SOLOtype = CHEMISTRY_SHEET["STAR.soloType"][tmp_chemistry]
        CB_WHITELIST = CHEMISTRY_SHEET["whitelist"][tmp_chemistry]

        shell(
            f"""
            mkdir -p {params.OUTDIR}/{wildcards.sample}

            {params.KB_EXEC} count \
            -i {params.KB_IDX} \
            -g {params.KB_T2G} \
            -x 10xv3 \
            -o {params.OUTDIR}/{wildcards.sample}/kb_wrapper/ \
            -t {threads} \
            -m {params.MEMLIMIT} \
            -w {CB_WHITELIST} \
            --mm \
            --filter \
            --h5ad \
            --report \
            --workflow standard \
            {input.FINAL_R1_FQ} {input.FINAL_R2_FQ} > {log}
            """
        )

        #TODO- gzip outputs to save some disk space
        # gzip -qf {params.OUTDIR}/{wildcards.sample}/kb/counts_filtered/*.txt
        # gzip -qf {params.OUTDIR}/{wildcards.sample}/kb/counts_filtered/*.mtx
        #
        # gzip -qf {params.OUTDIR}/{wildcards.sample}/kb/counts_unfiltered/*.txt
        # gzip -qf {params.OUTDIR}/{wildcards.sample}/kb/counts_unfiltered/*.mtx

# kallisto/bustools workflow (fastq to bus/txt)
#TODO- make chemistry passable
rule kallisto_bus_standard:
    input:
        FINAL_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz',
        FINAL_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz'
    output:
        BUSTEXT = temp('{OUTDIR}/{sample}/kb/output.corrected.bus'),
        TRANSCRIPTS = '{OUTDIR}/{sample}/kb/transcripts.txt',
        ECMAP = '{OUTDIR}/{sample}/kb/matrix.ec'
    params:
        OUTDIR = config['OUTDIR'],
        # KALLISTO_EXEC = config['KALLISTO_EXEC'],
        # KB_IDX = config['KB_IDX'],
        MEMLIMIT = config['MEMLIMIT']
    log:
        '{OUTDIR}/{sample}/kb/kallisto_bus_standard.log'
    threads:
        config['CORES']
    run:
        tmp_chemistry = CHEM_DICT[wildcards.sample] #TODO
        STAR_REF = REF_DICT[wildcards.sample]

        UMIlen = CHEMISTRY_SHEET["STAR.UMIlen"][tmp_chemistry]
        SOLOtype = CHEMISTRY_SHEET["STAR.soloType"][tmp_chemistry]
        CB_WHITELIST = CHEMISTRY_SHEET["whitelist"][tmp_chemistry]

        shell(f"""
            bash scripts/kb.sh {params.OUTDIR}/{wildcards.sample}/kb \
            {KB_IDX} \
            {input.CB_WHITELIST} \
            {params.CHEMISTRY} \
            {log} \
            {threads} \
            {params.MEMLIMIT} \
            {input.FINAL_R1_FQ} {input.FINAL_R2_FQ}
            """
        )

#TODO- multimapping or EM?
rule bus2mat_standard:
    input:
        BUS = '{OUTDIR}/{sample}/kb/output.corrected.bus',
        TRANSCRIPTS = '{OUTDIR}/{sample}/kb/transcripts.txt',
        ECMAP = '{OUTDIR}/{sample}/kb/matrix.ec'
    output:
        MAT = '{OUTDIR}/{sample}/kb/counts_unfiltered/output.mtx'
        # EC = '{OUTDIR}/{sample}/kb/counts_unfiltered/output.ec.txt'
    params:
        MATDIR = directory('{OUTDIR}/{sample}/kb/counts_unfiltered'),
        OUTDIR = config['OUTDIR'],
        # KB_T2G = config['KB_T2G']
    threads:
        1
    shell:
        """
        mkdir -p {params.MATDIR}

        {BUST_EXEC} count \
        --output {params.MATDIR}/ \
        --genemap {params.KB_T2G} \
        --ecmap {input.ECMAP} \
        --txnames {input.TRANSCRIPTS} \
        --genecounts \
        --umi-gene \
        --em \
        {input.BUS}
        """

        # --multimapping \
        # """
        # mkdir -p {params.OUTDIR}/{wildcards.sample}/kb/counts_unfiltered
        #
        # python scripts/bus2mat.py \
        # {params.OUTDIR}/{wildcards.sample}/kb \
        # {params.KB_T2G} \
        # {input.BUSTEXT}
        # """
