# Trimming for adapters/homopolymers and fastqc before/after

# FastQC on R2 before trimming
rule preTrim_FastQC_R2:
    input:
        MERGED_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2.fq.gz'
    output:
        fastqcDir = directory('{OUTDIR}/{sample}/preTrim_fastqc_R2'),
        # fastqcReport = ''
    params:
        adapters = config['FASTQC_ADAPTERS']
    threads:
        config['CORES']
        # min([config['CORES'],8]) # 8 core max based on recommendations from trim_galore authors
    run:
        shell(
            f"""
            mkdir -p {output.fastqcDir}

            {FASTQC_EXEC}  \
            --outdir {output.fastqcDir} \
            --threads {threads} \
            -a {params.adapters} \
            {input.MERGED_R2_FQ}
            """
        )

# FastQC on R1 before trimming
rule preTrim_FastQC_R1:
    input:
        MERGED_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1.fq.gz'
    output:
        fastqcDir = directory('{OUTDIR}/{sample}/preTrim_fastqc_R1'),
        # fastqcReport = ''
    params:
        adapters = config['FASTQC_ADAPTERS']
    threads:
        config['CORES']
        # min([config['CORES'],8]) # 8 core max based on recommendations from trim_galore authors
    run:
        shell(
            f"""
            mkdir -p {output.fastqcDir}

            {FASTQC_EXEC} \
            --outdir {output.fastqcDir} \
            --threads {threads} \
            -a {params.adapters} \
            {input.MERGED_R1_FQ}
            """
        )

# TSO, polyA, and polyG trimming
rule cutadapt:
    input:
        MERGED_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1.fq.gz',
        MERGED_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2.fq.gz'
    output:
        FINAL_R1_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz'),
        FINAL_R2_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz')
    params:
        HOMOPOLYMER_ERROR_RATE = 0.2,
        THREE_PRIME_R2_POLYA = "A"*100, # 100 A-mer
        THREE_PRIME_R2_POLYG = "G"*100, # 100 G-mer
        THREE_PRIME_R2_POLYT = "T"*100,
        FIVE_PRIME_R2_TSO = "CCCATGTACTCTGCGTTGATACCACTGCTT", #10x TSO sequence
        FIVE_PRIME_R2_rcTSO = "AAGCAGTGGTATCAACGCAGAGTACATGGG" # rev-comp of 10x TSO sequence
        # FIVE_PRIME_R2 = "TTCGTCACCATAGTTGCGTCTCATGTACCC" #rev 10x TSO sequence
    threads:
        config["CORES"]
    log:
        "{OUTDIR}/{sample}/log.cutadapt.json"
    run:
        shell(f"""
            {CUTADAPT_EXEC} \
            --minimum-length 18 \
            -A "{params.THREE_PRIME_R2_POLYA};max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \
            -A "{params.THREE_PRIME_R2_POLYG};max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \
            -A "{params.THREE_PRIME_R2_POLYT};max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \
            -G {params.FIVE_PRIME_R2_TSO} \
            -G {params.FIVE_PRIME_R2_rcTSO} \
            --pair-filter=any \
            -o {output.FINAL_R1_FQ} \
            -p {output.FINAL_R2_FQ} \
            --cores {threads} \
            --json={log} \
            {input.MERGED_R1_FQ} {input.MERGED_R2_FQ}
            """
        )

# FastQC on R1 after trimming
rule postTrim_FastQC_R1:
    input:
        FINAL_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz'
    output:
        fastqcDir = directory('{OUTDIR}/{sample}/postTrim_fastqc_R1'),
        # fastqcReport = ''
    params:
        adapters = config['FASTQC_ADAPTERS']
    threads:
        config['CORES']
        # min([config['CORES'],8]) # 8 core max based on recommendations from trim_galore authors
    run:
        shell(
            f"""
            mkdir -p {output.fastqcDir}

            {FASTQC_EXEC} \
            --outdir {output.fastqcDir} \
            --threads {threads} \
            -a {params.adapters} \
            {input.FINAL_R1_FQ}
            """
        )

# FastQC on R2 after trimming
rule postTrim_FastQC_R2:
    input:
        FINAL_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz'
    output:
        fastqcDir = directory('{OUTDIR}/{sample}/postTrim_fastqc_R2'),
        # fastqcReport = ''
    params:
        adapters = config['FASTQC_ADAPTERS']
    threads:
        config['CORES']
        # min([config['CORES'],8]) # 8 core max based on recommendations from trim_galore authors
    run:
        shell(
            f"""
            mkdir -p {output.fastqcDir}

            {FASTQC_EXEC} \
            --outdir {output.fastqcDir} \
            --threads {threads} \
            -a {params.adapters} \
            {input.FINAL_R2_FQ}
            """
        )
