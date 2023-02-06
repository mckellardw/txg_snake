
rule preTrim_FastQC_R2:
    input:
        MERGED_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2.fq.gz'
    output:
        fastqcDir = directory('{OUTDIR}/{sample}/preTrim_fastqc_R2_out'),
        # fastqcReport = ''
    params:
        adapters = config['FASTQC_ADAPTERS']
    threads:
        config['CORES']
        # min([config['CORES'],8]) # 8 core max based on recommendations from trim_galore authors
    shell:
        """
        mkdir -p {output.fastqcDir}

        {FASTQC_EXEC}  \
        --outdir {output.fastqcDir} \
        --threads {threads} \
        -a {params.adapters} \
        {input.MERGED_R2_FQ}
        """

rule preTrim_FastQC_R1:
    input:
        MERGED_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1.fq.gz'
    output:
        fastqcDir = directory('{OUTDIR}/{sample}/preTrim_fastqc_R1_out'),
        # fastqcReport = ''
    params:
        adapters = config['FASTQC_ADAPTERS']
    threads:
        config['CORES']
        # min([config['CORES'],8]) # 8 core max based on recommendations from trim_galore authors
    shell:
        """
        mkdir -p {output.fastqcDir}

        {FASTQC_EXEC} \
        --outdir {output.fastqcDir} \
        --threads {threads} \
        -a {params.adapters} \
        {input.MERGED_R1_FQ}
        """

# TSO, polyA, and polyG trimming
rule cutadapt_R2:
    input:
        MERGED_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1.fq.gz',
        MERGED_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2.fq.gz'
    output:
        FINAL_R1_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz'),
        FINAL_R2_FQ = temp('{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz')
    params:
        CUTADAPT_EXEC = config["CUTADAPT_EXEC"],
        THREE_PRIME_R2_POLYA = "A"*100, # 100 A-mer
        THREE_PRIME_R2_POLYG = "G"*100, # 100 G-mer
        FIVE_PRIME_R2_TSO = "CCCATGTACTCTGCGTTGATACCACTGCTT", #10x TSO sequence
        FIVE_PRIME_R2_rcTSO = "AAGCAGTGGTATCAACGCAGAGTACATGGG" # rev-comp of 10x TSO sequence
        # FIVE_PRIME_R2 = "TTCGTCACCATAGTTGCGTCTCATGTACCC" #rev 10x TSO sequence
    threads:
        config["CORES"]
    log:
        "{OUTDIR}/{sample}/log.cutadapt.json"
    shell:
        """
        {params.CUTADAPT_EXEC} \
        --minimum-length 18 \
        -A {params.THREE_PRIME_R2_POLYA} \
        -A {params.THREE_PRIME_R2_POLYG} \
 		-G {params.FIVE_PRIME_R2_TSO} \
 		-G {params.FIVE_PRIME_R2_rcTSO} \
        --pair-filter=any \
 		-o {output.FINAL_R1_FQ} \
        -p {output.FINAL_R2_FQ} \
        --cores {threads} \
        --json={log} \
        {input.MERGED_R1_FQ} {input.MERGED_R2_FQ}
        """

rule postTrim_FastQC_R2:
    input:
        FINAL_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz'
    output:
        fastqcDir = directory('{OUTDIR}/{sample}/postTrim_fastqc_R2_out'),
        # fastqcReport = ''
    params:
        adapters = config['FASTQC_ADAPTERS']
    threads:
        config['CORES']
        # min([config['CORES'],8]) # 8 core max based on recommendations from trim_galore authors
    shell:
        """
        mkdir -p {output.fastqcDir}

        {FASTQC_EXEC} \
        --outdir {output.fastqcDir} \
        --threads {threads} \
        -a {params.adapters} \
        {input.FINAL_R2_FQ}
        """
