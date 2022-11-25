# Generate pseudobam for IGV visualization, etc. 
rule kallisto_quant_standard:
    input:
        CB_WHITELIST = config['CB_WHITELIST'],
        FINAL_R1_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R1_final.fq.gz',
        FINAL_R2_FQ = '{OUTDIR}/{sample}/tmp/{sample}_R2_final.fq.gz'
    output:
        GENOMEBAM = '{OUTDIR}/{sample}/kb/quant/pseudoalignments.bam'
    params:
        OUTDIR = config['OUTDIR'],
        KALLISTO_EXEC = config['KALLISTO_EXEC'],
        KB_IDX = config['KB_IDX'],
        MEMLIMIT = config['MEMLIMIT'],
        GENES_GTF = config['GENES_GTF'],
        CHROMOSOMES = config['CHROMOSOMES']
    log:
        '{OUTDIR}/{sample}/kb_standard/kallisto_quant_standard.log'
    threads:
        config['CORES']
    shell:
        """
        mkdir -p {params.OUTDIR}/{wildcards.sample}/kb_standard/quant

        {params.KALLISTO_EXEC} quant \
        -i {params.KB_IDX} \
        -o {params.OUTDIR}/{wildcards.sample}/kb_standard/quant/ \
        -t {threads} \
        --fr-stranded \
        --single \
        -l 85 \
        -s 10 \
        --genomebam \
        --chromosomes {params.CHROMOSOMES} \
        --gtf {params.GENES_GTF} \
        {input.FINAL_R2_FQ}
        """
