# **txg_snake**
Flexible preprocessing, alignment, QC, and quantification workflow for 10x Genomics data (Chromium, Visium, & STRS)
**David W. McKellar**

***Contributions are welcome!***

## Dependencies & Sources:
- `cutadapt` [v4.1](https://cutadapt.readthedocs.io/en/stable/index.html)
- `fastqc` [v0.11.8](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- `STAR` [v2.7.10b](https://github.com/alexdobin/STAR) # Important!
- `kallisto` [v1.0.7](https://pachterlab.github.io/kallisto/)
- `bustools` [v0.1.0.dev2](https://bustools.github.io/)
- `umi-tools` [v1.1.2](https://umi-tools.readthedocs.io/en/latest/index.html)
- `qualimap` [v2.2.a](http://qualimap.conesalab.org/)
- `vsearch` [v2.17.0](https://github.com/torognes/vsearch)
- `BLAST`


## Format for `sample_sheet`:
|sampleID      |fastq_R1                                                               |fastq_R2                                                              |chemistry|STAR_rRNA_ref                 |STAR_ref                  |genes_gtf                |kb_idx                               |kb_t2g                                      |
|--------------|-----------------------------------------------------------------------|----------------------------------------------------------------------|---------|------------------------------|--------------------------|-------------------------|-------------------------------------|--------------------------------------------|
|sample1       | /path/to/sample1_L001_R1.fastq.gz /path/to/sample1_L002_R1.fastq.gz   | /path/to/sample1_L001_R2.fastq.gz /path/to/sample1_L002_R2.fastq.gz  | Visium    | /path/to/STAR_reference_rRNA | /path/to/STAR_reference  |/path/to/annotations.gtf | /path/to/kallisto/transcriptome.idx | /path/to/kallisto/transcripts_to_genes.txt |
|sample2      | /path/to/sample2_L001_R1.fastq.gz /path/to/sample2_L002_R1.fastq.gz   | /path/to/sample2_L001_R2.fastq.gz /path/to/sample2_L002_R2.fastq.gz  | STRS    | /path/to/STAR_reference_rRNA | /path/to/STAR_reference  |/path/to/annotations.gtf | /path/to/kallisto/transcriptome.idx | /path/to/kallisto/transcripts_to_genes.txt |


## Generating references:
#### rRNA STAR reference for in silico rRNA depletion/quantification
Ribosomal RNA (rRNA) molecules can make alignment/quantification very difficult because of the number of genomic copies of these genes. We added a first-pass-alignment just to rRNA sequences to enable stratified parameterization for these sequences, but maintain the ability to count and analyze them.  

Check out `scripts/GRCm39_GENCODEM31_STAR_rRNA.sh` for an example script showing how to generate a rRNA-only STAR reference using GENCODE annotations.  

#### Genomic STAR reference
This is a typical STAR reference that you would use for any other alignment job. Here is an example code snippet:
```
FASTA_GENOME="/path/to/GENCODE_M31/GRCm39.genome.fa"
GENES_DIR="/path/to//GENCODE_M31/gencode.vM31.annotation.gtf"

OUTDIR="/workdir/dwm269/genomes/mm39_all/STAR_GRCm39_GENCODEM31"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeDir ${OUTDIR} \
--genomeFastaFiles ${FASTA_DIR} \
--sjdbGTFfile ${GENES_DIR} \
--sjdbGTFfeatureExon exon
```
*You can find the reference files on [GENCODE's website](https://www.gencodegenes.org/mouse/)*


## Tree of Outputs:
```
{SAMPLE_ID}/
├── log.cutadapt.json
├── postTrim_fastqc_R2_out
│   ├── {SAMPLE_ID}_R2_final_fastqc.html
│   └── {SAMPLE_ID}_R2_final_fastqc.zip
├── preTrim_fastqc_R1_out
│   ├── {SAMPLE_ID}_R1_fastqc.html
│   └── {SAMPLE_ID}_R1_fastqc.zip
├── preTrim_fastqc_R2_out
│   ├── {SAMPLE_ID}_R2_fastqc.html
│   └── {SAMPLE_ID}_R2_fastqc.zip
├── qualimap_out
│   ├── css
│   │   ├── agogo.css
│   │   ├── ajax-loader.gif
│   │   ├── basic.css
│   │   ├── bgfooter.png
│   │   ├── bgtop.png
│   │   ├── comment-bright.png
│   │   ├── comment-close.png
│   │   ├── comment.png
│   │   ├── doctools.js
│   │   ├── down.png
│   │   ├── down-pressed.png
│   │   ├── file.png
│   │   ├── jquery.js
│   │   ├── minus.png
│   │   ├── plus.png
│   │   ├── pygments.css
│   │   ├── qualimap_logo_small.png
│   │   ├── report.css
│   │   ├── searchtools.js
│   │   ├── underscore.js
│   │   ├── up.png
│   │   ├── up-pressed.png
│   │   └── websupport.js
│   ├── images_qualimapReport
│   │   ├── Coverage Profile Along Genes (High).png
│   │   ├── Coverage Profile Along Genes (Low).png
│   │   ├── Coverage Profile Along Genes (Total).png
│   │   ├── Junction Analysis.png
│   │   ├── Reads Genomic Origin.png
│   │   └── Transcript coverage histogram.png
│   ├── qualimapReport.html
│   ├── raw_data_qualimapReport
│   │   ├── coverage_profile_along_genes_(high).txt
│   │   ├── coverage_profile_along_genes_(low).txt
│   │   └── coverage_profile_along_genes_(total).txt
│   └── rnaseq_qc_results.txt
├── rRNA_filtered_fastqc_out
│   ├── {SAMPLE_ID}_R1_final_filtered_fastqc.html
│   ├── {SAMPLE_ID}_R1_final_filtered_fastqc.zip
│   ├── {SAMPLE_ID}_R2_final_filtered_fastqc.html
│   └── {SAMPLE_ID}_R2_final_filtered_fastqc.zip
├── STARsolo
│   ├── Aligned.sortedByCoord.out.bam
│   ├── Aligned.sortedByCoord.out.bam.bai
│   ├── Log.final.out
│   ├── Log.out
│   ├── Log.progress.out
│   ├── SJ.out.tab
│   ├── Solo.out
│   │   ├── Barcodes.stats
│   │   ├── Gene
│   │   │   ├── Features.stats
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv
│   │   │   │   ├── features.tsv
│   │   │   │   └── matrix.mtx
│   │   │   ├── raw
│   │   │   │   ├── barcodes.tsv
│   │   │   │   ├── features.tsv
│   │   │   │   ├── matrix.mtx
│   │   │   │   └── UniqueAndMult-EM.mtx
│   │   │   ├── Summary.csv
│   │   │   └── UMIperCellSorted.txt
│   │   ├── GeneFull
│   │   │   ├── Features.stats
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv
│   │   │   │   ├── features.tsv
│   │   │   │   └── matrix.mtx
│   │   │   ├── raw
│   │   │   │   ├── barcodes.tsv
│   │   │   │   ├── features.tsv
│   │   │   │   ├── matrix.mtx
│   │   │   │   └── UniqueAndMult-EM.mtx
│   │   │   ├── Summary.csv
│   │   │   └── UMIperCellSorted.txt
│   │   ├── SJ
│   │   │   ├── Features.stats
│   │   │   ├── raw
│   │   │   │   ├── barcodes.tsv
│   │   │   │   ├── features.tsv 
│   │   │   │   └── matrix.mtx
│   │   │   └── Summary.csv
│   │   └── Velocyto
│   │       ├── Features.stats
│   │       ├── filtered
│   │       │   ├── ambiguous.mtx
│   │       │   ├── barcodes.tsv
│   │       │   ├── features.tsv
│   │       │   ├── spliced.mtx
│   │       │   └── unspliced.mtx
│   │       ├── raw
│   │       │   ├── ambiguous.mtx
│   │       │   ├── barcodes.tsv
│   │       │   ├── features.tsv
│   │       │   ├── spliced.mtx
│   │       │   └── unspliced.mtx
│   │       └── Summary.csv
│   ├── Unmapped.out.mate1.fastq.gz
│   └── Unmapped.out.mate2.fastq.gz
├── STARsolo_rRNA
│   ├── Aligned.sortedByCoord.out.bam
│   ├── Aligned.sortedByCoord.out.bam.bai
│   ├── Log.final.out
│   ├── Log.out
│   ├── Log.progress.out
│   ├── SJ.out.tab
│   └── Solo.out
│       ├── Barcodes.stats
│       └── GeneFull
│           ├── Features.stats
│           ├── filtered
│           │   ├── barcodes.tsv
│           │   ├── features.tsv
│           │   └── matrix.mtx
│           ├── raw
│           │   ├── barcodes.tsv
│           │   ├── features.tsv
│           │   ├── matrix.mtx
│           │   └── UniqueAndMult-EM.mtx
│           ├── Summary.csv
│           └── UMIperCellSorted.txt
├── tmp
│   ├── {SAMPLE_ID}_R1_final_filtered.fq.gz
│   ├── {SAMPLE_ID}_R1_final.fq.gz
│   ├── {SAMPLE_ID}_R1.fq.gz
│   ├── {SAMPLE_ID}_R2_final_filtered.fq.gz
│   ├── {SAMPLE_ID}_R2_final.fq.gz
│   └── {SAMPLE_ID}_R2.fq.gz
└── Unmapped_fastqc_out
    ├── Unmapped.out.mate1_fastqc.html
    ├── Unmapped.out.mate1_fastqc.zip
    ├── Unmapped.out.mate2_fastqc.html
    └── Unmapped.out.mate2_fastqc.zip
```