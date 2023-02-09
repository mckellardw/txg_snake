# txg_snake
Preprocessing, alignment, QC, and quantification workflow for 10x Genomics data (Chromium, Visium, & STRS)
**David W. McKellar**

#README TODO:
- Required packages & dependencies (add installation via .yml file)
- Write out pipeline details
- Info on sample_sheet format

## Dependencies & Sources:
- `cutadapt` [v4.1]()
- `fastqc` [v0.11.8]()
- `STAR` [v2.7.10b]() # Important!
- `kallisto` [v1.0.7]()
- `bustools` [v0.1.0.dev2]()
- `umi-tools` [v1.1.2]()
- `qualimap` [v2.2.a]()
- `vsearch` [v2.17.0](https://github.com/torognes/vsearch)
- `BLAST`

## Format for `sample_sheet`:
|sampleID      |fastq_R1                                                               |fastq_R2                                                              |chemistry|STAR_rRNA_ref                 |STAR_ref                  |genes_gtf                |kb_idx                               |kb_t2g                                      |
|--------------|-----------------------------------------------------------------------|----------------------------------------------------------------------|---------|------------------------------|--------------------------|-------------------------|-------------------------------------|--------------------------------------------|
|sample1       | /path/to/sample1_L001_R1.fastq.gz /path/to/sample1_L002_R1.fastq.gz   | /path/to/sample1_L001_R2.fastq.gz /path/to/sample1_L002_R2.fastq.gz  | Visium    | /path/to/STAR_reference_rRNA | /path/to/STAR_reference  |/path/to/annotations.gtf | /path/to/kallisto/transcriptome.idx | /path/to/kallisto/transcripts_to_genes.txt |
|sample2      | /path/to/sample2_L001_R1.fastq.gz /path/to/sample2_L002_R1.fastq.gz   | /path/to/sample2_L001_R2.fastq.gz /path/to/sample2_L002_R2.fastq.gz  | STRS    | /path/to/STAR_reference_rRNA | /path/to/STAR_reference  |/path/to/annotations.gtf | /path/to/kallisto/transcriptome.idx | /path/to/kallisto/transcripts_to_genes.txt |

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
│   │   │   │   ├── features.tsv -> /workdir/dwm269/totalRNA/data/align_out/uSTRS_rRNA/{SAMPLE_ID}/STARsolo/SJ.out.tab
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