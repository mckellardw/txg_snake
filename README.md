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
#TODO