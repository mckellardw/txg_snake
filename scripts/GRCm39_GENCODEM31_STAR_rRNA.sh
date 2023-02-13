#!/usr/bin/bash

# Script to extract transcript sequences and annotations for rRNAs, and generate a STAR reference

NCORES=24
FASTA_GENOME="/workdir/dwm269/genomes/mm39_all/GENCODE_M31/GRCm39.genome.fa"
FASTA_cDNA="/workdir/dwm269/genomes/mm39_all/GENCODE_M31/gencode.vM31.transcripts.fa.gz"
GENES="/workdir/dwm269/genomes/mm39_all/GENCODE_M31/gencode.vM31.annotation.gtf"

OUTDIR="/workdir/dwm269/genomes/mm39_all/STAR_GRCm39_GENCODEM31_rRNA"
FASTA_rRNA=${OUTDIR}/GRCm39_GENCODEM31_rRNA.fa
GENES_rRNA=${OUTDIR}/GRCm39_GENCODEM31_rRNA.gtf

mkdir -p ${OUTDIR}
cd ${OUTDIR}

# Extract rRNA sequences from cDNA/transcript fasta
zcat ${FASTA_cDNA} | \
 awk -v RS="\n>" -v FS="|" '$8=="rRNA" || $8=="Mt_rRNA" { print ">"$2 $9 }' > \
 ${FASTA_rRNA}

# Build custom gtf from the cDNA/transcript fasta
zcat ${FASTA_cDNA} | \
 awk -v RS="\n>" -v FS="|" '$8=="rRNA" || $8=="Mt_rRNA" {print $2 "\tCUSTOM\texon\t1\t" gsub(/A/,"",$9)+gsub(/C/,"",$9)+gsub(/G/,"",$9)+gsub(/T/,"",$9) "\t.\t+\t.\tgene_id " $2 "; transcript_id " $1 "; gene_type " $8 "; gene_name " $6 "; transcript_type " $8 "; transcript_name " $5 ";\n"}' > \
 ${GENES_rRNA}

# Can't build from the genomic .gtf b/c chromoseome names & start/stop positions are all mismatched
# cat ${GENES} | \
#  grep -E 'gene_type "rRNA"|gene_type "Mt_rRNA"' > \
#  awk -v FS="\t" '#TODO' \
#  ${GENES_rRNA}

STAR \
--runThreadN ${NCORES} \
--runMode genomeGenerate \
--genomeDir ${OUTDIR} \
--genomeFastaFiles ${FASTA_rRNA} \
--sjdbGTFfile ${GENES_rRNA} \
--genomeSAindexNbases 6 \
--sjdbGTFfeatureExon exon
