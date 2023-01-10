# sudo apt-get install fastqc bwa
# sudo apt install samtools bedtools vcftools
# pip install cutadapt

# wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
# unzip snpEff_latest_core.zip

OUTPUTFOLDER=output
DATAFOLDER=data
SAMPLE=510-7-BRCA_S8_L001

# Qualidade das amostras
sudo fastqc $DATAFOLDER/*.fastq.gz -o $OUTPUTFOLDER/fastqc/

# Cortando
cutadapt --minimum-length 75 --maximum-length 75 \
-o $OUTPUTFOLDER/cutadapt/$SAMPLE\_R1_001_cutadapt.fastq \
-p $OUTPUTFOLDER/cutadapt/$SAMPLE\_R2_001_cutadapt.fastq \
$DATAFOLDER/$SAMPLE\_R1_001.fastq.gz \
$DATAFOLDER/$SAMPLE\_R2_001.fastq.gz

# === Summary ===

# Total read pairs processed:             64,276

# == Read fate breakdown ==
# Pairs that were too short:               4,011 (6.2%)
# Pairs that were too long:               60,118 (93.5%)
# Pairs written (passing filters):           147 (0.2%)

# Total basepairs processed:    17,254,588 bp
#   Read 1:     8,626,127 bp
#   Read 2:     8,628,461 bp
# Total written (filtered):         22,050 bp (0.1%)
#   Read 1:        11,025 bp
#   Read 2:        11,025 bp

# Qualidade das amostras após processamento
INPUTFASTQCTRIMMED=output/cutadapt/*.fastq

sudo fastqc $INPUTFASTQCTRIMMED -o $OUTPUTFOLDER/fastqc/

### Decidi usar sem cortar mesmo

# Índice do genoma de ref
samtools faidx data/hg19.fasta

# BWA alinhamento
sudo bwa mem $DATAFOLDER/hg19.fasta \
$DATAFOLDER/$SAMPLE\_R1_001.fastq.gz \
$DATAFOLDER/$SAMPLE\_R2_001.fastq.gz > $OUTPUTFOLDER/bwa/$SAMPLE\-bwa-mem.sam

# [M::bwa_idx_load_from_disk] read 0 ALT contigs
# [M::process] read 74670 sequences (10000105 bp)...
# [M::process] read 53882 sequences (7254483 bp)...
# [M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR): (7, 35547, 0, 6)
# [M::mem_pestat] skip orientation FF as there are not enough pairs
# [M::mem_pestat] analyzing insert size distribution for orientation FR...
# [M::mem_pestat] (25, 50, 75) percentile: (121, 198, 322)
# [M::mem_pestat] low and high boundaries for computing mean and std.dev: (1, 724)
# [M::mem_pestat] mean and std.dev: (232.17, 136.35)
# [M::mem_pestat] low and high boundaries for proper pairs: (1, 925)
# [M::mem_pestat] skip orientation RF as there are not enough pairs
# [M::mem_pestat] skip orientation RR as there are not enough pairs
# [M::mem_process_seqs] Processed 74670 reads in 12.214 CPU sec, 12.120 real sec
# [M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR): (2, 25635, 0, 2)
# [M::mem_pestat] skip orientation FF as there are not enough pairs
# [M::mem_pestat] analyzing insert size distribution for orientation FR...
# [M::mem_pestat] (25, 50, 75) percentile: (124, 202, 327)
# [M::mem_pestat] low and high boundaries for computing mean and std.dev: (1, 733)
# [M::mem_pestat] mean and std.dev: (235.76, 137.18)
# [M::mem_pestat] low and high boundaries for proper pairs: (1, 936)
# [M::mem_pestat] skip orientation RF as there are not enough pairs
# [M::mem_pestat] skip orientation RR as there are not enough pairs
# [M::mem_process_seqs] Processed 53882 reads in 8.631 CPU sec, 8.551 real sec
# [main] Version: 0.7.17-r1188
# [main] CMD: bwa mem bio-env/data/hg19.fasta bio-env/data/510-7-BRCA_S8_L001_R1_001.fastq.gz bio-env/data/510-7-BRCA_S8_L001_R2_001.fastq.gz
# [main] Real time: 25.646 sec; CPU: 24.326 sec

# Conversão para BAM e sorting
samtools fixmate $OUTPUTFOLDER/bwa/$SAMPLE\-bwa-mem.sam $OUTPUTFOLDER/bwa/$SAMPLE\-bwa-mem.bam 
samtools sort -O bam -o $OUTPUTFOLDER/bwa/$SAMPLE\-bwa-mem-sorted.bam \
 $OUTPUTFOLDER/bwa/$SAMPLE\-bwa-mem.bam

# Visualizando
samtools view -H $OUTPUTFOLDER/bwa/$SAMPLE\-bwa-mem-sorted.bam

# Variant calling
samtools index $OUTPUTFOLDER/bwa/$SAMPLE-bwa-mem-sorted.bam 

# Aqui seria filtrando já
freebayes -f $DATAFOLDER/hg19.fasta \
--target $DATAFOLDER/BRCA.bed -F 0.3 -C 15 \
--pooled-continuous $OUTPUTFOLDER/bwa/$SAMPLE\-bwa-mem-sorted.bam \
> $OUTPUTFOLDER/$SAMPLE\-freebayes.vcf

# Aqui sem filtrar, selecionar só depois com vcftools
freebayes -f $DATAFOLDER/hg19.fasta \
--target $DATAFOLDER/BRCA.bed $OUTPUTFOLDER/bwa/$SAMPLE\-bwa-mem-sorted.bam \
> $OUTPUTFOLDER/freebayes/$SAMPLE\-freebayes.vcf

# Selecionar com vcftools
vcftools --vcf $OUTPUTFOLDER/freebayes/$SAMPLE\-freebayes.vcf --minQ 20 \
 --recode --recode-INFO-all --out $OUTPUTFOLDER/freebayes/$SAMPLE\-freebayes-q20.vcf

# Anotação funcional com snpeff
java -jar bio-env/snpEff/snpEff.jar -dataDir $OUTPUTFOLDER/snpeff -formatEff hg19 \
 $OUTPUTFOLDER/snpeff/$SAMPLE\-freebayes.vcf > $OUTPUTFOLDER/snpeff/$SAMPLE\-freebayes-anno.vcf

