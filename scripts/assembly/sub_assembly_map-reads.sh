#!/bin/bash

############
# MINIMAP2 #
############

module load biology bedtools samtools

minimap2 -ax map-hifi -t 16 --secondary=no -R '@RG\tID:crypul\tSM:pacbio_libraryprep' crypul_v2024.1_rmY.fa m64077_211121_052238.hifi_reads.q20.filt.fq.gz m64204e_211127_230633.hifi_reads.q20.filt.fq.gz | samtools sort -@ 16 -o crypul_v2024.1_rmY.fa.hifireads.sorted.bam

#########
# ARIMA #
#########

# same script as used for HiC scaffolding

#VARIABLES THAT NEED TO BE DEFINED BEFORE EACH RUN#
ml java biology samtools bwa bedtools

WORK_DIR='pulcher/map_hiC'
ARIMA_PATH=$WORK_DIR/mapping_pipeline

SRA1='P22657_101_S95_L004'

REF=$WORK_DIR/ref/crypul_v2024.1_rmY.fa
PREFIX=$WORK_DIR/ref/crypul_V1.0_rmY

IN_DIR=$WORK_DIR/raw
RAW_DIR=$WORK_DIR/raw_bams
FILT_DIR=$WORK_DIR/filt_bams
TMP_DIR=$WORK_DIR/temp_dir
PAIR_DIR=$WORK_DIR/pair_bams
REP_DIR=$WORK_DIR/rep_dir
MERGE_DIR=$WORK_DIR/merge_bams
REP_NUM=1 #number of the technical replicate set e.g. 1
REP_LABEL=$LABEL\_rep$REP_NUM

BWA='/share/software/user/open/bwa/0.7.17/bin/bwa'
PICARD='picard.jar'
SAMTOOLS='/share/software/user/open/samtools/1.16.1/bin/samtools'
BEDTOOLS='/share/software/user/open/bedtools/2.30.0/bin/bedtools'

FAIDX="$REF.fai"

FILTER=$ARIMA_PATH/filter_five_end.pl
COMBINER=$ARIMA_PATH/two_read_bam_combiner.pl
STATS=$ARIMA_PATH/get_stats.pl

REP_NUM=1 #number of the technical replicates (aka lanes that sequenced the same individual) 
REP_LABEL=$LABEL\_rep$REP_NUM
MAPQ_FILTER=10
CPU=48

####START PIPELINE#####

cd $WORK_DIR

echo "### Step 0: Index reference" # Run only once! Skip this step if you have already generated BWA index files
$BWA index -a bwtsw -p $PREFIX $REF
$SAMTOOLS faidx $REF

####LANE 1####

echo "### Step 1.A: FASTQ to BAM (1st)"
$BWA mem -t $CPU $REF $IN_DIR/$SRA1\_R1_001.fastq.gz | $SAMTOOLS view -@ $CPU -Sb - > $RAW_DIR/$SRA1\_1.bam

echo "### Step 1.B: FASTQ to BAM (2nd)"
$BWA mem -t $CPU $REF $IN_DIR/$SRA1\_R2_001.fastq.gz | $SAMTOOLS view -@ $CPU -Sb - > $RAW_DIR/$SRA1\_2.bam

echo "### Step 1.C: Filter 5' end (1st)"
$SAMTOOLS view -h $RAW_DIR/$SRA1\_1.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/$SRA1\_1.bam

echo "### Step 1.D: Filter 5' end (2nd)"
$SAMTOOLS view -h $RAW_DIR/$SRA1\_2.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/$SRA1\_2.bam

echo "### Step 1.E: Pair reads & mapping quality filter"
perl $COMBINER $FILT_DIR/$SRA1\_1.bam $FILT_DIR/$SRA1\_2.bam $SAMTOOLS $MAPQ_FILTER | $SAMTOOLS view -bS -t $FAIDX - | $SAMTOOLS sort -@ $CPU -o $TMP_DIR/$SRA1.bam -

echo "### Step 1.F: Add read group"
java -Xmx4G -Djava.io.tmpdir=temp/ -jar $PICARD AddOrReplaceReadGroups INPUT=$TMP_DIR/$SRA1.bam OUTPUT=$PAIR_DIR/$SRA1.bam ID=$SRA1 LB=$SRA1 SM=$LABEL PL=ILLUMINA PU=none

echo "### Step 2: Mark duplicates"
java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $PICARD MarkDuplicates INPUT=$PAIR_DIR/$SRA1.bam OUTPUT=$REP_DIR/$SRA1.dedup.bam METRICS_FILE=$REP_DIR/metrics.$SRA1.txt TMP_DIR=$TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

$SAMTOOLS index $REP_DIR/$SRA1.dedup.bam

perl $STATS $REP_DIR/$SRA1.dedup.bam > $REP_DIR/$SRA1.dedup.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"
