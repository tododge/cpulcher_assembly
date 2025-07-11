#!/bin/bash

##############################################
# ARIMA GENOMICS MAPPING PIPELINE 02/08/2019 #
#  MODIFIED JULY 2023 BY TRIS FOR CPULCHER   #
##############################################

# >48 hours, 24 threads (should use more for mapping), 160GB memory for yahs/juicer

################ README ######################

###################### THINGS YOU WILL NEED AT START OF EACH RUN ##########################
##	-fasta file in a directory called $WORK_DIR/ref/
##	-fastq.gz files in a directory called $WORK_DIR/raw/
##	-arima genomics scripts in a directory called $ARIMA_PATH. can be found on github.
##		-cd $WORK_DIR
##		-git clone https://github.com/ArimaGenomics/mapping_pipeline
##
##IF USING YAHS, YOUR DESIRED OUTPUT WILL BE A BAM FILE $REP_DIR/$REP_LABEL.sortedName.bam


#VARIABLES THAT NEED TO BE DEFINED BEFORE EACH RUN#
ml java biology samtools bwa bedtools

WORK_DIR='pulcher/yahs'
ARIMA_PATH=$WORK_DIR/mapping_pipeline

SRA1='P22657_101_S95_L004'
LABEL='crypul-M'

REF=$WORK_DIR/ref/crypul-M.q20.default.bp.p_ctg.purged.fa
PREFIX=$WORK_DIR/ref/crypul-M.q20.default.bp.p_ctg.purged.fa

#VARIABLES THAT SHOULD WORK AS IS ON SHERLOCK#

IN_DIR=$WORK_DIR/raw
RAW_DIR=$WORK_DIR/raw_bams
FILT_DIR=$WORK_DIR/filt_bams
TMP_DIR=$WORK_DIR/temp_dir
PAIR_DIR=$WORK_DIR/pair_bams
REP_DIR=$WORK_DIR/rep_dir
MERGE_DIR=$WORK_DIR/merge_bams
REP_NUM=1 #number of the technical replicate set e.g. 1
REP_LABEL=$LABEL\_rep$REP_NUM

#INPUTS_TECH_REPS=('INPUT=pulcher/yahs/pair_bams/DTG-CHI-584.bam' 'INPUT=pulcher/yahs/pair_bams/DTG-bHiC-696.bam')

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
CPU=62


####START PIPELINE#####

cd $WORK_DIR

echo "### Step 0: Check output directories exist & create them as needed"
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR
[ -d $REP_DIR ] || mkdir -p $REP_DIR
[ -d $MERGE_DIR ] || mkdir -p $MERGE_DIR

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


echo "### Step 4: Mark duplicates"
java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $PICARD MarkDuplicates INPUT=$PAIR_DIR/$SRA1.bam OUTPUT=$REP_DIR/$SRA1.dedup.bam METRICS_FILE=$REP_DIR/metrics.$SRA1.txt TMP_DIR=$TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

$SAMTOOLS index $REP_DIR/$SRA1.dedup.bam

perl $STATS $REP_DIR/$SRA1.dedup.bam > $REP_DIR/$SRA1.dedup.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"


####SORT BAM FILE BY NAME FOR YAHS####

echo "### Step 5: Sort bam file by name"

$SAMTOOLS sort -n $REP_DIR/$SRA1.dedup.bam -@ $CPU -o $REP_DIR/$SRA1.dedup.sortedName.bam


####RUN YAHS####

yahs $REF $REP_DIR/$SRA1.dedup.sortedName.bam -q 40 -o crypul-M_HiC --no-contig-ec


$SAMTOOLS faidx crypul-M_HiC_scaffolds_final.fa

cut -f1,2 crypul-M_HiC_scaffolds_final.fa.fai > crypul-M_HiC_scaffolds_final.chrom.sizes

(yahs/juicer pre crypul-M_HiC.bin crypul-M_HiC_scaffolds_final.agp $FAIDX | /usr/bin/sort -k2,2d -k6,6d -T ./ -S64G --parallel=48 | awk 'NF' > crypul-M_HiC_alignments_sorted.txt.part) && (mv crypul-M_HiC_alignments_sorted.txt.part crypul-M_HiC_alignments_sorted.txt)
(java -jar -Xmx64G juicer_tools_1.19.02.jar pre crypul-M_HiC_alignments_sorted.txt crypul-M_HiC_scaffolds_final.hic.part crypul-M_HiC_scaffolds_final.chrom.sizes) && (mv crypul-M_HiC_scaffolds_final.hic.part crypul-M_HiC_scaffolds_final.hic)

yahs/juicer pre -a -o crypul-M_HiC_out_JBAT crypul-M_HiC.bin crypul-M_HiC_scaffolds_final.agp $FAIDX >crypul-M_HiC_out_JBAT.log 2>&1
(java -jar -Xmx64G juicer_tools_1.19.02.jar pre crypul-M_HiC_out_JBAT.txt crypul-M_HiC_out_JBAT.hic.part <(cat crypul-M_HiC_out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv crypul-M_HiC_out_JBAT.hic.part crypul-M_HiC_out_JBAT.hic)

####made several manual edits to genome in juicebox ####
yahs/juicer post -o crypul-M_HiC_out_JBAT.review crypul-M_HiC_out_JBAT.review.assembly crypul-M_HiC_out_JBAT.liftover.agp $REF

cat crypul-M_HiC_out_JBAT.review.FINAL.fa.fai | cut -f1 | tail -n 3 > junk_scaf_list
seqkit grep -f junk_scaf_list -v crypul-M_HiC_out_JBAT.review.FINAL.fa -o crypul-M_HiC_out_JBAT.review.FINAL.edit.fa

#### with centromere locations unknown, scaffold order is arbitrary. Reverse complementing scaffolds 2, 7, 16 for conistency with 2023 MER manuscript.

seqkit grep -f revcomplist crypul-M_HiC_out_JBAT.review.FINAL.edit.fa | seqkit seq -r -p -t DNA -o revcomp.fa
seqkit grep -f revcomplist crypul-M_HiC_out_JBAT.review.FINAL.edit.fa -v -o norevcomp.fa
cat revcomp.fa norevcomp.fa | seqkit sort -n -o crypul-M_HiC_out_JBAT.review.FINAL.edit_revcomp.fa
sed -i 's/scaffold_/chr/g' crypul-M_HiC_out_JBAT.review.FINAL.edit_revcomp.fa
cp crypul-M_HiC_out_JBAT.review.FINAL.edit_revcomp.fa crypul_V1.0.fa
