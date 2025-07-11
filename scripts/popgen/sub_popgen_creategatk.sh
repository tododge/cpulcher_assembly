GENOME=crypul_v2024.1_rmY.fa

seqkit fx2tab -n -l $GENOME |  awk '$2>10000000' | cut -f1 > ctg.50mb.list

for ctg in $(cat ctg.50mb.list); do touch gatk.$ctg.submission.sh; cat gatk.header.file.48.txt >> gatk.$ctg.submission.sh; echo "#SBATCH --job-name="gatk.$ctg >> gatk.$ctg.submission.sh; echo CHR=$ctg >> gatk.$ctg.submission.sh; cat gatk.command.file.HaplotypeCaller_GenotypeGVCFs.txt >> gatk.$ctg.submission.sh; sbatch gatk.$ctg.submission.sh; done
