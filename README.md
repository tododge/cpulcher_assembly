# cpulcher_assembly

![Suboptimal sampling strategy for high-quality genome assembly: A comparative view on the demise of an Extinct-in-the-Wild reptile](/scripts/figures/pulcher_github.png "figure 2")

This workflow contains scripts used in:
\
\
**Suboptimal sampling strategy for high-quality genome assembly: A comparative view on the demise of an Extinct-in-the-Wild reptile**

Authors: Mario Ernst†, Tristram O. Dodge†, Paul Oliver & Mozes P.K. Blom
\
†These authors contributed equally to this study


## Computation

Computation for this project was conducted on Stanford's HPC Sherlock.

## Genome Assembly

### Filter ccs reads
Filter based on read quality with bamtools and adapter content with hifiadapterfilt.
\
\
```sbatch scripts/assembly/sub_assembly_filter-hifi.sh```
### Run assembly and duplicate purging
Assemble genome with hifiasm and duplicate purging with purge_dups.
\
\
```sbatch scripts/assembly/sub_assembly_hifiasm-purgedups.sh```
### Scaffold purged assembly
Scaffold purged C. pulcher contigs with Omni-C data.
\
\
```sbatch scripts/assembly/sub_assembly_scaffold-omniC.sh```
### Assemble mitochondrial genome
Assemble circular mitochondrial genome with mitohifi.
\
\
```sbatch scripts/assembly/sub_assembly_mitochondria.sh```
### QC stats of the pulcher assembly
Get completeness and contiguity metrics with BUSCO and seqkit.
\
\
```sbatch scripts/assembly/sub_assembly_QC.sh```
### Annotate repeats in C. pulcher assembly
Build repeat database with repeatmodeller and mask repeats with repeatmasker.
\
\
```sbatch scripts/assembly/sub_assembly_repeats.sh```
### Scaffold C. egeriae contigs
Scaffold C. egeriae contigs to C. pulcher chroms with RagTag.
\
\
```sbatch scripts/assembly/sub_assembly_scaffold-egeriae.sh```
### Align pulcher and egeriae genomes
Align genomes with minimap2 and mummer4.
\
\
```sbatch scripts/assembly/sub_assembly_align.sh```
### Map reads back to genome
Map hifi reads with minimap2, and map Omni-C reads using arima genomics pipeline.
\
\
```sbatch scripts/assembly/sub_assembly_map-reads.sh```
### Calculate read coverage
Calculate read depth with bedtools.
\
\
```sbatch scripts/assembly/sub_assembly_coverage.sh```

## Popgen

### Call variants
Call variants using GATK (by chromosome, skinks have long chromosomes)
```
GENOME=crypul_v2024.1_rmY.fa

seqkit fx2tab -n -l $GENOME |  awk '$2>10000000' | cut -f1 > ctg.50mb.list

for ctg in $(cat ctg.50mb.list); do touch gatk.$ctg.submission.sh; cat gatk.header.file.48.txt >> gatk.$ctg.submission.sh; echo "#SBATCH --job-name="gatk.$ctg >> gatk.$ctg.submission.sh; echo CHR=$ctg >> gatk.$ctg.submission.sh; cat scripts/popgen/gatk.command.file.HaplotypeCaller_GenotypeGVCFs.txt >> gatk.$ctg.submission.sh; sbatch gatk.$ctg.submission.sh; done
```

### Variant filtration depending on data type

Adjust depth filters for different datasets.
mean depth of C. pulcher autosomes is 26.4. (1.5x=40, 0.66x=17) for hifi
mean depth of C. pulcher autosomes is 73.2. (1.5x=110, 0.5x=37) for omniC
mean depth of C. egeriae autosomes is 42.4. (1.5x=64, 0.5x=21) for hifi

```
for chrom in $(seqkit fx2tab -n -l $GENOME | awk '$2>10000000' | cut -f1); do sbatch scripts/popgen/sub_gatk.filter.hifi.sh $chrom
for chrom in $(seqkit fx2tab -n -l $GENOME | awk '$2>10000000' | cut -f1); do sbatch scripts/popgen/sub_gatk.filter.omniC.sh $chrom
```
### Combine vcfs
Use gatk to get back to a single VCF for each data type.
\
\
```sbatch scripts/popgen/sub_popgen_mergevcf.sh```
### Define callable sites by both omniC and hifi
Use bedtools to find sites covered by both techologies.
\
\
```sbatch scripts/popgen/sub_popgen_callable-sites.sh```
### Population size inference
Infer historical population sizes using PSMC.
\
\
```sbatch scripts/popgen/sub_popgen_psmc.sh```
### Windowed heterozygosity
Calculate heterozygosity with pixy, accounting for missing data.
\
\
```sbatch scripts/popgen/sub_popgen_pixy.sh```
### Detect ROH
Use plink to detect ROH. Custom R script is in figures.
\
\
```sbatch scripts/popgen/sub_popgen_ROH.sh```
