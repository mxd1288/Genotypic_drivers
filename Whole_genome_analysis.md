# Whole Genome Analysis

# File pre-processing

1. Run FASTQC on raw sequences
2. Make a list of sample IDs: ls *.fastq.gz | grep -oe '.*_R' | uniq | wc -l > ids.txt

## Process raw data files

1. Trim (remove Nextera adapters and quality trim)
2. Flash (Combine overlapping reads, use a max overlap of 130bp)
3. Concatenate singletons into one file
4. bwa mem to map paired and single reads 
5. Samtools to convert sams to bams
6. Picard to add read groups (requires metadata file with the RG information, OCNJ_Library_Overview.txt), this step is needed before Picard.
7. Samtools to filter bam files
8. bamutil to softclip overlapping reads
9. Merge single and paired end bam files (get 1 bam per individual)
10. Mark and remove duplicates

## Fst Outlier Analysis

**Goal:** Identify outlier single nucleotide polymorphisms (SNPs) to examine genetic variation among populations. 

**Approach:** 

1. Use Freebayes and ANGSD to call SNPs from low-coverage whole genome sequencing (lcWGS) data
2. Calculate pairwise Fsts among three Oyster Creek (OC) triad populations using ANGSD (requires SAF and SFS calculations)
3. Permute individuals among populations and repeat per SNP Fst calculation 1000x
4. Calculate p-values as the number of permutations within a heterozygosity bin exceeding a given SNP empirical Fst value

## **Step 1: SNP calling and filtering**

### A) Identify high confidence SNPs by filtering in freebayes

```bash
#!/bin/bash
#BSUB -P fun_gen_1
#BSUB -J freebayes_par
#BSUB -q bigmem
#BSUB -W 120:00
#BSUB -n 15
#BSUB -e /projects2/rsmas/dcrawford/MKD/Whole_genomes/err_and_out/freebayes_par.err
#BSUB -o /projects2/rsmas/dcrawford/MKD/Whole_genomes/err_and_out/freebayes_par.out

#dependencies:
#freebayes
#vcflib
#python2

#you need to have your reference genome with an index

#set variables
REF=/projects2/rsmas/dcrawford/MKD/Funhe_ref
IN=/projects2/rsmas/dcrawford/MKD/Whole_genomes/processed_bams
OUT=/projects2/rsmas/dcrawford/MKD/Whole_genomes/freebayes_GL_all
SB=~/software/local/freebayes/scripts

$SB/freebayes-parallel <($SB/fasta_generate_regions.py $REF/GCF_011125445.2_MU-UCD_Fhet_4.1_genomic.fna.fai 100000) 15 -f $REF/GCF_011125445.2_MU-UCD_Fhet_4.1_genomic.fna $IN/*.bam > $OUT/freebayes_par_all.vcf
```

### B) Filter freebayes VCF for biallelic sites, 5% minor allele frequency, <5% individual missingness, <10% site missingness using VCFtools.

### C) Get genotype likelihoods for freebayes identified sites using ANGSD.

```bash
#!/bin/bash
#BSUB -P fun_gen_1
#BSUB -J angsd_GL
#BSUB -q bigmem
#BSUB -W 120:00
#BSUB -n 39
#BSUB -e /projects2/rsmas/dcrawford/MKD/Whole_genomes/err_and_out/angsd_GL_sites.err
#BSUB -o /projects2/rsmas/dcrawford/MKD/Whole_genomes/err_and_out/angsd_GL_sites.out

BAMS_dir=/projects2/rsmas/dcrawford/MKD/Whole_genomes/processed_bams
BAMS=/projects2/rsmas/dcrawford/MKD/Whole_genomes/processed_bams/bams.list
ANGSD=/home/mxd1288/software/angsd/angsd
OUT=/projects2/rsmas/dcrawford/MKD/Whole_genomes/angsd_all
REF=/projects2/rsmas/dcrawford/MKD/Funhe_ref/GCF_011125445.2_MU-UCD_Fhet_4.1_genomic.fna

# load env
module load anaconda3/biohpc

### First run angsd_depth.sh to get depth statistics and use that to set parameters ###
# don't use any of the set arguments -- just want a GL value for each variant
MINDP=30        # Minimum depth filter
MAXDP=5000        # Maximum depth filter
MININD=30         # Minimum individual filter
MINQ='20'       # Minimum quality filter
MINMAF='0.05'   # Minimum minor allele frequency filter
MINMAPQ='20'    # Minimum mapping quality (alignment score) filter, default value is 20
EXTRA_ARG='-remove_bads 1 -only_proper_pairs 1 -C 50' # Extra arguments when running ANGSD, default value is '-remove_bads 1 -only_proper_pairs 1 -C 50'
SITES=/projects2/rsmas/dcrawford/MKD/Whole_genomes/angsd_all/out.kept.sites #list from VCFtools --kept-sites

# PULL OUT ONLY SPECIFIC SITES - IDENTIFIED WITH FREEBAYES
$ANGSD sites index $SITES

cd $BAMS_dir
$ANGSD -bam $BAMS -ref $REF -out $OUT/ocnj_sites -sites $SITES -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -nThreads 39
```

## **Step 2: Pairwise Fst Calculations**

### A) Calculate Site Allele Frequency (-doSAF) and Minor Allele Frequency (-doMAF) for each group to be compared

```bash
#!/bin/bash
#BSUB -P fun_gen_1
#BSUB -J S19_saf_maf
#BSUB -q normal
#BSUB -W 72:00
#BSUB -n 10
#BSUB -e /projects2/rsmas/dcrawford/MKD/Whole_genomes/err_and_out/S19_pop_sfs_maf.err
#BSUB -o /projects2/rsmas/dcrawford/MKD/Whole_genomes/err_and_out/S19_pop_sfs_maf.out

# load environment
module load anaconda3/biohpc

# set variables
BAMS_dir=/projects2/rsmas/dcrawford/MKD/Whole_genomes/processed_bams
TE_BAMS=/projects2/rsmas/dcrawford/MKD/Whole_genomes/processed_bams/S19_TE_bams.list
S_BAMS=/projects2/rsmas/dcrawford/MKD/Whole_genomes/processed_bams/S19_S_bams.list
N_BAMS=/projects2/rsmas/dcrawford/MKD/Whole_genomes/processed_bams/S19_N_bams.list
ANGSD=/home/mxd1288/software/angsd/angsd
OUT=/projects2/rsmas/dcrawford/MKD/Whole_genomes/pop_saf_maf/
REF=/projects2/rsmas/dcrawford/MKD/Funhe_ref/GCF_011125445.2_MU-UCD_Fhet_4.1_genomic.fna

# ALREADY INDEXED -- MUST RUN angsd index $SITES first
SITES=/projects2/rsmas/dcrawford/MKD/Whole_genomes/angsd_all/out.kept.sites

# Submit one job for each population and timepoint
# SFS will be used for Fst 
# MAF will be used for pop AF

cd $BAMS_dir

# TE Pop
$ANGSD -bam $TE_BAMS -anc $REF -out $OUT/S19_TE_sites -sites $SITES -doSaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -nThreads 15

# South Pop
$ANGSD -bam $S_BAMS -anc $REF -out $OUT/S19_S_sites -sites $SITES -doSaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -nThreads 15

# North Pop
$ANGSD -bam $N_BAMS -anc $REF -out $OUT/S19_N_sites -sites $SITES -doSaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -nThreads 15
```

### B) Calculate 1D Site Frequency Spectrum (SFS) for each group to be compared, 2D SFS for each pairwise comparison, and Fst from the 2D SFS.

```bash
#!/bin/bash
#BSUB -J f19_sfs
#BSUB -P fun_gen_1
#BSUB -W 120:00
#BSUB -n 12
#BSUB -q normal
#BSUB -o /projects2/rsmas/dcrawford/MKD/Whole_genomes/err_and_out/F19_fst.out
#BSUB -e /projects2/rsmas/dcrawford/MKD/Whole_genomes/err_and_out/F19_fst.err

# Good tutorial: http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/angsd-activity-sfs-fst-pbs/ 

# load env
module load anaconda3

# variables
p1=TE
p2=NR
p3=SR

# one timepoint
POP1=/projects2/rsmas/dcrawford/MKD/Whole_genomes/pop_saf_maf/F19_TE_sites.saf.idx
POP2=/projects2/rsmas/dcrawford/MKD/Whole_genomes/pop_saf_maf/F19_N_sites.saf.idx 
POP3=/projects2/rsmas/dcrawford/MKD/Whole_genomes/pop_saf_maf/F19_S_sites.saf.idx
ANGSD=~/software/angsd/misc/realSFS
OUT=/projects2/rsmas/dcrawford/MKD/Whole_genomes/F19_pop_Fsts

# calculate population specific SFS
$ANGSD $POP1 -P 15 -fold 1 > $OUT/$p1\_out.sfs
$ANGSD $POP2 -P 15 -fold 1 > $OUT/$p2\_out.sfs
$ANGSD $POP3 -P 15 -fold 1 > $OUT/$p3\_out.sfs

# calculate 2d SFS
$ANGSD $POP1 $POP2 > $OUT/$p1\.$p2\.ml 
$ANGSD $POP1 $POP3 > $OUT/$p1\.$p3\.ml
$ANGSD $POP2 $POP3 > $OUT/$p2\.$p3\.ml

# calculting Fsts 
$ANGSD fst index $POP1 $POP2 -sfs $OUT/$p1\.$p2\.ml -fstout /projects2/rsmas/dcrawford/MKD/Whole_genomes/F19_pop_Fsts/$p1\.$p2
$ANGSD fst index $POP1 $POP3 -sfs $OUT/$p1\.$p3\.ml -fstout /projects2/rsmas/dcrawford/MKD/Whole_genomes/F19_pop_Fsts/$p1\.$p3
$ANGSD fst index $POP2 $POP3 -sfs $OUT/$p2\.$p3\.ml -fstout /projects2/rsmas/dcrawford/MKD/Whole_genomes/F19_pop_Fsts/$p2\.$p3

# print Fsts
$ANGSD fst print $OUT/TE_NR.fst.idx > $OUT/TE_NR.fst.txt
$ANGSD fst print $OUT/TE_SR.fst.idx > $OUT/TE_SR.fst.txt
$ANGSD fst print $OUT/NR_SR.fst.idx > $OUT/NR_SR.fst.txt
```

### C) Manually extract Fsts from the output *.fst.idx file (Columns from fst print are: chr, pos, a, b). have to add a+b and then get Fst by a/(a+b)

```bash
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $3+$4}' TE_NR.fst.txt > TE_NR_test1.fst.txt

awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $3/$5}' TE_NR_test1.fst.txt > TE_NR_test2.fst.txt
```

## **Step 3: Permutations to get null distribution of Fsts and p-values**

### A) Randomly assign individuals to populations and redo Fst calculations 1000x for each SNPs (parallelized across chromosomes to speed it up). This has to be done over for populations within each timepoint.

Below example is for a single timepoint (Fall 2019).

```bash
#!/bin/bash
#BSUB -P fun_gen_1
#BSUB -J F19_permut_SAF_SFS_FST
#BSUB -q normal
#BSUB -W 120:00
#BSUB -n 10
#BSUB -e /projects2/rsmas/dcrawford/MKD/Whole_genomes/err_and_out/F19_permut.err
#BSUB -o /projects2/rsmas/dcrawford/MKD/Whole_genomes/err_and_out/F19_permut.out

# load environment
module load anaconda3/biohpc

# set variables
BAMS_dir=/projects2/rsmas/dcrawford/MKD/Whole_genomes/processed_bams
BAMS_all=/projects2/rsmas/dcrawford/MKD/Whole_genomes/processed_bams/F19_bams.list
ANGSD=/home/mxd1288/software/angsd/angsd
mkdir /projects2/rsmas/dcrawford/MKD/Whole_genomes/F19_permut_SAF_SFS_FST 
OUT=/projects2/rsmas/dcrawford/MKD/Whole_genomes/F19_permut_SAF_SFS_FST
REF=/projects2/rsmas/dcrawford/MKD/Funhe_ref/GCF_011125445.2_MU-UCD_Fhet_4.1_genomic.fna

# ALREADY INDEXED -- MUST RUN angsd index $SITES first
SITES=/projects2/rsmas/dcrawford/MKD/Whole_genomes/angsd_all/out.kept.sites

# Make directory to store all the shuffled bam files
mkdir /projects2/rsmas/dcrawford/MKD/Whole_genomes/F19_bam_shuf
BAM_shuf=/projects2/rsmas/dcrawford/MKD/Whole_genomes/F19_bam_shuf
BAMS_dir=/projects2/rsmas/dcrawford/MKD/Whole_genomes/processed_bams

# Generate 1000 permuted lists of indivduals
# only do this once!
for i in {1..1000}
do
 #Shuffle bam list and set as a variable
shuf $BAMS_all > $BAM_shuf/BAMS_shuffled_"$i".list
new_bam_list=$BAM_shuf/BAMS_shuffled_"$i".list

 #Randomly assign individuals into each pop list
 #based on actual sample sizes
sed -n '1,100p' $new_bam_list > $BAM_shuf/TE_bams_"$i".list 
sed -n '101,197p' $new_bam_list > $BAM_shuf/SR_bams_"$i".list
sed -n '198,298p' $new_bam_list > $BAM_shuf/NR_bams_"$i".list

mv $BAM_shuf/* $BAMS_dir/
done

# Index sites files
# only do this once#
for i in {1..24}
do
ANGSD=/home/mxd1288/software/angsd/angsd
CHR_SITES=/projects2/rsmas/dcrawford/MKD/Whole_genomes/angsd_all
$ANGSD sites index $CHR_SITES/chr_"$i"_kept.sites
done

# Submit one job for each population and timepoint IN BATCHES of ~250 jobs
# SAF will be used for pairwise SFS, pairwise SFS will be used for Fst
# Fsts will be used to get null distribution for outlier analysis 
# MAF will be used for pop AF

mkdir /projects2/rsmas/dcrawford/MKD/Whole_genomes/F19_permut_jobs
PERMUT_SCRIPTS_DIR=/projects2/rsmas/dcrawford/MKD/Whole_genomes/F19_permut_jobs

# Submit a job for each permuted set to calculate pop specific SAF
# Edit the first line of the for loop to submit just 250 jobs at a time
for i in {1..250}
do
echo '#!/bin/bash' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '#BSUB -P fun_gen_1' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '#BSUB -J F19_permut'"$i"'' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '#BSUB -e /projects2/rsmas/dcrawford/MKD/Whole_genomes/err_and_out/'"$i"'_permut_F19.err' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '#BSUB -o /projects2/rsmas/dcrawford/MKD/Whole_genomes/err_and_out/'"$i"'_permut_F19.out' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '#BSUB -q normal' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '#BSUB -n 10' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '#BSUB -R "rusage[mem=30000]"' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '#BSUB -W 120:00' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '# set variables'>> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'ANGSD=/home/mxd1288/software/angsd/angsd' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'ANGSD_REALSFS=~/software/angsd/misc/realSFS' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'OUT=/projects2/rsmas/dcrawford/MKD/Whole_genomes/F19_permut_SAF_SFS_FST' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'REF=/projects2/rsmas/dcrawford/MKD/Funhe_ref/GCF_011125445.2_MU-UCD_Fhet_4.1_genomic.fna' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'BAM_shuf=/projects2/rsmas/dcrawford/MKD/Whole_genomes/F19_bam_shuf' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'BAMS_dir=/projects2/rsmas/dcrawford/MKD/Whole_genomes/processed_bams' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'POP_BAMS=/projects2/rsmas/dcrawford/MKD/Whole_genomes/processed_bams/{1}_bams_'"$i"'.list' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'POP_LIST=/projects2/rsmas/dcrawford/MKD/Whole_genomes/pop.list' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'CHR_LIST=/projects2/rsmas/dcrawford/MKD/Whole_genomes/chr.list' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'CHR_SITES=/projects2/rsmas/dcrawford/MKD/Whole_genomes/angsd_all' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'cd $BAMS_dir' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '# individual SAF calculations' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'parallel "$ANGSD -bam $POP_BAMS -anc $REF -out $OUT/F19_pop_{1}_chr_{2}_permut_'"$i"' -sites $CHR_SITES/chr_{2}_kept.sites -rf /projects2/rsmas/dcrawford/MKD/Whole_genomes/chr_list.\{2} -doSaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -nThreads 10" ::: TE NR SR ::: {1..24}' \
>> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '# calculate pairwise SFS' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '# variables' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'pop={1}' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'chr={2}' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh 
echo '# one timepoint' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'POP_SAF=$OUT/F19_pop_{1}_chr_{2}_permut_'"$i"'.saf.idx' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '# calculate population specific SFS' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'parallel "$ANGSD_REALSFS $POP_SAF -P 15 -fold 1 > $OUT/$pop\_$chr\_out_'"$i"'.sfs" ::: TE NR SR ::: {1..24}' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'POP1=$OUT/F19_pop_TE_chr_{1}_permut_'"$i"'.saf.idx' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'POP2=$OUT/F19_pop_NR_chr_{1}_permut_'"$i"'.saf.idx' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'POP3=$OUT/F19_pop_SR_chr_{1}_permut_'"$i"'.saf.idx' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'p1=TE' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'p2=NR' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'p3=SR' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'ml_1=$OUT/$p1\.$p2\_chr_{1}_permut_'"$i"'.ml' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'ml_2=$OUT/$p1\.$p3\_chr_{1}_permut_'"$i"'.ml' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'ml_3=$OUT/$p2\.$p3\_chr_{1}_permut_'"$i"'.ml' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '# calculate 2d SFS' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'parallel "$ANGSD_REALSFS $POP1 $POP2 > $ml_1" ::: {1..24}' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh 
echo 'parallel "$ANGSD_REALSFS $POP1 $POP3 > $ml_2" ::: {1..24}' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'parallel "$ANGSD_REALSFS $POP2 $POP3 > $ml_3" ::: {1..24}' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '# calculting Fsts' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh 
echo 'parallel "$ANGSD_REALSFS fst index $POP1 $POP2 -sfs $ml_1 -fstout $OUT/$p1\.$p2\_chr_{1}_'"$i"'" ::: {1..24}' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'parallel "$ANGSD_REALSFS fst index $POP1 $POP3 -sfs $ml_2 -fstout $OUT/$p1\.$p3\_chr_{1}_'"$i"'" ::: {1..24}' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'parallel "$ANGSD_REALSFS fst index $POP2 $POP3 -sfs $ml_3 -fstout $OUT/$p2\.$p3\_chr_{1}_'"$i"'" ::: {1..24}' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '# print Fsts' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'parallel --link "$ANGSD_REALSFS fst print $OUT/TE.NR_chr_{1}_'"$i"'.fst.idx > $OUT/TE.NR_chr_{1}_'"$i"'.fst.txt" ::: {1..24}' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'parallel --link "$ANGSD_REALSFS fst print $OUT/TE.SR_chr_{1}_'"$i"'.fst.idx > $OUT/TE.SR_chr_{1}_'"$i"'.fst.txt" ::: {1..24}' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'parallel --link "$ANGSD_REALSFS fst print $OUT/NR.SR_chr_{1}_'"$i"'.fst.idx > $OUT/NR.SR_chr_{1}_'"$i"'.fst.txt" ::: {1..24}' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '# SAVE SPACE' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '# delete BEAGLE, SAF, .ml, and intermediate fst files' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo 'rm $OUT/*_permut_'"$i"'.saf.idx $OUT\/*_permut_'"$i"'.beagle.gz $OUT/*_permut_'"$i"'.ml' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
echo '' >> $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh

chmod +x $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
bsub < $PERMUT_SCRIPTS_DIR\/"$i"_permut_F19.sh
done
```

*Note: Only need to do this for 1 population pair, not 3 pairwise (can only use one set because they are not independent).*

*Note: Best to run sets of ~250 scripts so as to not overwhelm the LFS HPC.* 

Note: A good way to check if the permutations are all successful is to run something like:

‘cat /projects2/rsmas/dcrawford/MKD/Whole_genomes/err_and_out/*_permut_F19.out' | grep ‘Successfully completed’ | wc -l

Number of lines should match the number of successful permutation jobs.

### B) Summarize the permutation data into a single file (output of the above script is one file per chromosome per permutation). Also need to calculate Fst (as above) and add heterozygosity for each SNP for later sorting.

```bash
# move into permut directory
cd F19_permut_SAF_SFS_FST/

#0 combine all the fsts from a single chromosome
cat *_chr_1_*.fst.txt > chr1.fst.txt

#1 count number of lines in each chr permut file and save to txt file
parallel 'wc -l *_chr_{1}_*fst.txt > chr{1}_order.txt' ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24

#2 Repeat each file name the correct number of times so it will line up with the concatenated chr specific permut file
#parallel 'awk '{ c=0; while ($1>c++) print $2}' chr{1}_order.txt > chr{1}_list.txt' ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
# parallel not working so just did individual for each chr.
awk '{ c=0; while ($1>c++) print $2 }' chr1_order.txt > chr1_list.txt

# but they all have weird "total" lines at the bottom so also remove those
parallel 'grep -v "total" chr{1}_list.txt > temp_{1} && mv temp_{1} chr{1}_list.txt' ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24

#3 combine fst with filelist for each chromosome
parallel 'paste -d '\t' chr{1}.fst.txt chr{1}_list2.txt > chr{1}_fst_filelist.txt' ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24

#4 concatenate all the chromosome specific files into one master permutation list
cat *_fst_filelist.txt > F19_fst_all_permuts.txt # sub in correct timepoint
```

### C) Prune permutations to just get 1000 independent samples and calculate Fst values for them

```bash
# subsample permutations
cat F18_fst_all_permuts.txt | grep -E "NR.SR_chr" > F18_fst_all_permuts_pruned.txt

# calculate heterozygosity 
zcat F18_maf.mafs.gz | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" 2*$6*(1-$6)}' > F18_heterozygosity.txt

# add heterozygosity to master file
awk 'NR==FNR {h[$1 $2] = $7; next} {print $1,$2,$3,$4,h[$1 $2]}' F18_heterozygosity.txt F18_fst_all_permuts_pruned.txt >  F18_fst_all_permuts_het_combined.txt

# add column for a+b
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $3+$4}' F18_fst_all_permuts_het_combined.txt > temp.txt

# add column for Fst (a/a+b)
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $3/$6}' temp.txt > F18_fst_all_het_CALC.txt
rm temp.txt
```

### D) Preparing empirical Fsts

```bash
# change the csv file to a tab seperated txt file
cat file | tr  ',' '\t' > fixed.txt
awk '{print $2, $3, $4, $5, $6, $7, $8 }' F19_fst_all.txt > F19_fst_all2.txt
rm F19_fst_all.txt
mv F19_fst_all2.txt F19_fst_all.txt

# then use bbedit to get rid of the "" around characters

# combine empirical fst list with per SNP heterozygosity
awk 'NR==FNR {h[$1 $2] = $7; next} {print $1,$2,$3,$4,$5,$6,$7,$8,h[$1 $2]}' F19_heterozygosity.txt F19_fst_all.txt > F19_emp_fst_het.txt
```

### E) Sort the master permutation file (from B above) into bins by heterozygosity level and calculate p-values.

```bash
#!/bin/bash
#BSUB -J fst_pval_calc1_par
#BSUB -q normal
#BSUB -W 120:00
#BSUB -n 15
#BSUB -R "rusage[mem=25000]"
#BSUB -e /projects2/rsmas/dcrawford/MKD/Whole_genomes/err_and_out/fst_pval_calc1_par.err 
#BSUB -o /projects2/rsmas/dcrawford/MKD/Whole_genomes/err_and_out/fst_pval_calc1_par.out

module load  anaconda3

#awk 'NR==FNR {h[$1 $2] = $7; next} {print $1,$2,$3,$4,h[$1 $2]}' F18_heterozygosity.txt F18_fst_all_permuts_pruned.txt >  F18_fst_all_permuts_het_combined.txt
#awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $3+$4}' F18_fst_all_permuts_het_combined.txt > temp.txt
#awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $3/$6}' temp.txt > F18_fst_all_het_CALC.txt
#rm temp.txt

cd /projects2/rsmas/dcrawford/MKD/Whole_genomes/F18_permut_summary

# split the master empirical Fst file by heterozygosity
awk '$8<0.000002' F18_emp_fst_het.txt > fst_emp_1.txt
awk '$8>0.000002 && $8<0.00000399999' F18_emp_fst_het.txt > fst_emp_2.txt
awk '$8>0.00000399999 && $8<0.0413881' F18_emp_fst_het.txt > fst_emp_3.txt
awk '$8>0.0413881 && $8<0.0625677' F18_emp_fst_het.txt > fst_emp_4.txt
awk '$8>0.0625677 && $8<0.0913522' F18_emp_fst_het.txt > fst_emp_5.txt
awk '$8>0.0913522 && $8<0.115206' F18_emp_fst_het.txt > fst_emp_6.txt
awk '$8>0.115206 && $8<0.138622' F18_emp_fst_het.txt > fst_emp_7.txt
awk '$8>0.138622 && $8<0.162181' F18_emp_fst_het.txt > fst_emp_8.txt
awk '$8>0.162181 && $8<0.186773' F18_emp_fst_het.txt > fst_emp_9.txt
awk '$8>0.186773 && $8<0.212903' F18_emp_fst_het.txt > fst_emp_10.txt
awk '$8<0.212903 && $8<0.24124' F18_emp_fst_het.txt > fst_emp_11.txt
awk '$8>0.24124 && $8<0.2722896' F18_emp_fst_het.txt > fst_emp_12.txt
awk '$8>0.2722896 && $8<0.305451' F18_emp_fst_het.txt > fst_emp_13.txt
awk '$8>0.305451 && $8<0.3410517' F18_emp_fst_het.txt > fst_emp_14.txt
awk '$8>0.3410517 && $8<0.3778227' F18_emp_fst_het.txt > fst_emp_15.txt
awk '$8>0.3778227 && $8<0.414117' F18_emp_fst_het.txt > fst_emp_16.txt
awk '$8>0.414117 && $8<0.4476487' F18_emp_fst_het.txt > fst_emp_17.txt
awk '$8>0.4476487 && $8<0.475088' F18_emp_fst_het.txt > fst_emp_18.txt
awk '$8>0.475088 && $8<0.493565' F18_emp_fst_het.txt > fst_emp_19.txt
awk '$8>0.493565' F18_emp_fst_het.txt > fst_emp_20.txt

# split again by PopPair
parallel 'cat fst_emp_{}.txt | grep -e 'NR.SR' > F18_emp_NR.SR_{}.txt ' ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
parallel 'cat fst_emp_{}.txt | grep -e 'TE.NR' > F18_emp_TE.NR_{}.txt ' ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
parallel 'cat fst_emp_{}.txt | grep -e 'TE.SR' > F18_emp_TE.SR_{}.txt ' ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20

# split the master FST_HET_COMB file by heterozygosity
awk '$5<0.000002' F18_fst_all_het_CALC.txt > fst_het_1.txt
awk '$5>0.000002 && $5<0.00000399999' F18_fst_all_het_CALC.txt > fst_het_2.txt
awk '$5>0.00000399999 && $5<0.0413881' F18_fst_all_het_CALC.txt > fst_het_3.txt
awk '$5>0.0413881 && $5<0.0625677' F18_fst_all_het_CALC.txt > fst_het_4.txt
awk '$5>0.0625677 && $5<0.0913522' F18_fst_all_het_CALC.txt > fst_het_5.txt
awk '$5>0.0913522 && $5<0.115206' F18_fst_all_het_CALC.txt > fst_het_6.txt
awk '$5>0.115206  && $5<0.138622' F18_fst_all_het_CALC.txt > fst_het_7.txt
awk '$5>0.138622  && $5<0.162181' F18_fst_all_het_CALC.txt > fst_het_8.txt
awk '$5>0.162181 && $5<0.186773' F18_fst_all_het_CALC.txt > fst_het_9.txt
awk '$5>0.186773 && $5<0.212903' F18_fst_all_het_CALC.txt > fst_het_10.txt
awk '$5<0.212903 && $5<0.24124' F18_fst_all_het_CALC.txt > fst_het_11.txt
awk '$5>0.24124  && $5<0.2722896' F18_fst_all_het_CALC.txt > fst_het_12.txt
awk '$5>0.2722896 && $5<0.305451' F18_fst_all_het_CALC.txt > fst_het_13.txt
awk '$5>0.305451 && $5<0.3410517' F18_fst_all_het_CALC.txt > fst_het_14.txt
awk '$5>0.3410517 && $5<0.3778227' F18_fst_all_het_CALC.txt > fst_het_15.txt
awk '$5>0.3778227 && $5<0.414117' F18_fst_all_het_CALC.txt > fst_het_16.txt
awk '$5>0.414117 && $5<0.4476487' F18_fst_all_het_CALC.txt > fst_het_17.txt
awk '$5>0.4476487 && $5<0.475088' F18_fst_all_het_CALC.txt > fst_het_18.txt
awk '$5>0.475088 && $5<0.493565' F18_fst_all_het_CALC.txt > fst_het_19.txt
awk '$5>0.493565' F18_fst_all_het_CALC.txt > fst_het_20.txt

## First had to prune down to keep just one Fst per SNP for each permutation
# cat F18_fst_all_permuts.txt | grep -E "NR.SR_chr" > F18_fst_all_permuts_pruned.txt

FILES=/projects2/rsmas/dcrawford/MKD/Whole_genomes/F18_permut_summary/fst_het_*.txt

for file in $FILES
do
echo "#!/bin/bash" >> "$file".sh
echo "#BSUB -J "$file"" >> "$file".sh
echo "#BSUB -q normal" >> "$file".sh
echo "#BSUB -W 120:00" >> "$file".sh
echo "#BSUB -n 15" >> "$file".sh
echo "#BSUB -R "rusage[mem=25000]"" >> "$file".sh
echo "#BSUB -e "$file".err" >> "$file".sh
echo "#BSUB -o "$file".out" >> "$file".sh

# make a function
echo 'get_pval() {
i=$1;

# HET FILE
echo "Het File = '$file'"
EMP_HET=/projects2/rsmas/dcrawford/MKD/Whole_genomes/F18_permut_summary/F18_emp_het_NR.SR_1.txt 
PVAL_FILE=/projects2/rsmas/dcrawford/MKD/Whole_genomes/F18_permut_summary/F18_snp_pvals.txt
snp_chromosome=$(awk -v row="$i" '"'NR==row+1 {print "'$1'"}'"' $EMP_HET)
snp_pos=$(awk -v row="$i" '"'NR==row+1 {print "'$2'"}'"' $EMP_HET)
snp_fst=$(awk -v row="$i" '"'NR==row+1 {print "'$6'"}'"' $EMP_HET)
snp_poppair=$(awk -v row="$i" '"'NR==row+1 {print "'$7'"}'"' $EMP_HET)
snp_het=$(awk -v row="$i" '"'NR==row+1 {print "'$8'"}'"' $EMP_HET)

echo "SNP chr = $snp_chromosome"
echo "SNP position = $snp_pos"
echo "SNP fst = $snp_fst"
echo "SNP population pair = $snp_poppair"
echo "SNP heterozygosity = $snp_het"

# calculate p-values from permutations using binned heterozygosity
awk -v snp_chromosome="$snp_chromosome" -v snp_pos="$snp_pos" -v snp_fst="$snp_fst" -v snp_poppair="$snp_poppair" -v snp_het="$snp_het" '"'{ ++total; if ("'$7'" > snp_fst) numerator++ } END { print snp_chromosome, snp_pos, snp_poppair, snp_het, snp_fst, total+0, numerator+0, numerator/total}'"' '$file' >> "$PVAL_FILE"
}' >> "$file".sh

# this is fastest but it monopolizes the HPC so cant do it
#for i in "${array[@]}"; 
 # do get_pval $i & done
  
echo 'export -f get_pval' >> "$file".sh
echo 'EMP_HET=/projects2/rsmas/dcrawford/MKD/Whole_genomes/F18_permut_summary/F18_emp_het_NR.SR_1.txt' >> "$file".sh
echo 'NUM_COMPARISONS=$(awk '"'END{print NR}'"' $EMP_HET)' >> "$file".sh
echo 'seq 1 1 $NUM_COMPARISONS | xargs -n 1 -P 15 -I starting_i bash -c '"'get_pval starting_i'"'' >> "$file".sh

chmod +x "$file".sh
#bsub < "$file".sh
done
```

- Example of a single script produced from the loop in E above.
    
    ```bash
    #!/bin/bash
    #BSUB -J fst_NR.SR_10
    #BSUB -q normal
    #BSUB -W 120:00
    #BSUB -n 45
    #BSUB -e /projects2/rsmas/dcrawford/MKD/Whole_genomes/F18_permut_summary/fst_NR.SR_10.txt.err
    #BSUB -o /projects2/rsmas/dcrawford/MKD/Whole_genomes/F18_permut_summary/fst_NR.SR_10.txt.out
    get_pval() {
    i=$1;
    
    # HET FILE
    echo "Het File = /projects2/rsmas/dcrawford/MKD/Whole_genomes/F18_permut_summary/fst_het_10.txt"
    EMP_HET=/projects2/rsmas/dcrawford/MKD/Whole_genomes/F18_permut_summary/F18_emp_NR.SR_10.txt 
    PVAL_FILE=/projects2/rsmas/dcrawford/MKD/Whole_genomes/F18_permut_summary/F18_snp_pvals.txt
    snp_chromosome=$(awk -v row="$i" 'NR==row+1 {print $1}' $EMP_HET)
    snp_pos=$(awk -v row="$i" 'NR==row+1 {print $2}' $EMP_HET)
    snp_fst=$(awk -v row="$i" 'NR==row+1 {print $6}' $EMP_HET)
    snp_poppair=$(awk -v row="$i" 'NR==row+1 {print $7}' $EMP_HET)
    snp_het=$(awk -v row="$i" 'NR==row+1 {print $8}' $EMP_HET)
    
    echo "SNP chr = $snp_chromosome"
    echo "SNP position = $snp_pos"
    echo "SNP fst = $snp_fst"
    echo "SNP population pair = $snp_poppair"
    echo "SNP heterozygosity = $snp_het"
    
    # calculate p-values from permutations using binned heterozygosity
    awk -v snp_chromosome="$snp_chromosome" -v snp_pos="$snp_pos" -v snp_fst="$snp_fst" -v snp_poppair="$snp_poppair" -v snp_het="$snp_het" '{ ++total; if ($7 > snp_fst) numerator++ } END { print snp_chromosome, snp_pos, snp_poppair, snp_het, snp_fst, total+0, numerator+0, numerator/total}' /projects2/rsmas/dcrawford/MKD/Whole_genomes/F18_permut_summary/fst_het_10.txt >> "$PVAL_FILE"
    }
    export -f get_pval
    EMP_HET=/projects2/rsmas/dcrawford/MKD/Whole_genomes/F18_permut_summary/F18_emp_NR.SR_10.txt
    NUM_COMPARISONS=$(awk 'END{print NR}' $EMP_HET)
    seq 1 1 $NUM_COMPARISONS | xargs -n 1 -P 45 -I starting_i bash -c 'get_pval starting_i'
    ```