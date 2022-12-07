# Association testing markdown

- using ANGSD association
    - Association between genotype likelihood (beagle file) and RNA MEs or quantitative whole animal/tissue level traits
    - Covariates for quantitative whole animal/tissue level traits: acclimation order and acclimation temp
    - Covariates for RNA MEs: acclimation temp
    - Use score test (-doAsso 2) in ANGSD

## Results

### RNA ME Associations

1. doAsso 2 with temperature as covariate and filters -minHigh 4 -minCount 1 (relaxed compared to default)

## Single RNA Associations

Test the top 10 mRNAs from each module (410 for brain,  390 for heart)

- Make the lists of genes: take one of the trait association files and sort by MM for each module column, then filter to just genes in that module and keep the top 10.
- use the associations data prep R script to pull out and write to .txt files the mRNA expression data for each of the genes, must be ordered in same order as the beagle file!
    - This will introduce a bunch of NAs that have to be replaced by -999 for the association test to work in ANGSD (do this using text editor)

```r
## Get phenotype data ready for association
# has to be one pheno file with -999 coding for missing data 
# has to be in same order as beagle genotype likelihood file

setwd("/Users/melissadrown/Documents/School/RSMAS/Research/Spring2021/Writing/Physio_Plasticity_Manuscript/Scripts and Data")
data_wam <- read.csv("wam_plasticity.csv")
data_ctmax <- read.csv("ctmax_plasticity.csv")

# get file order using FishID
samp_order <- read.csv("~/Desktop/wg_analysis/F18_bams.list.txt", sep="\t", header=FALSE)
samp_order$FishID <- sapply(strsplit(samp_order$V1, "-"), function(x){x[1]})
samp_order

# keep just the columns we want for phenotypic associations and covariates
data_all <- read.csv("1_DLC_MKD_Physiol_trait_Oct2020_JMP.csv")
data_all_ctmax <- subset(data_all, data_all$Date_CtMax!="04/15/2019")

data_all <- data_all[,c(2,4,6,8,9,39,40,41,42,43,44)]
head(data_all)

data_wam <- data_wam[,c(2,19,40)]
head(data_wam)

data_ctmax <- data_ctmax[,c(2,16,33)]
head(data_ctmax)

data_all_12 <- subset(data_all, data_all$Acc_Temp=="12")
rownames(data_all_12) <- data_all_12$Fish_ID

data_all_28 <- subset(data_all, data_all$Acc_Temp=="28")
rownames(data_all_28) <- data_all_28$Fish_ID

# Reorder
file_idx <- match(samp_order$FishID, data_all_12$Fish_ID)
data_all_12 <- data_all_12[file_idx,]

file_idx <- match(samp_order$FishID, data_all_28$Fish_ID)
data_all_28 <- data_all_28[file_idx,]

file_idx <- match(samp_order$FishID, data_wam$FishID)
data_wam <- data_wam[file_idx,]

file_idx <- match(samp_order$FishID, data_ctmax$FishID)
data_ctmax <- data_ctmax[file_idx,]

file_idx <- match(samp_order$FishID, data_all$Fish_ID)
data_all <- data_all[file_idx,]

qqnorm(data_all$Log_Log_WAM_Mass_residuals)
qqnorm(data_all_ctmax$CtMax_.C_Residual_WAMMASS)
qqnorm(data_all$FA_CamMass_Residuals)
qqnorm(data_all$LKA__CamMass_Residuals)
qqnorm(data_all$GLU_CamMass_Residuals)
qqnorm(data_all$END_CamMass_Residuals)
qqnorm(data_ctmax$ctmax_log2ratio)
qqnorm(data_wam$log2_ratio_accl)

file_idx <- match(samp_order$FishID, data_all_ctmax$Fish_ID)
data_all_ctmax <- data_all_ctmax[file_idx,]

# had to remove that one day of individuals
data_all_ctmax_28 <- subset(data_all_ctmax, data_all_ctmax$Acc_Temp=="28")
file_idx <- match(samp_order$FishID, data_all_ctmax_28$Fish_ID)
data_all_ctmax_28 <- data_all_ctmax_28[file_idx,]

data_ctmax
#write.table(data_all[,3:4], "~/Desktop/wg_analysis/associations/temp_cov_pheno_data_all.txt")
#write.table(data_all$Log_Log_WAM_Mass_residuals, "~/Desktop/wg_analysis/associations/wam_pheno_data_all.txt")
#write.table(data_all_ctmax$CtMax_.C_Residual_WAMMASS, "~/Desktop/wg_analysis/associations/ctmax_pheno_data_all.txt")
#write.table(data_all$Acc_Temp, "~/Desktop/wg_analysis/associations/cam_temp_pheno_data_all.txt")
#write.table(data_all$FA_CamMass_Residuals, "~/Desktop/wg_analysis/associations/FAcam_pheno_data_all.txt")
#write.table(data_all$GLU_CamMass_Residuals, "~/Desktop/wg_analysis/associations/GLUcam_pheno_data_all.txt")
#write.table(data_all$LKA__CamMass_Residuals, "~/Desktop/wg_analysis/associations/LKAcam_pheno_data_all.txt")
#write.table(data_all$END_CamMass_Residuals, "~/Desktop/wg_analysis/associations/ENDcam_pheno_data_all.txt")

# write out phenotype files for association
#write.table(data_all_12$Log_Log_WAM_Mass_residuals, "~/Desktop/wg_analysis/associations/wam_pheno_data_all_12.txt")
#write.table(data_all_12$CtMax_.C_Residual_WAMMASS, "~/Desktop/wg_analysis/associations/ctmax_pheno_data_all_12.txt")
#write.table(data_all_12$FA_CamMass_Residuals, "~/Desktop/wg_analysis/associations/FAcam_pheno_data_all_12.txt")
#write.table(data_all_12$LKA__CamMass_Residuals, "~/Desktop/wg_analysis/associations/LKAwam_pheno_data_all_12.txt")
#write.table(data_all_12$GLU_CamMass_Residuals, "~/Desktop/wg_analysis/associations/GLUwam_pheno_data_all_12.txt")
#write.table(data_all_12$END_CamMass_Residuals, "~/Desktop/wg_analysis/associations/ENDwam_pheno_data_all_12.txt")

#write.table(data_all_28$Log_Log_WAM_Mass_residuals, "~/Desktop/wg_analysis/associations/wam_pheno_data_all_28.txt")

write.table(data_all_ctmax_28$CtMax_.C_Residual_WAMMASS, "~/Desktop/wg_analysis/associations/ctmax_pheno_data_all_28.txt")
#write.table(data_all_28$FA_CamMass_Residuals, "~/Desktop/wg_analysis/associations/FAcam_pheno_data_all_28.txt")
#write.table(data_all_28$LKA__CamMass_Residuals, "~/Desktop/wg_analysis/associations/LKAwam_pheno_data_all_28.txt")
#write.table(data_all_28$GLU_CamMass_Residuals, "~/Desktop/wg_analysis/associations/GLUwam_pheno_data_all_28.txt")
#write.table(data_all_28$END_CamMass_Residuals, "~/Desktop/wg_analysis/associations/ENDwam_pheno_data_all_28.txt")

#write.table(data_all$SEX, "~/Desktop/wg_analysis/associations/Sex_pheno_data_all.txt")

#write.table(data_all_12$Acc_Order, "~/Desktop/wg_analysis/associations/accl_order_pheno_data_all_12.txt")
#write.table(data_all_28$Acc_Order, "~/Desktop/wg_analysis/associations/accl_order_pheno_data_all_28.txt")

#write.table(data_wam$log2_ratio_accl, "~/Desktop/wg_analysis/associations/pheno_wam_plasticity.txt")
#write.table(data_ctmax$ctmax_log2ratio, "~/Desktop/wg_analysis/associations/pheno_ctmax_plasticity.txt")

#write.table(data_wam$accl_order.x, "~/Desktop/wg_analysis/associations/accl_order_wam_plasticity.txt")
#write.table(data_ctmax$accl_order.x, "~/Desktop/wg_analysis/associations/accl_order_ctmax_plasticity.txt")

##########################################################################################################################################
# Prep RNA ME Data
heart_ME <- read.csv("~/Documents/School/RSMAS/Research/Spring2021/Writing/RNA_seq_F18_OCNJ/OCNJ_F18_RNA/July-21_WGCNA/ME_values/heart_MEs.csv")
brain_ME <- read.csv("~/Documents/School/RSMAS/Research/Spring2021/Writing/RNA_seq_F18_OCNJ/OCNJ_F18_RNA/July-21_WGCNA/ME_values/brain_MEs.csv")

file_idx <- match(samp_order$FishID, heart_ME$X.1)
heart_ME <- heart_ME[file_idx,]

file_idx <- match(samp_order$FishID, brain_ME$X.1)
brain_ME <- brain_ME[file_idx,]

#write.csv(heart_ME, "ordered_heart_MEs.csv")
#write.csv(brain_ME, "ordered_brain_MEs.csv")

########################################################################################################################################
##
# Prep single mRNA data
h_mrna <- read.csv("~/Documents/School/RSMAS/Research/Spring2021/Writing/RNA_seq_F18_OCNJ/OCNJ_F18_RNA/heart_counts_filtered.csv")

b_mrna <- read.csv("~/Documents/School/RSMAS/Research/Spring2021/Writing/RNA_seq_F18_OCNJ/OCNJ_F18_RNA/brain_counts_filtered.csv")

h_list <- read.csv("~/Desktop/wg_analysis/associations/single_heart_mRNAs.txt", header=FALSE)
b_list <- read.csv("~/Desktop/wg_analysis/associations/single_brain_mRNAs.txt", header=FALSE)

h_mrna <- subset(h_mrna, h_mrna$X %in% h_list$V1)
b_mrna <- subset(b_mrna, b_mrna$X %in% b_list$V1)

h_mrna <- t(h_mrna)
b_mrna <- t(b_mrna)

#write.csv(b_mrna, "~/Desktop/wg_analysis/associations/input_data/individual_mRNAs/brain_single_mRNA_expression2.csv")
#write.csv(h_mrna, "~/Desktop/wg_analysis/associations/input_data/individual_mRNAs/heart_single_mRNA_expression2.csv")

b_mrna_samps <- read.csv("~/Desktop/wg_analysis/associations/input_data/individual_mRNAs/brain_single_mRNA_expression2.csv")
h_mrna_samps <- read.csv("~/Desktop/wg_analysis/associations/input_data/individual_mRNAs/heart_single_mRNA_expression2.csv")

nrow(b_mrna_samps)
nrow(h_mrna_samps)

rownames(b_mrna_samps) <- b_mrna_samps$X
rownames(h_mrna_samps) <- h_mrna_samps$X
# reorder and save the new files
file_idx <- match(samp_order$FishID, rownames(b_mrna_samps))
b_ordered <- b_mrna_samps[file_idx,]

file_idx <- match(samp_order$FishID, rownames(h_mrna_samps))
h_ordered <- h_mrna_samps[file_idx,]

for (i in 1:length(colnames(b_ordered))){
  write.table(b_ordered[,i], paste0("~/Desktop/wg_analysis/associations/input_data/individual_mRNAs/",colnames(b_ordered)[i],"_mRNA_brain.txt"), sep="\t", row.names=FALSE)
}

for (i in 1:length(colnames(h_ordered))){
  write.csv(h_ordered[,i], paste0("~/Desktop/wg_analysis/associations/input_data/individual_mRNAs/",colnames(h_ordered)[i],"_mRNA_heart.txt"), sep="\t", row.names=FALSE)
}
```

- doAsso 2 (score test) with temperature as covariate and filters -minHigh 4 -minCount 1 (relaxed compared to default)

## Trait Association Testing

1. doAsso 2 (score test) with acclimation order and temperature as covariates and filtering as -minHigh 4 -minCount 1 (relaxed compared to default)

## Annotating the resulting SNPs that are important

```bash
# take the association list from ANGSD (*.lrt0.gz) and use R to get pvalues and corrected pvalues  
# save the csv list with the singificant SNPs in it from R
# Use the csv file as input for bedtools to merge with gtf file
awk 'BEGIN{FS=OFS=","}{if(NR>1){ print $2, $3-1 ,$3 }}' all_heart_signif.csv > heart_snp_bedfile.txt

# MAKE SURE to change the chromosome names back to NC_ and NW_ before doing the next step
# also the snp bedfile must be tab deliminated
~/software/bedtools2/bin/bedtools intersect -wb -a heart_snp_bedfile.txt -b *.gtf > heart_snp_annot.txt

# and then add a new column that has the file name in it so we know what mRNA was associated with what SNP
awk -F, -v OFS=, '
  NR==1{ sub(/\.csv$/, "", FILENAME) } # remove .csv suffix from FILENAME
  NR>1{ $1=FILENAME }                  # replace the first field with filename
  1                                    # print record
' *_mRNA_brain.txt_results_doAsso4_cov_default.lrt0.gz_results.csv | column -t > all_brain_signif.csv

# also annotate the list of mRNAs that are associated to eQTLs so that we can look at cis/trans (need to get coordinates for a given mRNA from gtf file)
grep -f brain_mRNA_coordinates.txt *.gtf > brain_mRNA_coordinates_annot.txt
```

Then use R to merge datasets and see how module~trait correlations and module~eQTL correlations match up. Identifies possible QTLs that impact traits via expression.
