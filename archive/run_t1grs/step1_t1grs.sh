#!/bin/bash

# SDS 20250714

# Step 1 is based on GitHub: https://github.com/Gaulton-Lab/extract-TOPMed-Michigan-HLA/tree/main
# Before running, run the following code once:

# Use Apptainer to get Docker image
    # Note that repo documentation is off - t1d-grs-analysis:latest  runs both steps combined
# cd /Users/slacksa/repos/explore_matprot/genetics/risk_scores/run_t1grs
# apptainer pull docker://kgaultonlab/t1grs:latest  # this runs both steps combined

# Set inputs
tm_plink="${RKJCOLLAB}/Maternal_Protection/data/teddy/genetics/imputation/tm_r3_imp/hwe_in_all/imputed_clean_maf0_rsq0.3/chr_all_concat"
hla_raw="${RKJCOLLAB}/Maternal_Protection/data/teddy/genetics/imputation/mich_hla_imp_v1/imputed"
tmp_dir="/Users/slacksa/temp_data/maternal_protection/t1grs"
t1grs_dir="/Users/slacksa/repos/explore_matprot/genetics/risk_scores/run_t1grs"  # point to single copy of repo
out_dir="${RKJCOLLAB}/Maternal_Protection/data/teddy/genetics/risk_scores/t1grs"

# Designed to run on raw TOPMed imputation output, so need to remake chr#.dose.vcf.gz files
# from the cleaned files to input
for chr in {1..22}; do
    plink2 --pfile "$tm_plink" \
        --chr $chr --recode vcf bgz \
        --out "${tmp_dir}/chr${chr}.dose"
done

# Modify given variant file to remove prefix of "chr" on variant IDs to match
# standard TOPMed output without "chr"
sed 's/^chr//' "${t1grs_dir}/t1grs/data/ALL5_199_TOPMED_SUSIE_HLA_T1D_signals_updateID_r3.vcf.alleles" > \
    "${t1grs_dir}/t1grs/data/ALL5_199_TOPMED_SUSIE_HLA_T1D_signals_updateID_r3_no_chr.vcf.alleles"

bash ${t1grs_dir}/extract-TOPMed-Michigan-HLA/extract_TOPMed_R3_Michigan_HLA.sh \
    /usr/local/bin/plink2 \
    "${tmp_dir}" \
    "${hla_raw}" \
    "${t1grs_dir}/t1grs/data/ALL5_199_TOPMED_SUSIE_HLA_T1D_signals_updateID_r3_no_chr.vcf.alleles" \
    "${out_dir}/step1_t1grs.vcf"

# Move outputs from code dir to out dir
mv chr* "$out_dir"
mv Mich* "$out_dir"

# Update TOPMed variants from chr:pos:ref:alt to RSIDs
Rscript ${t1grs_dir}/t1grs/add_t1grs_rsids.R \
  "${out_dir}/step1_t1grs.vcf" "${t1grs_dir}/t1grs/data/T1GRS_allele_order.txt"

# Clean up
rm -f ${tmp_dir}/tmp*
rm -f ${out_dir}/tmp*
