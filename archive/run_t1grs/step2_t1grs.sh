#!/bin/bash

# Step 2 is based on GitHub: https://github.com/Gaulton-Lab/t1grs
# Need to have run the code in step 1 to get the container for analysis.

# Set inputs
RKJCOLLAB="/Users/slacksa/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver"
t1grs_dir="/Users/slacksa/repos/explore_matprot/genetics/risk_scores/run_t1grs"  # point to single copy of repo
out_dir="${RKJCOLLAB}/Maternal_Protection/data/teddy/genetics/risk_scores/t1grs"

# Run the Docker container
cd ${t1grs_dir}/t1grs  # so T1GRS run script can be found
apptainer run --bind ${out_dir}:/out_dir \
  ${t1grs_dir}/t1grs/t1grs_latest.sif \
  --vcf_path /out_dir/step1_t1grs_id_fix.vcf \
  --r3_variants_path /data/T1GRS_allele_order.txt \
  --r2_snps_path /data/ALL5_199_TOPMED_SUSIE_HLA_T1D_signals_updateID_r3.vcf.alleles \
  --allele_order_path /data/ALL5_199_TOPMED_SUSIE_HLA_T1D_signals_updateID.vcf.alleles \
  --xgb_all_model_path /data/ALL_NoPCs_Final.ubj \
  --xgb_nohla_model_path /data/NOHLA_NoPCs_Final.ubj \
  --xgb_hlaonly_model_path /data/HLA_ONLY_NoPCs_Final.ubj \
  --all_columns_path /data/ALL_columns.txt \
  --hla_columns_path /data/HLA_columns.txt \
  --nonhla_columns_path /data/nonHLA_columns.txt \
  --total_percentiles_path /data/T1GRS_total_percentile_risk.txt \
  --mhc_percentiles_path /data/T1GRS_MHC_percentile_risk.txt \
  --nonmhc_percentiles_path /data/T1GRS_nonMHC_percentile_risk.txt \
  --output_path /out_dir/T1GRS_probabilities_r3.csv \
  2>&1 | tee ${out_dir}/T1GRS_probabilities_r3.log
