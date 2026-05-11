#!/bin/bash

# Use KING to get PCs and ancestry proportions. Based off of PMID 38453145,
# and GitHub repo https://github.com/chenlab-uva/AncestryInference_KING.

# Get args
study="$1"
study_data="$2"
ref_data="$3"
out_prefix="$4"
cpus="$SLURM_NTASKS"

# Change PLINK1.9 .fam column 6 to be 2
awk '{ $6=2; print }' "${study_data}.fam" > "${study_data}_mod.fam"
cp "${study_data}.bed" "${study_data}_mod.bed"
cp "${study_data}.bim" "${study_data}_mod.bim"

# Get PCs from KING PCA projection
# king -b "${ref_data},${study_data}_mod" --pca --projection \
#    --prefix "${out_prefix}" --cpus "$cpus"

# Run Rscript for ancestry inference
# Rscript Ancestry_Inference.R \
	# "${out_prefix}_pc.txt" "${out_prefix}_popref.txt" "${out_prefix}_script"

# Get ancestry projection directly from KING
king -b "${ref_data},${study_data}_mod" --pca --projection --pngplot \
	--cpus "$cpus"
