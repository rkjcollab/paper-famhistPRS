#!/bin/bash

if [ "$#" -eq 0 ]
then
   echo "Usage: ${0##*/} <input_vcf> <crossmap_chain>"
   echo "       <ref_fasta> <output_dir>"
   echo "For given chromosome, script uses CrossMap.py for"
   echo "liftover according to the given chain and reference"
   echo "files. Outputs VCF and PLINK2 files."
   exit
fi

input_vcf=$1
crossmap_chain=$2
ref_fasta=$3
output_dir=$4

# Get chr# from input VCF file
chr=$(basename $input_vcf | grep -o '[0-9]*')
# Get to build from chain file
to_build=$(basename $crossmap_chain | grep -o '[0-9]*' | tail -n 1)

# Get basename of VCF file
vcf_name_tmp=$(basename ${input_vcf%.*})
vcf_name=$(basename ${vcf_name_tmp%.*})  # removes both .vcf.gz if needed

# Lift over from hg38 to hg19
CrossMap.py vcf \
   $crossmap_chain \
   $input_vcf \
   $ref_fasta \
   "${output_dir}/${vcf_name}_hg${to_build}.vcf.gz"

# Convert to PLINK
plink2 --vcf "${output_dir}/${vcf_name}_hg${to_build}.vcf.gz" \
   --make-pgen --keep-allele-order --sort-vars \
   --id-delim \
   --chr "$chr" \
   --out "${output_dir}/${vcf_name}_${to_build}"

# Cleanup
rm ${output_dir}/tmp_*
