#!/bin/bash

#$ -N SMuRF
#$ -cwd
#$ -pe threaded 8
#$ -l h_rt=1:0:0

BGZIP=/hpc/local/CentOS7/pmc_vanboxtel/bin/samtools-bcftools-htslib-1.0_x64-linux/bin/bgzip
TABIX=/hpc/local/CentOS7/pmc_vanboxtel/bin/samtools-bcftools-htslib-1.0_x64-linux/bin/tabix

VCF=/path/to/file.vcf{.gz}

if [[ $VCF == *.vcf ]]; then
  VCF_BGZIP=$VCF.gz
  $BGZIP $VCF
  $TABIX $VCF_BGZIP
else
  VCF_BGZIP=$VCF
fi

SMURF_VCF=${VCF_BGZIP%.vcf.gz}_SMuRF.vcf
SMURF_VCF_BGZIP=$SMURF_VCF.gz

. /hpc/pmc_vanboxtel/tools/SMuRF/venv_3.6/bin/activate

python /hpc/pmc_vanboxtel/tools/SMuRF/SMuRF.py \
-i $VCF_BGZIP \
-b /path/to/*.bam \
-t 8 \
-bl /hpc/pmc_vanboxtel/data/Mutation_blacklists/MSC_healthyBM_raw_variants_hg38.bed \
-c CONTROL_NAME_1 \
-c CONTROL_NAME_2

$BGZIP $SMURF_VCF
$TABIX $SMURF_VCF_BGZIP

python /hpc/pmc_vanboxtel/tools/SMuRF/scripts/driver_mutations_filter.py \
-i $SMURF_VCF_BGZIP \
-g /hpc/pmc_vanboxtel/data/hotspot_genes/cosmic_cancer_gene_census_09052019.txt
