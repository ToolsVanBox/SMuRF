#!/bin/bash

VCF=$1

NCOL=$(grep "^#CHROM" $VCF | awk '{print NF}')

for NC in $(seq 10 $NCOL); do
  SAMPLE_NAME=$(grep "^#CHROM" $VCF | cut -f $NC)
  SAMPLE_OUT_VCF=${VCF%.SMuRF.filtered.vcf}_${SAMPLE_NAME}.SMuRF.filtered.vcf
  grep -P "^##" $VCF > $SAMPLE_OUT_VCF
  grep -P "^#CHROM" $VCF | cut -f -9,$NC >> $SAMPLE_OUT_VCF
  grep -P ";CLONAL_SAMPLE_NAMES([=\w-,]+|=)$SAMPLE_NAME[,\s;]+" $VCF | grep -P "PASS_QC_SAMPLE_NAMES([=\w-,]+|=)$SAMPLE_NAME[,\s;]+" | cut -f -9,$NC >> $SAMPLE_OUT_VCF
done
