#!/bin/bash

VCF=$1
SNV_OUT=${VCF%.vcf}_SNV.vcf
INDEL_OUT=${VCF%.vcf}_INDEL.vcf

grep -P "^#" $VCF > $SNV_OUT
grep -v "^#" $VCF | awk '{ if ( length($4) == 1 && length($5) == 1 ) { print $0 } }' >> $SNV_OUT

grep -P "^#" $VCF > $INDEL_OUT
grep -v "^#" $VCF | awk '{ if ( length($4) > 1 || length($5) > 1 ) { print $0 } }' >> $INDEL_OUT
