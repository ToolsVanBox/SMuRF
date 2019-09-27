#!/bin/bash

#$ -N TOOL
#$ -cwd
#$ -pe threaded 8
#$ -l h_rt=1:0:0

VCF=/path/to/file.vcf.gz

. /hpc/pmc_vanboxtel/tools/TOOL/venv_3.6/bin/activate

python /home/pmc_research/mroosmalen/pmc_vanboxtel/tools/TOOL/tool.py \
-i $VCF \
-b /path/to/*.bam \
-t 8 \
-bl /hpc/pmc_vanboxtel/data/Mutation_blacklists/MSC_healthyBM_raw_variants.vcf \
-c CONTROL_NAME_1 \
-c CONTROL_NAME_2
