#!/bin/bash

usage() {
  echo "
  Required parameters:
    -v|--vcf                      Path to vcf file
    -b|--bam                      Input bam file. Multiple bam files can be used

  Optional parameters:
    -h|--help                     Shows help
    -t|--threads                  Number of threads [default: $THREADS]
    -c|--control                  Control sample name. Multiple sample names can be used
    -bl|--blacklist               Blacklist vcf or bed file. Multiple files can be used
    -q|--qual                     Flag variants with a low QUAL value [default: $QUAL]
    -mq|--mq                      Flag variants with a low MQ value [default: $MQ]
  "
exit
}

POSITIONAL=()

THREADS=8
QUAL=100
MQ=60

while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
    -h|--help)
    usage
    shift
    ;;
    -v|--vcf)
    VCF="$2"
    shift
    shift
    ;;
    -b|--bed)
    BEDDIR="$2"
    shift
    shift
    ;;
    -a|--age)
    AGE="$2"
    shift
    shift
    ;;
    -s|--shape)
    SHAPE="$2"
    shift
    shift
    ;;
    -c|--color)
    COLOR="$2"
    shift
    shift
    ;;
    -f|--fill)
    FILL="$2"
    shift
    shift
    ;;
    -z|--size)
    SIZE="$2"
    shift
    shift
    ;;
    -l|--lme)
    LME="$2"
    shift
    shift
    ;;
    -d|--donor)
    DONOR="$2"
    shift
    shift
    ;;
    -t|--bedtools)
    BEDTOOLS="$2"
    shift
    shift
    ;;
    -g|--autosomal_genome_size)
    AUTOSOMAL="$2"
    shift
    shift
    ;;
    *)
    POSITIONAL+=("$1")
    shift
    ;;
  esac
done
set -- "${POSITIONAL[@]}"
if [ -z $VCF ]; then
  echo "Missing -v|--vcf parameter"
  usage
elif [ -z $BEDDIR ]; then
  echo "Missing -b|--bed parameter"
  usage
elif [ -z $AGE ]; then
  echo "Missing -a|--age parameter"
  usage
fi

if [ ! -f $VCF ]; then
  echo "VCF file \"$VCF\" does not exists."
  exit
fi
