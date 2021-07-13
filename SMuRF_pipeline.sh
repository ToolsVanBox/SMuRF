#!/bin/bash

usage() {
  echo "
  Required parameters:
    -i|--input                     Path to vcf file
    -b|--bam                      Input bam file. Multiple bam files can be used. Note, use quotes if you make use of wildcard, eg: '/path/to/*/*bam'

  Optional parameters:
    -h|--help                     Shows help
    -o|--outputdir                Output directory ( default: $OUTDIR )
    -n|--normal                   Normal sample name. Multiple sample names can be used
    -q|--queue                    Queuing system (slurm or sge) (default: $QUEUE)
    -c|--config                   Give the full path to your own ini file (default: $CONFIG)
    -d|--driver                   Run driver genes script (true|false) ( default: $DRIVER )
    -m|--split_mut_type           Run split mutation type script (true|false) (default: $SPLIT_MUT_TYPE)
    -s|--split_single_sample      Run split single sample script (true|false) (default: $SPLIT_SINGLE_SAMPLE)
    -b|--bgzip                    Path to bgzip (default: $BGZIP)
    -t|--tabix                    Patho to tabix (default: $TABIX)

    -mem|--mem                    Memory (default: $MEM)
    -time|--time                  Time (default: $TIME)
  "
exit
}

POSITIONAL=()

SOURCE=${BASH_SOURCE[0]}
SOURCE_DIR=$(dirname $SOURCE)
CONFIG=$SOURCE_DIR/config.ini
VENV=$SOURCE_DIR/venv_3.6/bin/activate
OUTDIR=./
QUEUE=slurm
DRIVER=true
SPLIT_MUT_TYPE=true
SPLIT_SINGLE_SAMPLE=true
BGZIP=/hpc/local/CentOS7/pmc_vanboxtel/bin/samtools-bcftools-htslib-1.0_x64-linux/bin/bgzip
TABIX=/hpc/local/CentOS7/pmc_vanboxtel/bin/samtools-bcftools-htslib-1.0_x64-linux/bin/tabix

MEM=20G
TIME='2:0:0'

SMURF_SH=SMuRF.sh
SMURF_ERR=SMuRF.err
SMURF_LOG=SMuRF.log
PWD=`pwd`

BAM=()
NORMAL=()

while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
    -h|--help)
    usage
    shift
    ;;
    -i|--input)
    INPUT=$(realpath $2)
    shift
    shift
    ;;
    -b|--bam)
    BAM+=$(realpath $2)
    shift
    shift
    ;;
    -n|--normal)
    NORMAL+=("$2")
    shift
    shift
    ;;
    -o|--outputdir)
    OUTDIR=$(realpath $2)
    shift
    shift
    ;;
    -q|--queue)
    QUEUE="$2"
    shift
    shift
    ;;
    -c|--config)
    CONFIG=$(realpath $2)
    shift
    shift
    ;;
    -d|--driver)
    DRIVER="$2"
    shift
    shift
    ;;
    -m|--split_mut_type)
    SPLIT_MUT_TYPE="$2"
    shift
    shift
    ;;
    -s|--split_single_sample)
    SPLIT_SINGLE_SAMPLE="$2"
    shift
    shift
    ;;
    -mem|--mem)
    MEM="$2"
    shift
    shift
    ;;
    -time|--time)
    TIME="$2"
    shift
    shift
    ;;
    -b|--bgzip)
    BGZIP=$(realpath $2)
    shift
    shift
    ;;
    -t|--tabix)
    TABIX=$(realpath $2)
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
if [ -z $INPUT ]; then
  echo "Missing -i|--input parameter"
  usage
elif [ ${#BAM[@]} == 0 ]; then
  echo "Missing -b|--bam parameter"
  usage
fi

# Variables check
if [ ! -f $INPUT ]; then
  echo "VCF file \"$INPUT\" does not exists."
  exit
elif [ ! -f $CONFIG ]; then
  echo "Config file \"$CONFIG\" does not exists."
  exit
elif [ ! -f $BGZIP ]; then
  echo "BGZIP \"$BGZIP\" does not exists."
  exit
elif [ ! -f $TABIX ]; then
  echo "TABIX \"$TABIX\" does not exists."
  exit
fi
if [[ $QUEUE != slurm && $QUEUE != sge ]]; then
  echo "QUEUE parameter must be slur or sge"
  exit
fi
if [[ $DRIVER != true && $DRIVER != false ]]; then
  echo "DRIVER parameter must be true or false"
  exit
elif [[ $SPLIT_MUT_TYPE != true && $SPLIT_MUT_TYPE != false ]]; then
  echo "SPLIT_MUT_TYPE parameter must be true or false"
  exit
elif [[ $SPLIT_SINGLE_SAMPLE != true && $SPLIT_SINGLE_SAMPLE != false ]]; then
  echo "SPLIT_SINGLE_SAMPLE parameter must be true or false"
  exit
fi

for B in ${BAM[@]}; do
  if [ ! -f $B ]; then
    echo "BAM file \"$B\" does not exists."
    exit
  fi
done

if [ ! -d $OUTDIR ]; then
  mkdir -p $OUTDIR
fi


THREADS=`grep "SMuRF" -A 19 ${CONFIG} | grep -m 1 threads | cut -f 3 -d ' '`

cat << EOF > $SMURF_SH
#!/bin/bash
EOF

if [ $QUEUE == 'sge' ]; then
cat << EOF >> $SMURF_SH
#$ -N SMuRF
#$ -cwd
#$ -pe threaded $THREADS
#$ -l h_vmem=$MEM
#$ -l h_rt=$TIME
#$ -e $SMURF_ERR
#$ -o $SMURF_LOG
EOF
elif [ $QUEUE == 'slurm' ]; then
cat << EOF >> $SMURF_SH
#SBATCH --job-name=SMuRF
#SBATCH -c $THREADS
#SBATCH --mem=$MEM
#SBATCH --time=$TIME
#SBATCH -e $SMURF_ERR
#SBATCH -o $SMURF_LOG
EOF
fi

cat << EOF >> $SMURF_SH
cd $(realpath $OUTDIR)

. $VENV
EOF

# Create bgzip vcf if needed
if [[ $INPUT == *.vcf ]]; then
  INPUT_BGZIP=$INPUT.gz
  INPUT_BGZIP=$(basename $INPUT_BGZIP)
cat << EOF >> $SMURF_SH
  $BGZIP -c $(realpath $INPUT) > $(realpath $INPUT_BGZIP)
  $TABIX $(realpath $INPUT_BGZIP)
EOF
else
  INPUT_BGZIP=$(basename $INPUT)
  if [ ! -f $INPUT_BGZIP ]; then
cat << EOF >> $SMURF_SH
      ln -s $(realpath $INPUT) $INPUT_BGZIP
EOF
  fi
  if [ -f $INPUT.tbi ]; then
cat << EOF >> $SMURF_SH
      ln -s $(realpath $INPUT.tbi) $INPUT_BGZIP.tbi
EOF
  else
cat << EOF >> $SMURF_SH
     $TABIX $INPUT_BGZIP
EOF
  fi
fi

SMURF_VCF=${INPUT_BGZIP%.vcf.gz}_SMuRF.vcf
SMURF_VCF_BGZIP=$SMURF_VCF.gz
SMURF_VCF_FILTERED=${SMURF_VCF%.vcf}_filtered.vcf

# # Create bgzip vcf if needed
# if [[ $INPUT == *.vcf ]]; then
#   $BGZIP -c $INPUT > $INPUT_BGZIP
#   $TABIX $INPUT_BGZIP
# else
#   ln -s $INPUT $INPUT_BGZIP
# fi
#
cat << EOF >> $SMURF_SH
if [ ! -f $(realpath $SMURF_VCF) ]; then
  python $SOURCE_DIR/SMuRF.py \\
  -i $(realpath $INPUT_BGZIP) \\
EOF

for B in ${BAM[@]}; do
cat << EOF >> $SMURF_SH
  -b $(realpath $B) \\
EOF
done

for N in ${NORMAL[@]}; do
cat << EOF >> $SMURF_SH
  -n $N \\
EOF
done

cat << EOF >> $SMURF_SH
  -c $(realpath $CONFIG)
fi
EOF

if [ $SPLIT_MUT_TYPE == true ]; then
cat << EOF >> $SMURF_SH

bash $SOURCE_DIR/scripts/split_mutation_type.sh $(realpath $SMURF_VCF_FILTERED)
EOF
fi

if [ $SPLIT_SINGLE_SAMPLE == true ]; then
cat << EOF >> $SMURF_SH

bash $SOURCE_DIR/scripts/split_in_single_sample_vcfs.sh $(realpath $SMURF_VCF_FILTERED)
EOF
fi

if [ $DRIVER == true ]; then
cat << EOF >> $SMURF_SH

$BGZIP -c $(realpath $SMURF_VCF) > $(realpath $SMURF_VCF_BGZIP)
$TABIX $(realpath $SMURF_VCF_BGZIP)
python $SOURCE_DIR/scripts/driver_mutations_filter.py -i $(realpath $SMURF_VCF_BGZIP) -c $(realpath $CONFIG)
EOF
fi

#qsub $SMURF_SH
