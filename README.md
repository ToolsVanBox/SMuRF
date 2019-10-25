# TOOL
Tool without a name

## Name ideas
- Mutations
- Variants
- SNV
- Indel

- Annotation
- Flag
- Filter
- Genotyping
- VAF

- Checker
- Tool
- Pipeline

- SMuRF

# How to run

## Input
The vcf file must be zipped by bgzip and indexed by tabix
```
> /hpc/local/CentOS7/pmc_vanboxtel/bin/samtools-bcftools-htslib-1.0_x64-linux/bin/bgzip /path/to/<file.vcf>
```
```
/hpc/local/CentOS7/pmc_vanboxtel/bin/samtools-bcftools-htslib-1.0_x64-linux/bin/tabix /path/to/<file.vcf.gz>
```

## Virtual environment

### Create new venv
If you need to create a new virtual environment
```
> virtualenv venv_3.6 -p /hpc/local/CentOS7/common/lang/python/3.6.1/bin/python
```

### Load venv
Before you run the tool, you need load the virtual environment
```
> . /hpc/pmc_vanboxtel/tools/TOOL/venv_3.6/bin/activate
```

### Install python modules
If you created a new virtual environment, install the required modules
```
> pip install -r requirements.txt
```

## Run tool
Run the tool
```
> python /hpc/pmc_vanboxtel/tools/TOOL/tool.py -i /path/to/<file.vcf.gz> -b /path/to/*.bam -c <NAME_OF_CONTROL> -bl /hpc/pmc_vanboxtel/data/Mutation_blacklists/MSC_healthyBM_raw_variants_hg38.bed

```
