#!/usr/bin/python

import vcf as pyvcf
import argparse
import os
import time
import collections
import sys
import multiprocessing as mp
import queue
import re
import numpy as np
import pysam

# Set arguments
parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-i', '--input', type=str, help='Input indexed vcf.gz file', required=True)
parser.add_argument('-t', '--threads', default=8, type=int, help="Number of threads (default: %(default)s)")
parser.add_argument('-f', '--effect', default="MODERATE|HIGH", type=str, help="Filter on effect (default: %(default)s)")
parser.add_argument('-g', '--genes', default=[], action='append', nargs="*", type=str, help='Input genes list for germline filtering')
parser.add_argument('-q', '--qual', default=60, type=int, help="Flag variants with a low QUAL value (default: %(default)s)")
parser.add_argument('-mq', '--mq', default=30, type=int, help="Flag variants with a low MQ value (default: %(default)s)")
parser.add_argument('-ct', '--clonal_threshold', default=0.3, type=float, help="Sample reported as subclonal if VAF is lower (default: %(default)s)")
parser.add_argument('-at', '--absent_threshold', default=0.0, type=float, help="Sample reported as absent if VAF is lower(default: %(default)s)")
parser.add_argument('-sgq11','--sample_gq_homozygous', default=10, type=int, help="Minimal Genome Quality of a homozygous SNV in a sample (default: %(default)s)")
parser.add_argument('-sgq01','--sample_gq_heterozygous', default=20, type=int, help="Minimal Genome Quality of a heterozygous SNV in a sample (default: %(default)s)")
parser.add_argument('-cgq11','--control_gq_homozygous', default=10, type=int, help="Minimal Genome Quality of a homozygous SNV in a control (default: %(default)s)")
parser.add_argument('-cgq01','--control_gq_heterozygous', default=10, type=int, help="Minimal Genome Quality of a heterozygous SNV in a control (default: %(default)s)")
parser.add_argument('-igq11','--indel_gq_homozygous', default=60, type=int, help="Minimal Genome Quality of a homozygous indel (default: %(default)s)")
parser.add_argument('-igq01','--indel_gq_heterozygous', default=60, type=int, help="Minimal Genome Quality of a heterozygous indel (default: %(default)s)")
parser.add_argument('-igq00','--indel_gq_homref', default=20, type=int, help="Minimal Genome Quality of a homozygous reference indel (default: %(default)s)")
parser.add_argument('-cov', '--coverage', default=10, type=int, help="Flag variants with a low COV value (default: %(default)s)")
parser.add_argument('-m', '--mapq', default=0, type=int, help="Include only reads with a minimal mapq (default: %(default)s)")
parser.add_argument('-p', '--base_phred_quality', default=0, type=int, help="Include only bases with a minimal base phred quality (default: %(default)s)")
parser.add_argument('-indel', '--indel', default=True, action='store_false', help="Exclude indels")
args = parser.parse_args()

args.genes = [x for l in args.genes for x in l]

vcf_reader = pyvcf.Reader(filename=args.input, encoding='utf-8')
vcf_name = os.path.basename(args.input)
vcf_name = vcf_name.replace(".vcf.gz","")
contig_list = []
genes_list = []
bam_sample_names = collections.defaultdict(dict)
genes_table = {}


args.control = []
for control in re.findall(r"--control\s+(\S+)", vcf_reader.metadata['SMuRFCmd'][0]):
    args.control.append(control)
args.bam = []
for bam in re.findall(r"--bam\s+(\S+)", vcf_reader.metadata['SMuRFCmd'][0]):
    args.bam.append(bam)

# Create tmp directory if it does not exists
try:
    os.stat('./SMuRF_tmp')
except:
    os.mkdir('./SMuRF_tmp')

def main():
    global vcf_reader, vaf_df, blacklist, genes_list

    create_genes_list()

    vcf_reader = fix_vcf_header(vcf_reader)
    vcf_reader = add_vcf_header(vcf_reader)

    for contig in vcf_reader.contigs:
        contig_list.append(contig)

    # Create an input queue with the contigs and an empty output queue
    q = mp.Queue()
    q_out = mp.Queue()
    for contig in contig_list:
        q.put(contig)

    # Create number of processes to parse the vcf file
    processes = [mp.Process(target=parse_chr_vcf, args=(q, q_out, vcf_reader)) for x in range(args.threads)]

    for p in processes:
        p.start()
    liveprocs = list(processes)
    while liveprocs:
        time.sleep(5)
        try:
            while 1:
                contig_genes_table = q_out.get(block=False, timeout=1)
                if contig_genes_table:
                    genes_table.update( contig_genes_table )
        except queue.Empty:
            pass
    # Give tasks a chance to put more data in
        time.sleep(10)
        if not q.empty():
            continue
        liveprocs = [p for p in liveprocs if p.is_alive()]

    for p in processes:
        p.join()

    genes_table_file = open("genes_table.txt",'w')
    genes_table_file.write( "\t".join(['gene','germline','somatic','other', 'FP']) + "\n" )
    for gene_name in genes_table:
        outline = [gene_name]
        for mut_type in ['GL','SOM','OTHER', 'FP']:
            outline.append(str(genes_table[gene_name][mut_type]))
        genes_table_file.write( "\t".join(outline) + "\n")
    genes_table_file.close()

def parse_chr_vcf(q, q_out, contig_vcf_reader):
    """
    Function to parse the vcf per contig.
    Write the new record to a vcf file.
    Input: Queue object
    Input: Queue out object
    Input: VCF reader object
    Input: List with the bam names
    """
    while True:
        try:
            # Get contig one by one from the queue
            contig = q.get(block=False,timeout=1)
            contig_vcf_other_drivers_writer = pyvcf.Writer(open('./SMuRF_tmp/{}_drivers_others.vcf'.format(contig),'w', encoding='utf-8'), contig_vcf_reader)
            contig_vcf_FP_drivers_writer = pyvcf.Writer(open('./SMuRF_tmp/{}_drivers_FP.vcf'.format(contig),'w', encoding='utf-8'), contig_vcf_reader)
            contig_vcf_GL_drivers_writer = pyvcf.Writer(open('./SMuRF_tmp/{}_drivers_GL.vcf'.format(contig),'w', encoding='utf-8'), contig_vcf_reader)
            contig_vcf_SOM_drivers_writer = pyvcf.Writer(open('./SMuRF_tmp/{}_drivers_SOM.vcf'.format(contig),'w', encoding='utf-8'), contig_vcf_reader)
            contig_genes_table = {}
            try:
                # Try to parse the specific contig from the vcf
                contig_vcf_reader.fetch(contig)
            except:
                # Skip contig if it is not present in the vcf file
                continue
            for record in contig_vcf_reader.fetch(contig):
                write = True
                for fil in ['BadQual','BadMQ','NoClonalSample']:
                    if fil in record.FILTER:
                        record.FILTER.remove(fil)
                for fil in ['ControlEvidence','ControlClonal','ControlSubclonal']:
                    if fil in record.FILTER:
                        record.FILTER.remove(fil)
                    gl = True
                if not record.FILTER:
                    if record.QUAL < args.qual:
                        record.FILTER.append('BadQual')
                        write = False
                    elif record.INFO['MQ'] < args.mq:
                        record.FILTER.append('BadMQ')
                        write = False
                    elif sample_quality_control( record ):
                        if gl and not ann_check(record):
                            write = False
                            continue
                        if "CLONAL_SAMPLES" not in record.INFO:
                            if not calculate_vaf( record ):
                                write = False
                                continue
                        if write:
                            gl, som = check_gl_som( record )
                            if gl and not ann_check(record):
                                write = False
                                continue
                            if gl and not som:
                                gene_name = get_gene_name( record )
                                if gene_name not in contig_genes_table:
                                    contig_genes_table[gene_name] = {"GL":0, "SOM":0, "OTHER":0, "FP":0 }
                                contig_genes_table[gene_name]['GL'] += 1
                                contig_vcf_GL_drivers_writer.write_record(record)
                            elif not gl and som:
                                gene_name = get_gene_name( record )
                                if gene_name not in contig_genes_table:
                                    contig_genes_table[gene_name] = {"GL":0, "SOM":0, "OTHER":0, "FP":0 }
                                contig_genes_table[gene_name]['SOM'] += 1
                                contig_vcf_SOM_drivers_writer.write_record(record)
                            elif not gl and not som:
                                gene_name = get_gene_name( record )
                                if gene_name not in contig_genes_table:
                                    contig_genes_table[gene_name] = {"GL":0, "SOM":0, "OTHER":0, "FP":0 }
                                contig_genes_table[gene_name]['FP'] += 1
                                contig_vcf_FP_drivers_writer.write_record(record)
                            else:
                                gene_name = get_gene_name( record )
                                if gene_name not in contig_genes_table:
                                    contig_genes_table[gene_name] = {"GL":0, "SOM":0, "OTHER":0, "FP":0 }
                                contig_genes_table[gene_name]['OTHER'] += 1
                                contig_vcf_other_drivers_writer.write_record(record)
            q_out.put( contig_genes_table )

        # Break the loop if the queue is empty
        except queue.Empty:
            break

def get_gene_name( record ):
    if 'ANN' in record.INFO:
        annlist = record.INFO['ANN'][0].split("|")
        gene_name = annlist[3]
        return( gene_name )
    return( False )

def sample_quality_control( record ):
    """
    Function to check if the sample pass the quality control.
    Input: VCF record object
    Return: True of False if all samples or controls failed the quality control
    """
    qc = collections.defaultdict(dict)
    noSampleEvidence = 0
    controlEvidence = False
    indel = False
    check = True

    # Check if variant is an indel
    if (len(record.REF) > 1 or len(record.ALT[0]) > 1):
        indel = True

    for call in (record.samples):
        sample = True
        # update_call_data(call, ['FT'], [None], vcf_reader)

        # Check is the sample is in the control list
        if call.sample in args.control:
            sample = False
        # QC fails if there is no genotype
        if call['GT'] == './.':
            qc[sample][call.sample] = 'NoGenoType'
        # QC fails if the coverage is too low
        elif (call['DP'] == None or call['DP'] < args.coverage):
            qc[sample][call.sample] = 'LowCov'
        elif sample:
            # If sample is homozygous reference
            if call['GT'] == '0/0':
                noSampleEvidence += 1
                if indel and (not call['GQ'] or call['GQ'] < args.indel_gq_homref):
                    qc[sample][call.sample] = 'LowGQ'
                elif not indel and (not call['GQ'] or call['GQ'] < args.sample_gq_homozygous):
                    qc[sample][call.sample] = 'LowGQ'
                else:
                    qc[sample][call.sample] = 'PASS'
            # Check QC for homozygous variant
            elif call['GT'] == '1/1':
                if indel and (not call['GQ'] or call['GQ'] < args.indel_gq_homozygous):
                    qc[sample][call.sample] = 'LowGQ'
                elif not indel and (not call['GQ'] or call['GQ'] < args.sample_gq_homozygous):
                    qc[sample][call.sample] = 'LowGQ'
                else:
                    qc[sample][call.sample] = 'PASS'
            elif call['GT'] == '0/1':
                if indel and (not call['GQ'] or call['GQ'] < args.indel_gq_heterozygous):
                    qc[sample][call.sample] = 'LowGQ'
                elif not indel and (not call['GQ'] or call['GQ'] < args.sample_gq_heterozygous):
                    qc[sample][call.sample] = 'LowGQ'
                else:
                    qc[sample][call.sample] = 'PASS'
        else:
            # If variant is also found in a control
            if call['GT'] == '0/1':
                controlEvidence = True
                if indel and (not call['GQ'] or call['GQ'] < args.indel_gq_heterozygous):
                    qc[sample][call.sample] = 'LowGQ'
                elif not indel and (not call['GQ'] or call['GQ'] < args.control_gq_heterozygous):
                    qc[sample][call.sample] = 'LowGQ'
                else:
                    qc[sample][call.sample] = 'PASS'
            elif call['GT'] == '1/1':
                controlEvidence = True
                if indel and (not call['GQ'] or call['GQ'] < args.indel_gq_homozygous):
                    qc[sample][call.sample] = 'LowGQ'
                elif not indel and (not call['GQ'] or call['GQ'] < args.control_gq_homozygous):
                    qc[sample][call.sample] = 'LowGQ'
                else:
                    qc[sample][call.sample] = 'PASS'
            elif call['GT'] == '0/0':
                if indel and (not call['GQ'] or call['GQ'] < args.indel_gq_homref):
                    qc[sample][call.sample] = 'LowGQ'
                elif not indel and (not call['GQ'] or call['GQ'] < args.control_gq_homozygous):
                    qc[sample][call.sample] = 'LowGQ'
                else:
                    qc[sample][call.sample] = 'PASS'
        if call.sample in qc[sample]:
            update_call_data(call,['FT'],[qc[sample][call.sample]], vcf_reader)
    format_list = list(vcf_reader.formats.keys())
    format_list.remove('GT')
    format_list.insert(0,'GT')

    # Add VAF information to the format field of each sample
    if 'AllControlsFailedQC' in record.FILTER:
        record.FILTER.remove('AllControlsFailedQC')
    if 'AllSamplesFailedQC' in record.FILTER:
        record.FILTER.remove('AllSamplesFailedQC')
    record.FORMAT = ":".join(format_list)

    if len(qc[False].keys()) > 0 and list(qc[False].values()).count('PASS') == 0:
        record.FILTER.append('AllControlsFailedQC')
        check = False
    elif list(qc[True].values()).count('PASS') == 0:
        record.FILTER.append('AllSamplesFailedQC')
        check = False
    else:
        record.INFO['PASS_QC_SAMPLES'] = list(qc[True].values()).count('PASS')
        record.INFO['FAIL_QC_SAMPLES'] = len(qc[True])-list(qc[True].values()).count('PASS')
        record.INFO['PASS_QC_SAMPLE_NAMES'] = list(np.array(list(qc[True].keys()))[list(np.where(np.array(list(qc[True].values())) == 'PASS')[0])])
        record.INFO['FAIL_QC_SAMPLE_NAMES'] = list(np.array(list(qc[True].keys()))[list(np.where(np.array(list(qc[True].values())) != 'PASS')[0])])
        if ( len(qc[False].keys()) > 0 ):
            record.INFO['PASS_QC_CONTROLS'] = list(qc[False].values()).count('PASS')
            record.INFO['FAIL_QC_CONTROLS'] = len(qc[False])- list(qc[False].values()).count('PASS')
            record.INFO['PASS_QC_CONTROL_NAMES'] = list(np.array(list(qc[False].keys()))[list(np.where(np.array(list(qc[False].values())) == 'PASS')[0])])
            record.INFO['FAIL_QC_CONTROL_NAMES'] = list(np.array(list(qc[False].keys()))[list(np.where(np.array(list(qc[False].values())) != 'PASS')[0])])

    return( check )

def create_genes_list():
    for genes_file in args.genes:
        with(open(genes_file,'r')) as f:
            for line in f:
                line = line.rstrip()
                genes_list.append(line)

def ann_check( record ):
    ann_check = False
    if 'ANN' in record.INFO:
        for ann in record.INFO['ANN']:
            annlist = ann.split("|")
            if re.search(args.effect, annlist[2]) is not None:
                if len(args.genes) > 0:
                    if annlist[3] in genes_list:
                        ann_check = True
                else:
                    ann_check = True
    return( ann_check )

def get_sample_name( bamfile ):
    """
    Function to get the sample name from the bam file
    Input: An AlignmentFile object of the bam file
    Return: The sample or False if there is no SM tag in the bam header
    """
    header = bamfile.header
    sample_name = False
    if 'RG' in header:
        if type(header['RG']) is list:
            sample_name = header['RG'][0]['SM']
        else:
            sample_name = header['RG']['SM']

    return( sample_name )

def check_pileupread( pileupread ):
    """
    Function to check a pileup read.
    Returns True if the read needs to be kept and returns False if read can be skipped.
    Input: Pileupread object
    Return: True or False
    """
    check = True
    if pileupread.alignment.is_duplicate:
        check = False
    elif pileupread.is_del:
        check = False
    elif pileupread.is_refskip:
        check = False
    elif not pileupread.query_position:
        check = False
    elif pileupread.alignment.mapq < args.mapq:
        check = False
    elif pileupread.alignment.query_qualities[pileupread.query_position] < args.base_phred_quality:
        check = False

    return( check )

def check_gl_som( record ):
    germline = None
    somatic = None
    if record.INFO['SUBCLONAL_CONTROLS'] == 0 and record.INFO['CLONAL_CONTROLS'] == 0 and record.INFO['SUBCLONAL_SAMPLES'] == 0 and record.INFO['CLONAL_SAMPLES'] == 0:
        somatic = False
        germline = False
    if ( record.INFO['ABSENT_CONTROLS'] > 0 and record.INFO['CLONAL_CONTROLS'] == 0 and record.INFO['SUBCLONAL_CONTROLS'] == 0 ):
        if ( record.INFO['CLONAL_SAMPLES'] > 0 ):
        # if ( record.INFO['CLONAL_SAMPLES'] > 0 or record.INFO['SUBCLONAL_SAMPLES'] > 0 ):
            somatic = True
    if ( record.INFO['ABSENT_CONTROLS'] == 0 and record.INFO['CLONAL_CONTROLS'] > 0 ):
    # if ( record.INFO['ABSENT_CONTROLS'] == 0 and ( record.INFO['CLONAL_CONTROLS'] > 0 or record.INFO['SUBCLONAL_CONTROLS'] > 0 ) ):
        if ( record.INFO['ABSENT_SAMPLES'] == 0 and record.INFO['CLONAL_SAMPLES'] > 0 ) :
        # if ( record.INFO['ABSENT_SAMPLES'] == 0 and ( record.INFO['CLONAL_SAMPLES'] > 0 or record.INFO['SUBCLONAL_SAMPLES'] > 0 ) ) :
            germline = True

    return( germline, somatic)

def update_call_data( call, edit_keys, edit_values, vcf_reader ):
    """
    Function to add or update a field to the format field.
    This will be automatically update in the call object
    Input: A call object
    Input: A list with format fields
    Input: A list with format values
    """
    f_keys = list(vcf_reader.formats.keys())
    d = dict(call.data._asdict())
    f_vals = []
    for key in f_keys:
        if key in edit_keys:
            f_vals.append(edit_values[edit_keys.index(key)] )
        elif key in d:
            f_vals.append(d[key] )
        else:
            f_vals.append(None)
    handy_dict = dict(zip(f_keys, f_vals))
    f_keys.remove('GT')
    f_keys.insert(0,'GT')
    call.data = collections.namedtuple('CallData',f_keys)(**handy_dict)


def calculate_vaf( record ):
    """
    Function to calculate the VAF in the bam files
    Input: VCF record object
    Return: A tuple with VAF of each sample
    """
    record_vaf = {}
    vaf_info = collections.defaultdict(lambda: collections.defaultdict(list))
    qc = collections.defaultdict(dict)

    for call in (record.samples):
        # Add empty VAF and CAD tag to the record
        update_call_data(call, ['VAF','CAD'], [None, None], vcf_reader)
    for bam in args.bam:
        F=pysam.AlignmentFile(bam,'rb')
        if bam not in bam_sample_names:
            sample_name = get_sample_name(F)
            bam_sample_names[bam] = sample_name
        else:
            sample_name = bam_sample_names[bam]
        dv = 0
        dr = 0
        vaf = 0.0
        # Loop through each reads that is overlapping the position of the variant
        for pileupcolumn in F.pileup(record.CHROM, int(record.POS)-1, int(record.POS), truncate=True, stepper='nofilter',min_base_quality=args.base_phred_quality):
            for pileupread in pileupcolumn.pileups:
                # QC the read
                if ( check_pileupread( pileupread) ):
                    alt = record.ALT[0]
                    # If variant is a SNV
                    if (len(record.REF) == 1 and len(alt) == 1):
                        # Read has the reference
                        if pileupread.alignment.query_sequence[pileupread.query_position] == record.REF:
                            dr+=1
                        # Read has the variant
                        elif pileupread.alignment.query_sequence[pileupread.query_position] == alt:
                            dv+=1
                    # If variant is deletion
                    elif (len(record.REF) > 1 and len(alt) == 1):
                        # Read has the deletion
                        if ( pileupread.indel*-1 == len(record.REF)-1 ):
                            dv+=1
                        # Read has no deletion
                        elif pileupread.indel == 0:
                            dr+=1
                    # If variant is an insertion
                    elif ( len(record.REF) == 1 and len(alt) > 1 ):
                        # Read has the insertion
                        if ( pileupread.indel == len(alt)-1 ):
                            dv+=1
                        # Read has no insertion
                        elif pileupread.indel == 0:
                            dr+=1
                    # If variant is an INDEL
                    else:
                        # Read has the INDEL
                        if ( pileupread.indel == (len(alt)-len(record.REF)) ):
                            dv+=1
                        # Read has no INDEL
                        elif pileupread.indel == 0:
                            dr+=1
        # Calculate the VAF
        try:
            vaf = float("{0:.2f}".format(dv/float(dv+dr)))
        except ZeroDivisionError:
            continue
        # Loop through each sample in the vcf file
        for call in (record.samples):
            sample = True
            if call.sample in args.control:
                sample = False
            # Check if the sample name in the vcf file is the same as a sample name in the bam file
            if call.sample == sample_name:
                # Add the VAF and sample name to the output tuple
                if vaf > args.absent_threshold and call.sample not in args.control:
                    record_vaf[call.sample] = vaf
                # Update the format field for this sample
                ft = call['FT']
                if float(dv+dr) < args.coverage:
            	     ft = 'LowCov'
            	     qc[sample][call.sample] = 'LowCov'
                else:
                    qc[sample][call.sample] = call['FT']
                update_call_data(call, ['VAF','CAD','FT'], [vaf, [dr, dv], ft], vcf_reader)

                # Set absent, subclonal or clonal based on the VAF and threshold
                if vaf <= args.absent_threshold:
                    vaf_info[sample]['ABSENT'].append(call.sample)
                elif vaf < args.clonal_threshold:
                    vaf_info[sample]['SUBCLONAL'].append(call.sample)
                else:
                    vaf_info[sample]['CLONAL'].append(call.sample)

    format_list = list(vcf_reader.formats.keys())
    format_list.remove('GT')
    format_list.insert(0,'GT')

    # Add VAF information to the format field of each sample
    record.FORMAT = ":".join(format_list)
    # Add QC information to the INFO field
    check = True

    if len(qc[False].keys()) > 0 and list(qc[False].values()).count('PASS') == 0:
        record.FILTER.append('AllControlsFailedQC')
        check = False
    elif list(qc[True].values()).count('PASS') == 0:
        record.FILTER.append('AllSamplesFailedQC')
        check = False
    else:
        record.INFO['PASS_QC_SAMPLES'] = list(qc[True].values()).count('PASS')
        record.INFO['FAIL_QC_SAMPLES'] = len(qc[True])-list(qc[True].values()).count('PASS')
        record.INFO['PASS_QC_SAMPLE_NAMES'] = list(np.array(list(qc[True].keys()))[list(np.where(np.array(list(qc[True].values())) == 'PASS')[0])])
        record.INFO['FAIL_QC_SAMPLE_NAMES'] = list(np.array(list(qc[True].keys()))[list(np.where(np.array(list(qc[True].values())) != 'PASS')[0])])
        if ( len(qc[False].keys()) > 0 ):
            record.INFO['PASS_QC_CONTROLS'] = list(qc[False].values()).count('PASS')
            record.INFO['FAIL_QC_CONTROLS'] = len(qc[False])- list(qc[False].values()).count('PASS')
            record.INFO['PASS_QC_CONTROL_NAMES'] = list(np.array(list(qc[False].keys()))[list(np.where(np.array(list(qc[False].values())) == 'PASS')[0])])
            record.INFO['FAIL_QC_CONTROL_NAMES'] = list(np.array(list(qc[False].keys()))[list(np.where(np.array(list(qc[False].values())) != 'PASS')[0])])

        # Add clonal information to the INFO field
        record.INFO['ABSENT_SAMPLES'] = len(vaf_info[True]['ABSENT'])
        record.INFO['SUBCLONAL_SAMPLES'] = len(vaf_info[True]['SUBCLONAL'])
        record.INFO['CLONAL_SAMPLES'] = len(vaf_info[True]['CLONAL'])
        record.INFO['ABSENT_CONTROLS'] = len(vaf_info[False]['ABSENT'])
        record.INFO['SUBCLONAL_CONTROLS'] = len(vaf_info[False]['SUBCLONAL'])
        record.INFO['CLONAL_CONTROLS'] = len(vaf_info[False]['CLONAL'])

        record.INFO['ABSENT_SAMPLE_NAMES'] = vaf_info[True]['ABSENT']
        record.INFO['SUBCLONAL_SAMPLE_NAMES'] = vaf_info[True]['SUBCLONAL']
        record.INFO['CLONAL_SAMPLE_NAMES'] = vaf_info[True]['CLONAL']
        record.INFO['ABSENT_CONTROL_NAMES'] = vaf_info[False]['ABSENT']
        record.INFO['SUBCLONAL_CONTROL_NAMES'] = vaf_info[False]['SUBCLONAL']
        record.INFO['CLONAL_CONTROL_NAMES'] = vaf_info[False]['CLONAL']

    return( check )

def get_command_line():
    """
    Function to get the commandline arguments
    Return: A string with the actual command.
    """
    cmdline = [sys.argv[0]]
    for arg in vars(args):
        if type(getattr(args,arg)) is list:
            for a in getattr(args,arg):
                    cmdline.append('--{} {}'.format(arg,str(a)))
        else:
            cmdline.append('--{} {}'.format(arg,str(getattr(args,arg))))
    return( '"{}"'.format(" ".join(cmdline)) )

def fix_vcf_header( vcf_reader ):
    """
    Function to fix fields in the vcf header
    Input: A vcf reader object
    Return: The vcf reader object with fixed headers
    """
    #dbNSFP_clinvar_clnsig has a Integer type but sometimes it is a String, e.g. 2|2
    vcf_reader.infos['dbNSFP_clinvar_clnsig'] = pyvcf.parser._Info("dbNSFP_clinvar_clnsig",1,"String","Field 'clinvar_clnsig' from dbNSFP", None, None)
    #dbNSFP_clinvar_golden_stars has a Integer type but sometimes it is a String, e.g. 0|1
    vcf_reader.infos['dbNSFP_clinvar_golden_stars'] = pyvcf.parser._Info("dbNSFP_clinvar_golden_stars",1,"String","Field 'clinvar_golden_stars' from dbNSFP", None, None)
    vcf_reader.infos['dbNSFP_hg18_chr'] = pyvcf.parser._Info("dbNSFP_hg18_chr",1,"String","Field 'hg18_chr' from dbNSFP", None, None)
    vcf_reader.infos['dbNSFP_hg19_chr'] = pyvcf.parser._Info("dbNSFP_hg19_chr",1,"String","Field 'hg19_chr' from dbNSFP", None, None)
    return( vcf_reader )

def add_vcf_header( vcf_reader ):
    """
    Function to add a new field to the vcf header
    Input: A vcf reader object
    Return: The vcf reader object with new headers added
    """
    # Metadata
    vcf_reader.metadata['DriverCmd'] = [get_command_line()]

    return( vcf_reader )

def merge_tmp_vcfs():
    """
    Function to merge all the tmp contig vcf files
    """
    start = time.time()
    header = False
    # Loop through all chromomsomes
    for contig in contig_list:
        if not header:
            os.system('cat SMuRF_tmp/{}_drivers_others.vcf > {}_drivers_others.vcf'.format(contig, vcf_name))
            os.system('cat SMuRF_tmp/{}_drivers_FP.vcf > {}_drivers_FP.vcf'.format(contig, vcf_name))
            os.system('cat SMuRF_tmp/{}_drivers_GL.vcf > {}_drivers_GL.vcf'.format(contig, vcf_name))
            os.system('cat SMuRF_tmp/{}_drivers_SOM.vcf > {}_drivers_SOM.vcf'.format(contig, vcf_name))
            header = True
        else:
            os.system('grep -v \'^#\' SMuRF_tmp/{}_drivers_others.vcf >> {}_drivers_others.vcf'.format(contig, vcf_name))
            os.system('grep -v \'^#\' SMuRF_tmp/{}_drivers_FP.vcf >> {}_drivers_FP.vcf'.format(contig, vcf_name))
            os.system('grep -v \'^#\' SMuRF_tmp/{}_drivers_GL.vcf >> {}_drivers_GL.vcf'.format(contig, vcf_name))
            os.system('grep -v \'^#\' SMuRF_tmp/{}_drivers_SOM.vcf >> {}_drivers_SOM.vcf'.format(contig, vcf_name))

if __name__ == "__main__":
    #get_command_line()
    main()
    merge_tmp_vcfs()
