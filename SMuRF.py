#!/usr/bin/python

# Import modules
import vcf as pyvcf
import pysam
import argparse
import multiprocessing as mp
import queue
import time
import sys
import collections
import subprocess
import os
import glob
import pandas as pd
#import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import numpy as np
import re
from sklearn.mixture import GaussianMixture
from matplotlib.backends.backend_pdf import PdfPages
import warnings
from configparser import SafeConfigParser

# Get version from git
#__version__ = subprocess.check_output(["git", "describe"]).strip().decode('UTF-8')
__version__ = 'v3.0.2'

# Set arguments
parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-i', '--input', type=str, help='Input indexed vcf.gz file', required=True)
parser.add_argument('-b', '--bam', action='append', nargs="*", type=str, help='Input bam file', required=True)
parser.add_argument('-n', '--normal', action='append', type=str, help='Normal sample name')
parser.add_argument('-t', '--threads', type=int, help='Number of threads',default=8, required=False)
parser.add_argument('-c','--config', default=os.path.dirname(os.path.abspath(__file__))+"/config.ini",type=str,help='Give the full path to your own ini file (default: %(default)s)')
parser.add_argument('-v', '--version', action='version', version=__version__)
args = parser.parse_args()

# Flatten input list of bam files
args.bam = [x for l in args.bam for x in l]
# Set default control None if no control is given at the command line
if not args.normal:
    args.normal = [ None ]

# cfg = configparser.ConfigParser()
cfg = SafeConfigParser(os.environ)
if not os.path.exists(args.config):
    sys.stderr.write("Config file "+args.config+" does not exists\n")
    sys.exit()

if not args.config == os.path.dirname(os.path.abspath(__file__))+"/config.ini":
    cfg.read([os.path.dirname(os.path.abspath(__file__))+"/config.ini", args.config])
else:
    cfg.read(args.config)

with open('SMuRF_config.ini', 'w') as configfile:
    cfg.write(configfile)

# Read the vcf, fix and add fields to the header
vcf_reader = pyvcf.Reader(filename=args.input, encoding='utf-8')
vcf_name = os.path.basename(args.input)
vcf_name = vcf_name.replace(".vcf.gz","")

# Create tmp directory if it does not exists
try:
    os.stat('./SMuRF_tmp')
except:
    os.mkdir('./SMuRF_tmp')

# Define global variables
vaf_dict = collections.defaultdict(list)
responsibilities_dict = collections.defaultdict(dict)
contig_dict = {}
bam_sample_names = collections.defaultdict(dict)

def main():
    global vcf_reader, vaf_df, blacklist
    vcf_reader = fix_vcf_header(vcf_reader)
    vcf_reader = add_vcf_header(vcf_reader)

    blacklist = create_blacklist()

    for contig in vcf_reader.contigs:
        contig_dict[contig] = vcf_reader.contigs[contig][1]

    contig_list = dict(sorted(contig_dict.items(), key=lambda contig_dict: contig_dict[1], reverse=True)).keys()

    # Create an input queue with the contigs and an empty output queue
    q = mp.Queue()
    q_out = mp.Queue()
    for contig in contig_list:
        q.put(contig)

    # Create number of processes to parse the vcf file
    processes = [mp.Process(target=parse_chr_vcf, args=(q, q_out, vcf_reader, args.bam)) for x in range(int(args.threads))]

    for p in processes:
        p.start()
    liveprocs = list(processes)
    while liveprocs:
        time.sleep(5)
        try:
            while 1:
                for s,v in q_out.get(block=False, timeout=1).items():
                    vaf_dict[s].extend(v)
        except queue.Empty:
            pass
    # Give tasks a chance to put more data in
        time.sleep(10)
        if not q.empty():
            continue
        liveprocs = [p for p in liveprocs if p.is_alive()]

    for p in processes:
        p.join()

def parse_chr_vcf(q, q_out, contig_vcf_reader, bams):
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
            contig_vaf = collections.defaultdict(list)
            contig_vcf_flag_writer = pyvcf.Writer(open('./SMuRF_tmp/{}.SMuRF.vcf'.format(contig),'w', encoding='utf-8'), contig_vcf_reader)
            try:
                # Try to parse the specific contig from the vcf
                contig_vcf_reader.fetch(contig)
            except:
                # Skip contig if it is not present in the vcf file
                continue
            # print( "blacklist", blacklist )
            for record in contig_vcf_reader.fetch(contig):
                known = False
                if type(record.FILTER) != list:
                    record.FILTER = []
                if not record.FILTER:
                    chr = record.CHROM
                    chr = chr.lower()
                    chr = re.sub("chr|chrom", "", chr)
             
#                    if (len(record.ALT[0]) > 1):
#                        check_flanking_indels(record, contig_vcf_reader)
                    if 'known_variant_flag_ids' in cfg['SMuRF'] and 'known_variant_flag_values' in cfg['SMuRF']:
                        flag_ids = cfg['SMuRF']['known_variant_flag_ids'].split(" ")
                        flag_values = cfg['SMuRF']['known_variant_flag_values'].split(" ")
                        if 'CSQ' in record.INFO:
                            if check_VEP_CSQ(contig_vcf_reader, flag_ids, flag_values, record):
 #                               record.FILTER.append("KnownVariant")
                                known = True
                                                    
                        for i in range(0, len(flag_ids)):
                            flag_id = flag_ids[i]
                            flag_value = float(flag_values[i])
                            if flag_id in record.INFO and record.INFO[flag_id][0] > flag_value:
#                                record.FILTER.append("KnownVariant")
                                known = True
                                break
 #                       if next:
#                            continue
                    if known:
                        record.FILTER.append("KnownVariant")
                    elif "MQ" not in record.INFO:
                        record.FILTER.append("NoMQtag")
                    elif record.ID and "COSM" not in record.ID and cfg['SMuRF']['known_variant_flag_ids'] == '':
                        record.FILTER.append('KnownVariant')
                    elif record.QUAL < int(cfg['SMuRF']['qual']):
                        record.FILTER.append('BadQual')
                    elif len(record.ALT) > 1:
                        record.FILTER.append('MultiAllelic')
                    elif record.INFO['MQ'] < int(cfg['SMuRF']['mq']):
                        record.FILTER.append('BadMQ')
                    elif chr in blacklist and record.POS in blacklist[chr]:
                        record.FILTER.append("BlackList")
                    elif (len(record.ALT[0]) > 1 or len(record.REF) > 1) and not cfg['SMuRF']['indel']:
                        record.FILTER.append("Indel")
                    elif sample_quality_control( record ):
                        for s,v in calculate_vaf( record ).items():
                            if ( not record.FILTER or record.FILTER == ['NoClonalSample'] ):
                                contig_vaf[s].append(v)
                contig_vcf_flag_writer.write_record(record)
            q_out.put( contig_vaf )

        # Break the loop if the queue is empty
        except queue.Empty:
            break

def check_VEP_CSQ( tmp_reader, f_ids, f_values, record):
    description = tmp_reader.infos['CSQ'][3].split("|")
    f_idx = [None] * len(f_ids)
    for i in range(0, len(f_ids)):
        if f_ids[i] in description:
            f_idx[i] = description.index( f_ids[i] )
    
    csq = record.INFO['CSQ'][0].split("|")
    for i in range(0, len(f_values)):
        if f_idx[i]:
            if csq[f_idx[i]] and float(csq[f_idx[i]]) > float(f_values[i]):
                return( True )

    return( False )

def check_flanking_indels( record, contig_vcf_reader):
    for record2 in contig_vcf_reader.fetch(record.CHROM, record.POS-int(cfg['SMuRF']['indel_flank']), record.POS+int(cfg['SMuRF']['indel_flank'])):
        if (len(record2.ALT[0]) > 1):
            if (record2.CHROM == record.CHROM and record2.POS == record.POS):
                continue
            for call in (record2.samples):
                sample = True
                if call.sample in args.normal:
                    sample = False
                if not sample and (call['GT'] == '0/1' or call['GT'] == '1/1' or call['GT'] == '0|1' or call['GT'] == '1|1'):
                    record.FILTER.append("FlankingControlEvidence")
                    return( 1 )
    return( 1 )


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
    vcf_reader.metadata['SMuRFCmd'] = [get_command_line()]

    # Formats
    vcf_reader.formats['VAF'] = pyvcf.parser._Format('VAF',None,'Float','Variant Allele Frequency calculated from the BAM file')
    vcf_reader.formats['CAD'] = pyvcf.parser._Format('CAD',None,'Integer','Calculated Allelic Depth, used for VAF calculation')
    vcf_reader.formats['FT'] = pyvcf.parser._Format('FT',None,'String','Sample filter')

    # Filters
    vcf_reader.filters['KnownVariant'] = pyvcf.parser._Filter('KnownVariant','Variant has already an ID, excluding COSMIC_IDs')
    vcf_reader.filters['BadMQ'] = pyvcf.parser._Filter('BadMQ', 'Variant with MQ <'+str(cfg['SMuRF']['mq']))
    vcf_reader.filters['BadQual'] = pyvcf.parser._Filter('BadQual','Variant with a QUAL <'+str(cfg['SMuRF']['qual']))
    vcf_reader.filters['MultiAllelic'] = pyvcf.parser._Filter('MultiAllelic', 'Variant has multiple alternative alleles')
    vcf_reader.filters['BlackList'] = pyvcf.parser._Filter('BlackList', 'Variant exists in a blacklist')
    vcf_reader.filters['Indel'] = pyvcf.parser._Filter('Indel','Variant is an indel')
    vcf_reader.filters['ControlEvidence'] = pyvcf.parser._Filter('ControlEvidence','Variant is also found in a control based on the GT')
    vcf_reader.filters['NoSampleEvidence'] = pyvcf.parser._Filter('NoSampleEvidence','Variant is not found in any of the samples based on the GT')
    vcf_reader.filters['AllSamplesFailedQC'] = pyvcf.parser._Filter('AllSamplesFailedQC', 'All samples failed the quality control')
    vcf_reader.filters['AllControlsFailedQC'] = pyvcf.parser._Filter('AllControlsFailedQC', 'All controls failed the quality control')
    vcf_reader.filters['ControlSubclonal'] = pyvcf.parser._Filter('ControlSubclonal', 'Variant is found as subclonal in a control based on the recalculated VAF')
    vcf_reader.filters['ControlClonal'] = pyvcf.parser._Filter('ControlClonal', 'Variant is found as clonal in a control based on the recalculated VAF')
    vcf_reader.filters['NoClonalSample'] = pyvcf.parser._Filter('NoClonalSample', 'Variant is not found as clonal in any of the samples based on the recalculated VAF')
    vcf_reader.filters['OverlappingIndelInControl'] = pyvcf.parser._Filter('OverlappingIndelInControl', 'Overlapping indel found in control')
    vcf_reader.filters['NoMQtag'] = pyvcf.parser._Filter('NoMQtag', 'No MQ tag found in info field')
    vcf_reader.filters['FlankingControlEvidence'] = pyvcf.parser._Filter('FlankingControlEvidence', 'Flanking indel found in control')

    # Sample filters
    vcf_reader.filters['LowCov'] = pyvcf.parser._Filter('LowCov', 'Variant has a coverage <'+str(cfg['SMuRF']['coverage'])+' in this sample/control')
    vcf_reader.filters['NoGenoType'] = pyvcf.parser._Filter('NoGenoType', 'Genotype is empty for this sample/control')
    vcf_reader.filters['isRef'] = pyvcf.parser._Filter('isRef', 'Genotype is a reference (i.e. reference 0/0)')
    vcf_reader.filters['isVariant'] = pyvcf.parser._Filter('isVariant', 'Genotype is a variant (i.e. not reference 0/0)')
    vcf_reader.filters['LowGQ'] = pyvcf.parser._Filter('LowGQ', 'Variant has a low genome quality for this sample/control')

    # Infos
    vcf_reader.infos['ABSENT_SAMPLES'] = pyvcf.parser._Info('ABSENT_SAMPLES',1,'Integer','Number of samples without the variant', None, None)
    vcf_reader.infos['SUBCLONAL_SAMPLES'] = pyvcf.parser._Info('SUBCLONAL_SAMPLES',1,'Integer','Number of samples with a subclonal variant', None, None)
    vcf_reader.infos['CLONAL_SAMPLES'] = pyvcf.parser._Info('CLONAL_SAMPLES',1,'Integer','Number of samples with a clonal variant', None, None)
    vcf_reader.infos['ABSENT_CONTROLS'] = pyvcf.parser._Info('ABSENT_CONTROLS',1,'Integer','Number of controls without the variant', None, None)
    vcf_reader.infos['SUBCLONAL_CONTROLS'] = pyvcf.parser._Info('SUBCLONAL_CONTROLS',1,'Integer','Number of controls with a subclonal variant', None, None)
    vcf_reader.infos['CLONAL_CONTROLS'] = pyvcf.parser._Info('CLONAL_CONTROLS',1,'Integer','Number of controls with a clonal variant', None, None)
    vcf_reader.infos['ABSENT_SAMPLE_NAMES'] = pyvcf.parser._Info('ABSENT_SAMPLE_NAMES',None,'String','Samples without the variant', None, None)
    vcf_reader.infos['SUBCLONAL_SAMPLE_NAMES'] = pyvcf.parser._Info('SUBCLONAL_SAMPLE_NAMES',None,'String','Samples with a subclonal variant', None, None)
    vcf_reader.infos['CLONAL_SAMPLE_NAMES'] = pyvcf.parser._Info('CLONAL_SAMPLE_NAMES',None,'String','Samples with a clonal variant', None, None)
    vcf_reader.infos['ABSENT_CONTROL_NAMES'] = pyvcf.parser._Info('ABSENT_CONTROL_NAMES',None,'String','Controls without the variant', None, None)
    vcf_reader.infos['SUBCLONAL_CONTROL_NAMES'] = pyvcf.parser._Info('SUBCLONAL_CONTROL_NAMES',None,'String','Controls with a subclonal variant', None, None)
    vcf_reader.infos['CLONAL_CONTROL_NAMES'] = pyvcf.parser._Info('CLONAL_CONTROL_NAMES',None,'String','Controls with a clonal variant', None, None)
    vcf_reader.infos['PASS_QC_SAMPLES'] = pyvcf.parser._Info('PASS_QC_SAMPLES',1,'Integer','Number of samples which pass all quality control filters', None, None)
    vcf_reader.infos['PASS_QC_CONTROLS'] = pyvcf.parser._Info('PASS_QC_CONTROLS',1,'Integer','Number of controls which pass all quality control filters', None, None)
    vcf_reader.infos['FAIL_QC_SAMPLES'] = pyvcf.parser._Info('FAIL_QC_SAMPLES',1,'Integer','Number of samples which failed one or multiple quality control filters', None, None)
    vcf_reader.infos['FAIL_QC_CONTROLS'] = pyvcf.parser._Info('FAIL_QC_CONTROLS',1,'Integer','Number of controls which failed one or multiple quality control filters', None, None)
    vcf_reader.infos['PASS_QC_SAMPLE_NAMES'] = pyvcf.parser._Info('PASS_QC_SAMPLE_NAMES',None,'String','Samples which pass all quality control filters', None, None)
    vcf_reader.infos['PASS_QC_CONTROL_NAMES'] = pyvcf.parser._Info('PASS_QC_CONTROL_NAMES',None,'String','Controls which pass all quality control filters', None, None)
    vcf_reader.infos['FAIL_QC_SAMPLE_NAMES'] = pyvcf.parser._Info('FAIL_QC_SAMPLE_NAMES',None,'String','Samples which failed one or multiple quality control filters', None, None)
    vcf_reader.infos['FAIL_QC_CONTROL_NAMES'] = pyvcf.parser._Info('FAIL_QC_CONTROL_NAMES',None,'String','Controls which failed one or multiple quality control filters', None, None)

    return( vcf_reader )

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
    elif pileupread.alignment.mapq < int(cfg['SMuRF']['mapq']):
        check = False
    elif pileupread.alignment.query_qualities[pileupread.query_position] < int(cfg['SMuRF']['base_phred_quality']):
        check = False

    return( check )

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

def gmm( X ):
    """
    Function to fit a GMM model with two components
    Input: Array with VAF values
    Return: x values
    Return: List with probability values for all components
    Return: List with probability values of each component
    Return: List of the means of each component
    Return: List of standard deviations of each component
    """
    X = np.array(X)
    X = np.reshape(X, [len(X), 1])

    # ------------------------------------------------------------
    # Learn the best-fit GMM models
    #  Here we'll use GMM in the standard way: the fit() method
    #  uses an Expectation-Maximization approach to find the best
    #  mixture of Gaussians for the data

    # fit models with 2 components
    N = np.arange(int(cfg['SMuRF']['min_components']), int(cfg['SMuRF']['max_components']))

    models = [None for i in range(len(N))]
    for i in range(len(N)):
        models[i] = GaussianMixture(N[i]).fit(X)

    # compute the AIC and the BIC
    AIC = [m.aic(X) for m in models]
    BIC = [m.bic(X) for m in models]

    x = np.linspace(0.01, 1, 100)
    x = x.reshape([len(x), 1])
    M_best = models[np.argmin(AIC)]
    means = M_best.means_
    N = M_best.n_components
    std_devs = [ np.sqrt(  np.trace(M_best.covariances_[i])/N) for i in range(0,N) ]

    logprob = M_best.score_samples(x)
    responsibilities = M_best.predict_proba(x)

    p = np.exp(logprob)
    p_individual = responsibilities * p[:, np.newaxis]

    r = responsibilities
    r = r[:, np.argsort([u for m in means for u in m]) ]

    return( x, p, p_individual, means, std_devs, r )

def solve(m1,m2,std1,std2):
    """
    Function to calculate the intersection points of two distributions
    Input: Mean of one distribution
    Input: Mean of another distribution
    Input: Standard deviation of a distribution
    Input: Standard deviation of another distribution
    Return: Intersection points
    """
    a = 1/(2*std1**2) - 1/(2*std2**2)
    b = m2/(std2**2) - m1/(std1**2)
    c = m1**2 /(2*std1**2) - m2**2 / (2*std2**2) - np.log(std2/std1)
    return np.roots([a,b,c])

def create_vaf_plot():
    """
    Function to plot the VAF values
    """
    # Open a multipage pdf file
    with PdfPages(vcf_name+'.SMuRF.vafplot.pdf') as pdf:
        for sample in vaf_dict:
            if ( len(vaf_dict[sample]) <= 1):
                continue
            plt.figure(figsize=(30,10))

            x, p, p_individual, means, std_devs, respons = gmm( vaf_dict[sample] )

            for idx in range(0,len(x)):
                x_vaf = '{0:.2f}'.format(x[idx][0])
                responsibilities_dict[sample][x_vaf] = respons[idx].tolist()

            ax0 = plt.subplot(1,2,1)
            ax0.hist(vaf_dict[sample],bins=50)

            ax1 = plt.subplot(1,2,2)
            ax1.plot(x, p_individual, '--k', color='lightgray',linewidth=0.5)
            ax1.plot(x, p, '-k')
            ax1.axvline(x=0.5)
            ax1.axvline(x=0.3)
            mu = []

            # Calculate the intersection points between the distributions
            result = solve(means[0],means[1], std_devs[0], std_devs[1])
            for r in result:
                if r > 0.0 and r <= 1.0:
                    # ax1.axvline(x=r,linestyle='dashed',color='lightgray',linewidth=0.5)
                    means = np.append(means,r)


            # Plot lines at the mean of each distribution
            for m in means:
                mu.append(str(("{0:.2f}".format(float(m)))))
                ax1.axvline(x=m,linestyle='dashed',color='lightgray',linewidth=0.5)

            # Plot formatting
            plt.title(sample)
            # plt.xlabel('VAF')
            plt.ylabel('p(VAF)')

            # Set second x axis for the means of each component
            ax2 = ax1.twiny()
            ax2.xaxis.set_tick_params(length=15)
            ax2.set_xticks(means)
            ax2.set_xticklabels(mu)
            ax2.xaxis.set_ticks_position('bottom')
            ax2.xaxis.set_label_position('bottom')
            ax2.set_xlim(ax1.get_xlim())

            pdf.savefig()
            plt.close()

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
        update_call_data(call, ['FT'], [None], vcf_reader)

        # Check is the sample is in the control list
        if call.sample in args.normal:
            sample = False
        # QC fails if there is no genotype
        if call['GT'] == './.' or call['GT'] == '.|.':
            qc[sample][call.sample] = 'NoGenoType'
        # QC fails if the coverage is too low
        elif (call['DP'] == None or call['DP'] < int(cfg['SMuRF']['coverage'])):
            qc[sample][call.sample] = 'LowCov'
        elif sample:
            # If sample is homozygous reference
            if call['GT'] == '0/0' or call['GT'] == '0|0':
                noSampleEvidence += 1
                if indel and (not call['GQ'] or call['GQ'] < int(cfg['SMuRF']['indel_gq_homref'])):
                    qc[sample][call.sample] = 'LowGQ'
                elif not indel and (not call['GQ'] or call['GQ'] < int(cfg['SMuRF']['sample_gq_homref'])):
                    qc[sample][call.sample] = 'LowGQ'
                else:
                    qc[sample][call.sample] = 'PASS'
            # Check QC for homozygous variant
            elif call['GT'] == '1/1' or call['GT'] == '1|1':
                if indel and (not call['GQ'] or call['GQ'] < int(cfg['SMuRF']['indel_gq_homozygous'])):
                    qc[sample][call.sample] = 'LowGQ'
                elif not indel and (not call['GQ'] or call['GQ'] < int(cfg['SMuRF']['sample_gq_homozygous'])):
                    qc[sample][call.sample] = 'LowGQ'
                else:
                    qc[sample][call.sample] = 'PASS'
            elif call['GT'] == '0/1' or call['GT'] == '0|1':
                if indel and (not call['GQ'] or call['GQ'] < int(cfg['SMuRF']['indel_gq_heterozygous'])):
                    qc[sample][call.sample] = 'LowGQ'
                elif not indel and (not call['GQ'] or call['GQ'] < int(cfg['SMuRF']['sample_gq_heterozygous'])):
                    qc[sample][call.sample] = 'LowGQ'
                else:
                    qc[sample][call.sample] = 'PASS'
        else:
            # If variant is also found in a control
            if call['GT'] == '0/1' or call['GT'] == '0|1':
                controlEvidence = True
                if indel and (not call['GQ'] or call['GQ'] < int(cfg['SMuRF']['indel_gq_heterozygous'])):
                    qc[sample][call.sample] = 'LowGQ'
                elif not indel and (not call['GQ'] or call['GQ'] < int(cfg['SMuRF']['control_gq_heterozygous'])):
                    qc[sample][call.sample] = 'LowGQ'
                else:
                    qc[sample][call.sample] = 'isVariant'
            elif call['GT'] == '1/1' or call['GT'] == '1|1':
                controlEvidence = True
                if indel and (not call['GQ'] or call['GQ'] < int(cfg['SMuRF']['indel_gq_homozygous'])):
                    qc[sample][call.sample] = 'LowGQ'
                elif not indel and (not call['GQ'] or call['GQ'] < int(cfg['SMuRF']['control_gq_homozygous'])):
                    qc[sample][call.sample] = 'LowGQ'
                else:
                    qc[sample][call.sample] = 'isVariant'
            elif call['GT'] == '0/0' or call['GT'] == '0|0':
                if indel and (not call['GQ'] or call['GQ'] < int(cfg['SMuRF']['indel_gq_homref'])):
                    qc[sample][call.sample] = 'LowGQ'
                elif not indel and (not call['GQ'] or call['GQ'] < int(cfg['SMuRF']['control_gq_homref'])):
                    qc[sample][call.sample] = 'LowGQ'
                else:
                    qc[sample][call.sample] = 'PASS'
        if call.sample in qc[sample]:
            update_call_data(call,['FT'],[qc[sample][call.sample]], vcf_reader)
    format_list = list(vcf_reader.formats.keys())
    format_list.remove('GT')
    format_list.insert(0,'GT')
    # Add VAF information to the format field of each sample
    record.FORMAT = ":".join(format_list)
    # Flag variant if it is found in one of the controls
    if controlEvidence:
        record.FILTER.append('ControlEvidence')
        check = False
    # Flag variant if it is not found in one of the samples
    elif len(qc[True]) == noSampleEvidence :
        record.FILTER.append('NoSampleEvidence')
        check = False
    elif len(qc[False].keys()) > 0 and list(qc[False].values()).count('PASS') == 0:
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

def calculate_vaf( record ):
    """
    Function to calculate the VAF in the bam files
    Input: VCF record object
    Return: A tuple with VAF of each sample
    """
    record_vaf = {}
    vaf_info = collections.defaultdict(lambda: collections.defaultdict(list))
    qc = collections.defaultdict(dict)
    sample = "UNKNOWN"

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
        dx = 0
        vaf = 0.0
        # Loop through each reads that is overlapping the position of the variant
        for pileupcolumn in F.pileup(record.CHROM, int(record.POS)-1, int(record.POS), truncate=True, stepper='nofilter',min_base_quality=int(cfg['SMuRF']['base_phred_quality'])):
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
                        else:
                            dx+=1
                    # If variant is an insertion
                    elif ( len(record.REF) == 1 and len(alt) > 1 ):
                        # Read has the insertion
                        if ( pileupread.indel == len(alt)-1 ):
                            dv+=1
                        # Read has no insertion
                        elif pileupread.indel == 0:
                            dr+=1
                        else:
                            dx+=1
                    # If variant is an INDEL
                    else:
                        # Read has the INDEL
                        if ( pileupread.indel == (len(alt)-len(record.REF)) ):
                            dv+=1
                        # Read has no INDEL
                        elif pileupread.indel == 0:
                            dr+=1
                        else:
                            dx+=1
        # Calculate the VAF
        try:
            vaf = float("{0:.2f}".format(dv/float(dv+dr)))
        except ZeroDivisionError:
            continue
        # Loop through each sample in the vcf file
        for call in (record.samples):
            sample = True
            if call.sample in args.normal:
                sample = False
            # Check if the sample name in the vcf file is the same as a sample name in the bam file
            if call.sample == sample_name:
                # Add the VAF and sample name to the output tuple
                if vaf > float(cfg['SMuRF']['absent_threshold']) and call.sample not in args.normal:
                    record_vaf[call.sample] = vaf
                # Update the format field for this sample
                ft = call['FT']
                if float(dv+dr) < int(cfg['SMuRF']['coverage']):
            	     ft = 'LowCov'
            	     qc[sample][call.sample] = 'LowCov'
                else:
                    qc[sample][call.sample] = call['FT']
                update_call_data(call, ['VAF','CAD','FT'], [vaf, [dr, dv], ft], vcf_reader)

                # Set absent, subclonal or clonal based on the VAF and threshold
                if vaf <= float(cfg['SMuRF']['absent_threshold']):
                    vaf_info[sample]['ABSENT'].append(call.sample)
                elif vaf < float(cfg['SMuRF']['clonal_threshold']):
                    vaf_info[sample]['SUBCLONAL'].append(call.sample)
                else:
                    vaf_info[sample]['CLONAL'].append(call.sample)

    format_list = list(vcf_reader.formats.keys())
    format_list.remove('GT')
    format_list.insert(0,'GT')
    # Add VAF information to the format field of each sample
    record.FORMAT = ":".join(format_list)
    if not sample and dx > 0:
        record.FILTER.append('OverlappingIndelInControl')
    # Add QC information to the INFO field
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

        # Flag variant if it is found subclonal in a control
        if len(vaf_info[False]['SUBCLONAL']) > 0:
            record.FILTER.append('ControlSubclonal')
        # Flag variant if it is found clonal in a control
        elif len(vaf_info[False]['CLONAL']) > 0:
            record.FILTER.append('ControlClonal')
        # Flag variant if it is not found clonal in one of the samples
        elif len(vaf_info[True]['CLONAL']) == 0:
            record.FILTER.append('NoClonalSample')

    return( record_vaf )

def create_blacklist():
    """
    Function to fill the blacklist dictionary
    """
    # global blacklist
    # Loop through every blacklist file given on the command line
    blacklists = []
    for bl_vcf in cfg['SMuRF']['blacklist'].split(","):
        blacklist_single = pd.read_csv(bl_vcf,
                                           sep="\t",
                                           comment="#",
                                           header=None,
                                           names=["chr", "loc"],
                                           usecols=["chr", "loc"],
                                           dtype={0: "str", 1: "int"})
        if bl_vcf.endswith(".bed"):
            blacklist_single["loc"] += 1 #Bed files are 0-based and are converted to 1-based.
        elif not bl_vcf.endswith(".vcf"):
            warnings.warn("""The blacklist: {0} is not a .vcf or .bed file. Continuing with the following assumptions:\n
                          Column 1: Chromosome\n
                          Column 2: 1-based position\n
                          Header/Comments: #""".format(bl_vcf))

        blacklists.append(blacklist_single)

    blacklist = pd.concat(blacklists)
    blacklist["chr"] = blacklist["chr"].str.lower().str.replace("chr|chrom", "")
    blacklist = {k: g["loc"].tolist() for k, g in blacklist.groupby("chr")}

    return( blacklist )

def merge_tmp_vcfs():
    """
    Function to merge all the tmp contig vcf files
    """
    start = time.time()
    header = False
    # Loop through all chromomsomes
    for contig in contig_dict.keys():
        if not os.path.exists('SMuRF_tmp/{}.SMuRF.vcf'.format(contig)):
            continue
        if not header:
            os.system('cat SMuRF_tmp/{}.SMuRF.vcf > {}.SMuRF.vcf'.format(contig, vcf_name))
            header = True
        else:
            os.system('grep -v \'^#\' SMuRF_tmp/{}.SMuRF.vcf >> {}.SMuRF.vcf'.format(contig, vcf_name))
    os.system("grep -P '^#|\s+PASS\s+' "+vcf_name+".SMuRF.vcf > SMuRF_tmp/filter.vcf")

def add_responsibilities():
    with open('SMuRF_tmp/filter.vcf', 'r', encoding='utf-8') as vcf_file:
        vcf_reader = pyvcf.Reader(vcf_file)
    #vcf_reader = pyvcf.Reader(filename="SMuRF_tmp/filter.vcf", encoding='utf-8')
        vcf_reader.formats['PC'] = pyvcf.parser._Format('PC',None,'Float','Probability of each component')

        vcf_writer =  pyvcf.Writer(open(vcf_name+'.SMuRF.filtered.vcf','w', encoding='utf-8'), vcf_reader)

        for record in vcf_reader:
            for call in (record.samples):
                vaf = call['VAF']
                sample = call.sample
                if vaf != None and vaf > float(cfg['SMuRF']['absent_threshold']) and vaf in responsibilities_dict[sample]:
                    vaf = '{0:.2f}'.format(vaf)
                    update_call_data(call, ['PC'], [responsibilities_dict[sample][str(vaf)]], vcf_reader)
                else:
                    update_call_data(call, ['PC'], [None], vcf_reader)
            format_list = list(vcf_reader.formats.keys())
            format_list.remove('GT')
            format_list.insert(0,'GT')
            # Add VAF information to the format field of each sample
            record.FORMAT = ":".join(format_list)
            vcf_writer.write_record(record)
    # os.system("rm -rf SMuRF_tmp/*")
    time.sleep(10)
    os.system("rm -rf SMuRF_tmp")

if __name__ == "__main__":
    #get_command_line()
    main()
    create_vaf_plot()
    merge_tmp_vcfs()
    add_responsibilities()
