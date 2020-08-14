#!/usr/bin/python

import vcf as pyvcf
import collections
import argparse
import multiprocessing as mp
import queue
#import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pysam
import time


# Set arguments
parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-i', '--input', type=str, help='Input indexed vcf.gz file', required=True)
parser.add_argument('-b', '--bam', action='append', nargs="*", type=str, help='Input bam file', required=True)
parser.add_argument('-t', '--threads', default=8, type=int, help="Number of threads (default: %(default)s)")
parser.add_argument('-f', '--filter', default=[], action='append', nargs="*", type=str, help='Plot only variants with FILTER flag (default: %(default)s)')
parser.add_argument('-v', '--vaf', default=0, type=float, help="Plot only variants with a VAF > threshold (default: %(default)s)")
parser.add_argument('-m', '--mapq', default=0, type=int, help="Include only reads with a minimal mapq (default: %(default)s)")
parser.add_argument('-p', '--base_phred_quality', default=0, type=int, help="Include only bases with a minimal base phred quality (default: %(default)s)")
parser.add_argument('-indel', '--indel', default=True, action='store_false', help="Exclude indels")
args = parser.parse_args()

vcf_reader = pyvcf.Reader(filename=args.input, encoding='utf-8')

vaf_dict = collections.defaultdict(list)
bam_sample_names = collections.defaultdict(dict)
contig_list = []

args.bam = [x for l in args.bam for x in l]
args.filter = [x for l in args.filter for x in l]

def main():
    global vcf_reader, vaf_df, blacklist
    vcf_reader = fix_vcf_header(vcf_reader)

    for contig in vcf_reader.contigs:
        # if contig != '21':
        #     continue
        contig_list.append(contig)

    # Create an input queue with the contigs and an empty output queue
    q = mp.Queue()
    q_out = mp.Queue()
    for contig in contig_list:
        q.put(contig)

    # Create number of processes to parse the vcf file
    processes = [mp.Process(target=parse_chr_vcf, args=(q, q_out, vcf_reader, args.bam)) for x in range(args.threads)]

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

def parse_chr_vcf(q, q_out, contig_vcf_reader, bams):
    """
    Function to parse the vcf per contig.
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
            try:
                # Try to parse the specific contig from the vcf
                contig_vcf_reader.fetch(contig)
            except:
                # Skip contig if it is not present in the vcf file
                continue
            for record in contig_vcf_reader.fetch(contig):
                filter_flag = True
                for filter in record.FILTER:
                    if filter not in args.filter:
                        filter_flag = False
                        break
                if not filter_flag:
                    continue
                vaf = False
                for call in record.samples:
                    try:
                        if call['VAF'] is not None:
                            contig_vaf[call.sample].append(call['VAF'])
                            vaf = True
                    except:
                        break
                if not vaf:
                    for s,v in calculate_vaf( record ).items():
                        contig_vaf[s].append(v)
            q_out.put( contig_vaf )

        # Break the loop if the queue is empty
        except queue.Empty:
            break


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

def calculate_vaf( record ):
    """
    Function to calculate the VAF in the bam files
    Input: VCF record object
    Return: A tuple with VAF of each sample
    """
    record_vaf = {}
    vaf_info = collections.defaultdict(lambda: collections.defaultdict(list))
    qc = collections.defaultdict(dict)

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
            # Check if the sample name in the vcf file is the same as a sample name in the bam file
            if call.sample == sample_name:
                # Add the VAF and sample name to the output tuple
                if vaf > args.vaf:
                    record_vaf[call.sample] = vaf

    return( record_vaf )

def create_vaf_plot():
    """
    Function to plot the VAF values
    """
    # Open a multipage pdf file
    with PdfPages('VAFplot.pdf') as pdf:
        for sample in vaf_dict:
            plt.figure(figsize=(30,10))

            plt.hist(vaf_dict[sample],bins=50)

            # Plot formatting
            plt.title(sample)
            # plt.xlabel('VAF')
            plt.ylabel('p(VAF)')

            pdf.savefig()
            plt.close()



if __name__ == "__main__":
    #get_command_line()
    main()
    create_vaf_plot()
