#!/usr/bin/python

import vcf as pyvcf
import argparse
import re
import os
import sys
import collections
import multiprocessing as mp
import queue
import time

contig_list = []


def main():
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
                done = q_out.get(block=False, timeout=1)
        except queue.Empty:
            pass
    # Give tasks a chance to put more data in
        time.sleep(10)
        if not q.empty():
            continue
        liveprocs = [p for p in liveprocs if p.is_alive()]

    for p in processes:
        p.join()


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

# Set arguments
parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-i', '--input', type=str, help='Input indexed vcf.gz file', required=True)
parser.add_argument('-c','--clonal', default=0.3,type=float,help='Clonal threshold (default: %(default)s)')
parser.add_argument('-a','--absent', default=0,type=float,help='Clonal threshold (default: %(default)s)')
parser.add_argument('-t','--threads', default=8,type=int,help='Number of threads')
args = parser.parse_args()

vcf_reader = pyvcf.Reader(filename=args.input, encoding='utf-8')
args.normal = []
for control in re.findall(r"--normal\s+(\S+)", vcf_reader.metadata['SMuRFCmd'][0]):
    args.normal.append(control)
vcf_name = os.path.basename(args.input)
vcf_name = vcf_name.replace(".vcf.gz","")
# Create tmp directory if it does not exists
try:
    os.stat('./SMuRF_tmp')
except:
    os.mkdir('./SMuRF_tmp')

def parse_chr_vcf(q, q_out, contig_vcf_reader):
    while True:
        try:
            contig = q.get(block=False,timeout=1)
            contig_vcf_flag_writer = pyvcf.Writer(open('./SMuRF_tmp/{}_SMuRF_reannotate.vcf'.format(contig),'w', encoding='utf-8'), contig_vcf_reader)
            vcf_reader.metadata['SMuRF_reannotate'] = [get_command_line()]
            try:
                # Try to parse the specific contig from the vcf
                contig_vcf_reader.fetch(contig)
            except:
                # Skip contig if it is not present in the vcf file
                continue
            for record in contig_vcf_reader.fetch(contig):
                for fil in ['NoClonalSample','ControlClonal','ControlSubclonal']:
                    if fil in record.FILTER:
                        record.FILTER.remove(fil)
                if not record.FILTER:
                    vaf_info = collections.defaultdict(lambda: collections.defaultdict(list))
                    for call in (record.samples):
                        sample = True
                        if call.sample in args.normal:
                            sample = False
                        vaf = call['VAF']
                        if vaf <= float(args.absent):
                            vaf_info[sample]['ABSENT'].append(call.sample)
                        elif vaf < float(args.clonal):
                            vaf_info[sample]['SUBCLONAL'].append(call.sample)
                        else:
                            vaf_info[sample]['CLONAL'].append(call.sample)
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
                contig_vcf_flag_writer.write_record(record)
            q_out.put( "Done" )

        # Break the loop if the queue is empty
        except queue.Empty:
            break

def merge_tmp_vcfs():
    """
    Function to merge all the tmp contig vcf files
    """
    start = time.time()
    header = False
    # Loop through all chromomsomes
    for contig in contig_list:
        if not header:
            os.system('cat SMuRF_tmp/{}_SMuRF_reannotate.vcf > {}_reannotate.vcf'.format(contig, vcf_name))
            header = True
        else:
            os.system('grep -v \'^#\' SMuRF_tmp/{}_SMuRF_reannotate.vcf >> {}_reannotate.vcf'.format(contig, vcf_name))
    os.system("grep -P '^#|\s+PASS\s+' "+vcf_name+"_reannotate.vcf > "+vcf_name+"_reannotate_filtered.vcf")
    time.sleep(5)
    os.system("rm -rf SMuRF_tmp")

if __name__ == "__main__":
    #get_command_line()
    main()
    merge_tmp_vcfs()
