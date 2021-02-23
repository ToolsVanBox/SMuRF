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
    for vcf_reader in vcf_readers:
        for contig in vcf_reader.contigs:
            if contig not in contig_list:
                contig_list.append(contig)

    # Create an input queue with the contigs and an empty output queue
    q = mp.Queue()
    q_out = mp.Queue()
    for contig in contig_list:
        q.put(contig)

    # Create number of processes to parse the vcf file
    processes = [mp.Process(target=parse_chr_vcfs, args=(q, q_out, vcf_readers)) for x in range(args.threads)]

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
parser.add_argument('-i', '--input', action='append', type=str, help='Input indexed SMuRF.vcf.gz file', required=True)
parser.add_argument('-o', '--output', default="SMuRF", type=str, help='Name of output file (default: %(default)s)')
parser.add_argument('-t','--threads', default=8,type=int,help='Number of threads (default: %(default)s)')
args = parser.parse_args()

vcf_readers = []
header_samples = []
sample_names = []
normal_names = []
for input_file in args.input:
    vcf_reader = pyvcf.Reader(filename=input_file, encoding='utf-8')
    sample_names.extend( vcf_reader.samples )
    header_samples.extend( vcf_reader.samples )
    for control in re.findall(r"--normal\s+(\S+)", vcf_reader.metadata['SMuRFCmd'][0]):
        if control != "None":
            normal_names.append(control)
        if control in sample_names:
            sample_names.remove(control)
    vcf_readers.append(vcf_reader)
vcf_name = os.path.basename(args.output)
#vcf_name = vcf_name.replace(".vcf.gz","")
# Create tmp directory if it does not exists
try:
    os.stat('./SMuRF_tmp')
except:
    os.mkdir('./SMuRF_tmp')

def parse_chr_vcfs(q, q_out, contig_vcf_readers):
    while True:
        try:
            contig = q.get(block=False,timeout=1)
            #contig_vcf_flag_writer = pyvcf.Writer(open('./SMuRF_tmp/{}_SMuRF_reannotate.vcf'.format(contig),'w', encoding='utf-8'), contig_vcf_reader)
            #vcf_reader.metadata['SMuRF_merge'] = [get_command_line()]
            merged_records = {}
            for contig_vcf_reader in contig_vcf_readers:
                record_dict = get_record_dict( contig, contig_vcf_reader)
                if not record_dict:
                    continue
                if not merged_records:
                    merged_records = record_dict
                else:
                    merged_records = mergeDict( merged_records, record_dict )
            if merged_records:
                merge_records(merged_records, contig, contig_vcf_readers[0])

            q_out.put( "Done" )

        # Break the loop if the queue is empty
        except queue.Empty:
            break

def merge_records( merged_records, contig, contig_vcf_reader):
    contig_vcf_reader.samples = header_samples
    contig_vcf_flag_writer = pyvcf.Writer(open('./SMuRF_tmp/{}_SMuRF_merged.vcf'.format(contig),'w', encoding='utf-8'), contig_vcf_reader)

    for key in merged_records:
        if type(merged_records[key]) is list:
            record = merged_records[key][0]
            for record2 in merged_records[key][1:]:
                if 'MERGED' not in record.INFO:
                    record.INFO['MERGED'] = 0
                record.INFO['MERGED'] += 1
                if record2.INFO['ABSENT_SAMPLE_NAMES'] != ['']:
                    record.INFO['ABSENT_SAMPLE_NAMES'].extend(record2.INFO['ABSENT_SAMPLE_NAMES'])
                    record.INFO['ABSENT_SAMPLES'] += record2.INFO['ABSENT_SAMPLES']
                if record2.INFO['SUBCLONAL_SAMPLE_NAMES'] != ['']:
                    record.INFO['SUBCLONAL_SAMPLE_NAMES'].extend(record2.INFO['SUBCLONAL_SAMPLE_NAMES'])
                    record.INFO['SUBCLONAL_SAMPLES'] += record2.INFO['SUBCLONAL_SAMPLES']
                if record2.INFO['CLONAL_SAMPLE_NAMES'] != ['']:
                    record.INFO['CLONAL_SAMPLE_NAMES'].extend(record2.INFO['CLONAL_SAMPLE_NAMES'])
                    record.INFO['CLONAL_SAMPLES'] += record2.INFO['CLONAL_SAMPLES']
                if record2.INFO['ABSENT_CONTROL_NAMES'] != ['']:
                    record.INFO['ABSENT_CONTROL_NAMES'].extend(record2.INFO['ABSENT_CONTROL_NAMES'])
                    record.INFO['ABSENT_CONTROLS'] += record2.INFO['ABSENT_CONTROLS']
                if record2.INFO['SUBCLONAL_CONTROL_NAMES'] != ['']:
                    record.INFO['SUBCLONAL_CONTROL_NAMES'].extend(record2.INFO['SUBCLONAL_CONTROL_NAMES'])
                    record.INFO['SUBCLONAL_CONTROLS'] += record2.INFO['SUBCLONAL_CONTROLS']
                if record2.INFO['CLONAL_CONTROL_NAMES'] != ['']:
                    record.INFO['CLONAL_CONTROL_NAMES'].extend(record2.INFO['CLONAL_CONTROL_NAMES'])
                    record.INFO['CLONAL_CONTROLS'] += record2.INFO['CLONAL_CONTROLS']
                if record2.INFO['FAIL_QC_SAMPLE_NAMES'] != ['']:
                    record.INFO['FAIL_QC_SAMPLE_NAMES'].extend(record2.INFO['FAIL_QC_SAMPLE_NAMES'])
                    record.INFO['FAIL_QC_SAMPLES'] += record2.INFO['FAIL_QC_SAMPLES']
                if record2.INFO['PASS_QC_SAMPLE_NAMES'] != ['']:
                    record.INFO['PASS_QC_SAMPLE_NAMES'].extend(record2.INFO['PASS_QC_SAMPLE_NAMES'])
                    record.INFO['PASS_QC_SAMPLES'] += record2.INFO['PASS_QC_SAMPLES']
                if 'FAIL_QC_CONTROL_NAMES' in record2.INFO and 'FAIL_QC_CONTROL_NAMES' in record.INFO and record2.INFO['FAIL_QC_CONTROL_NAMES'] != ['']:
                    record.INFO['FAIL_QC_CONTROL_NAMES'].extend(record2.INFO['FAIL_QC_CONTROL_NAMES'])
                    record.INFO['FAIL_QC_CONTROLS'] += record2.INFO['FAIL_QC_CONTROLS']
                if 'PASS_QC_CONTROL_NAMES' in record2.INFO and 'PASS_QC_CONTROL_NAMES' in record.INFO and record2.INFO['PASS_QC_CONTROL_NAMES'] != ['']:
                    record.INFO['PASS_QC_CONTROL_NAMES'].extend(record2.INFO['PASS_QC_CONTROL_NAMES'])
                    record.INFO['PASS_QC_CONTROLS'] += record2.INFO['PASS_QC_CONTROLS']
                record.samples.extend(record2.samples)
        else:
            record = merged_records[key]
            format_field = record.FORMAT.split(":")
            genotype = {}
            for f in format_field:
                if f == "GT":
                    genotype[f] = './.'
                else:
                    genotype[f] = None
            samples_data = []
            for sample in header_samples:
                if sample in sample_names:
                    if sample not in record.INFO['CLONAL_SAMPLE_NAMES'] and sample not in record.INFO['SUBCLONAL_SAMPLE_NAMES'] and sample not in record.INFO['ABSENT_SAMPLE_NAMES']:
                        record.INFO['ABSENT_SAMPLE_NAMES'].append(sample)
                        record.INFO['ABSENT_SAMPLES'] += 1
                        call = pyvcf.model._Call('site', sample, collections.namedtuple('CallData', format_field)(**genotype))
                    else:
                        call = pyvcf.model._Call('site', sample, record.genotype(sample).data)
                    samples_data.append(call)
                if sample in normal_names:
                    if sample not in record.INFO['CLONAL_CONTROL_NAMES'] and sample not in record.INFO['SUBCLONAL_CONTROL_NAMES'] and sample not in record.INFO['ABSENT_CONTROL_NAMES']:
                        record.INFO['ABSENT_CONTROL_NAMES'].append(sample)
                        record.INFO['ABSENT_CONTROLS'] += 1
                        call = pyvcf.model._Call('site', sample, collections.namedtuple('CallData', format_field)(**genotype))
                    else:
                        call = pyvcf.model._Call('site', sample, record.genotype(sample).data)
                    samples_data.append(call)
            record.samples = samples_data
        if record.INFO['CLONAL_SAMPLES'] == 0:
            record.FILTER.append('NoClonalSample')
        if record.INFO['CLONAL_CONTROLS'] > 0:
            record.FILTER.append('ControlClonal')
        if record.INFO['SUBCLONAL_CONTROLS'] > 0:
            record.FILTER.append('ControlSubclonal')
            # if not record.FILTER:
        contig_vcf_flag_writer.write_record(record)

def mergeDict(dict1, dict2):
   ''' Merge dictionaries and keep values of common keys in list'''
   dict3 = {**dict1, **dict2}
   for key, value in dict3.items():
       if key in dict1 and key in dict2:
               dict3[key] = [value , dict1[key]]
   return dict3

def get_record_dict(contig, contig_vcf_reader):
    record_dict = {}
    try:
        # Try to parse the specific contig from the vcf
        contig_vcf_reader.fetch(contig)
    except:
        # Skip contig if it is not present in the vcf file
        return( False )
    for record in contig_vcf_reader.fetch(contig):
        for fil in ['NoClonalSample','ControlClonal','ControlSubclonal']:
            if fil in record.FILTER:
                record.FILTER.remove(fil)
        if not record.FILTER:
            alt = ",".join(list(map(str, record.ALT)))
            record_dict[(contig,record.POS)] = record
    return( record_dict )


def merge_tmp_vcfs():
    """
    Function to merge all the tmp contig vcf files
    """
    start = time.time()
    header = False
    # Loop through all chromomsomes
    for contig in contig_list:
        if not header:
            os.system('cat SMuRF_tmp/{}_SMuRF_merged.vcf > {}_merged.vcf'.format(contig, vcf_name))
            header = True
        else:
            os.system('grep -v \'^#\' SMuRF_tmp/{}_SMuRF_merged.vcf >> {}_merged.vcf'.format(contig, vcf_name))
    os.system("grep -P '^#|\s+PASS\s+' "+vcf_name+"_merged.vcf > "+vcf_name+"_merged_filtered.vcf")
    time.sleep(5)
    #os.system("rm -rf SMuRF_tmp")

if __name__ == "__main__":
    #get_command_line()
    main()
    merge_tmp_vcfs()
