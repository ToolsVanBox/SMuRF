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
parser.add_argument('-n', '--normal', action='append', type=str, help='Normal sample name, ignore for filtering')
parser.add_argument('-A', '--absent_pass_qc', default='n', type=str, help='Minimal number of PASS QC ABSENT samples (default: %(default)s)')
parser.add_argument('-C', '--clonal_pass_qc', default='n', type=str, help='Minimal number of PASS QC CLONAL samples (default: %(default)s)')
parser.add_argument('-a', '--absent', default=1, type=int, help='At least number of ABSENT samples (default: %(default)s)')
parser.add_argument('-c', '--clonal', default=1, type=int, help='At leat number of PASS samples (default: %(default)s)')

args = parser.parse_args()

vcf_reader = pyvcf.Reader(filename=args.input, encoding='utf-8')
total_samples = vcf_reader.samples
print("#"+get_command_line())
print("pos\t"+"\t".join(total_samples))

for record in vcf_reader:
    absent_sample_names = record.INFO['ABSENT_SAMPLE_NAMES']
    clonal_sample_names = record.INFO['CLONAL_SAMPLE_NAMES']
    subclonal_sample_names = record.INFO['SUBCLONAL_SAMPLE_NAMES']
    pass_qc_sample_names = record.INFO['PASS_QC_SAMPLE_NAMES']
    fail_qc_sample_names = record.INFO['FAIL_QC_SAMPLE_NAMES']

    if args.normal:
        for normal in args.normal:
            if normal in absent_sample_names:
                absent_sample_names.remove(normal)
            if normal in clonal_sample_names:
                clonal_sample_names.remove(normal)
            if normal in pass_qc_sample_names:
                pass_qc_sample_names.remove(normal)
            if normal in fail_qc_sample_names:
                fail_qc_sample_names.remove(normal)
            if normal in subclonal_sample_names:
                subclonal_sample_names.remove(normal)

    if '' in absent_sample_names:
        absent_sample_names.remove('')
    if '' in clonal_sample_names:
        clonal_sample_names.remove('')
    if '' in pass_qc_sample_names:
        pass_qc_sample_names.remove('')


    if len(absent_sample_names) < args.absent:
        continue
    if len(clonal_sample_names) < args.clonal:
        continue

    absent_pass_qc_sample_names = list(set(pass_qc_sample_names) & set(absent_sample_names))
    clonal_pass_qc_sample_names = list(set(pass_qc_sample_names) & set(clonal_sample_names))

    ap = False
    cp = False
    if args.absent_pass_qc == 'n':
        if len(absent_pass_qc_sample_names) == len(absent_sample_names):
            ap = True
    else:
        if len(absent_pass_qc_sample_names) >= int(args.absent_pass_qc):
            ap = True
    if args.clonal_pass_qc == 'n':
        if len(clonal_pass_qc_sample_names) == len(clonal_sample_names):
            cp = True
    else:
        if len(clonal_pass_qc_sample_names) >= int(args.clonal_pass_qc):
            cp = True

    if ap and cp:
        outline = [record.CHROM+":"+str(record.POS)]
        for sample in total_samples:
            label = "NA"
            if sample in absent_sample_names:
                label="A"
            elif sample in subclonal_sample_names:
                label="S"
            elif sample in clonal_sample_names:
                label="C"
            if sample in pass_qc_sample_names:
                label+="P"
            elif sample in fail_qc_sample_names:
                label+="F"
            outline.append(label)
        print( "\t".join(outline))
