#!/usr/bin/env python3
import argparse
import os
import glob
import sys

import spaTyper.spa_typing
import spaTyper.utils

#####################################################
parser = argparse.ArgumentParser(prog='spaTyper.py', formatter_class=argparse.RawDescriptionHelpFormatter, 
                                 description='''
spaTyper.py: Get spa types

Version: 0.2.0
License: GPLv3

USAGE: python spaTyper.py fasta_file.fa
Prints spa type to stdout

Will download sparepeats.fasta and spatypes.txt to repository directory if files not provided or already in directory.

''', epilog="Original code: mjsull. Modified by: JFSanchezHerrero")

#####################################################
parser.add_argument('-r', '--repeat_file', action='store', 
                    help='List of spa repeats (http://spa.ridom.de/dynamic/sparepeats.fasta)')

parser.add_argument('-o', '--repeat_order_file', action='store', 
                    help='List spa types and order of repeats (http://spa.ridom.de/dynamic/spatypes.txt)')

parser.add_argument('-d', '--folder', action='store', help='Folder to save downloaded files from Ridom/Spa server')


parser.add_argument('-f', '--fasta', action='store', nargs='+', 
                    help='List of one or more fasta files.')

parser.add_argument('-g', '--glob', action='store', 
                    help='Uses unix style pathname expansion to run spa typing on all files. '
                        'If your shell autoexpands wildcards use -f.')

parser.add_argument("-e", '--do_enrich', action="store_true", default=False, 
                    help="Do PCR product enrichment. [Default: False]")

parser.add_argument("-c", "--clean_output", action="store_true", default=False, 
                    help="Make output clean")

parser.add_argument('--version', action='version', version='%(prog)s 0.2.0')

args = parser.parse_args()
#####################################################

## Let's go
if not args.glob is None and not args.fasta is None:
    sys.exit('Please provide script with either a list of one or more fasta files, or a glob.')
elif not args.glob is None:
    fasta_list = glob.glob(args.glob)
elif not args.fasta is None:
    fasta_list = args.fasta
else:
    sys.exit('Please provide spaTyper.py with either a fasta file (-f) or glob (-g).')

print ("Start the identification of repeats in Spa protein:")

### Check if sparepeats and spatypes files are provided and available 
## or download them from SeqNet/Ridom Spa Server

## Sparepeats file
if (args.repeat_file):
    args.repeat_file = os.path.abspath(args.repeat_file)
    print ("+ Repeats fasta file provided: ", args.repeat_file)
else:
    if not (args.folder):
        args.folder =  os.path.dirname(os.path.realpath(__file__))
    else:
        args.folder = os.path.abspath(args.folder)
        print ("+ Create folder: ", args.folder)
        spaTyper.utils.create_folder(args.folder)
    
    # print message
    print ("+ Check or download repeats fasta file in folder: ", args.folder)
    
    ## download file repeats   
    args.repeat_file = spaTyper.utils.download_file_repeats(args.folder)

## spatypes file
if (args.repeat_order_file):
    args.repeat_order_file = os.path.abspath(args.repeat_order_file)
    print ("+ Repeat types file provided: ", args.repeat_order_file)
else:
    if not (args.folder):
        args.folder =  os.path.dirname(os.path.realpath(__file__))
    
    # print message
    print ("+ Check or download repeats types file in folder: ", args.folder)
    ## download file repeats   
    args.repeat_order_file = spaTyper.utils.download_file_types(args.folder)

## Get the SpaTypes in fasta sequences
## getSpaTypes
seqDict, letDict, typeDict, seqLengths = spaTyper.spa_typing.getSpaTypes(args.repeat_file, args.repeat_order_file)

## findPatterns
if not args.clean_output:
    sys.stdout.write('FILENAME\tEGENOMICS_SPA_TYPE\tRIDOM_SPA_TYPE\n')

for i in fasta_list:
    the_out = spaTyper.spa_typing.findPattern(i, seqDict, letDict, typeDict, seqLengths, args.do_enrich)
    if not args.clean_output:
        sys.stdout.write('Spa type:\t' + '\t'.join(the_out) + '\n')
    else:
        sys.stdout.write('\t'.join(the_out) + '\n')                
