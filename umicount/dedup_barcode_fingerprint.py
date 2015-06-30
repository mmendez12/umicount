"""
.. module:: fingerprint
   :platform: Unix
   :synopsis: Use UMI to count transcripts
   :version: 1.0

.. moduleauthor:: Mickael Mendez <mendez.mickael@gmail.com>

"""

import csv
import itertools
import subprocess
import argparse
import tempfile
import os
import shutil
from collections import defaultdict

import bed12


def get_fingerprint(read):
    """Get fingerprint id from the read's name. It assumes that the read's name
    contains the following pattern *FP:XXX;* where *XXX* is the fingerprint id.

    Args:
        read: A list of twelve elements where each element refers to a field in the BED format.

    Returns:
        A string containing the fingerprint id

    >>> read = ['chrX', '100', '200', 'BC:ATGC;FP:0012', '12', '+', '100', '110', '255,0,0', '2', '21,25', '0,75']
    >>> get_fingerprint(read)
    '0012'
    """
    return read[3].split('FP:')[1].split(';')[0]


def get_barcode(read):
    """Get barcode from the read's name. It assumes that the read's name
    contains the following pattern *BC:XXX;* where *XXX* is the actual barcode.

    Args:
        read: A list of twelve elements where each element refers to a field in the BED format.

    Returns:
        A string containing the barcode

    >>> read = ['chrX', '100', '200', 'BC:ATGC;FP:0012', '12', '+', '100', '110', '255,0,0', '2', '21,25', '0,75']
    >>> get_barcode(read)
    'ATGC'
    """
    return read[3].split('BC:')[1].split(';')[0]


def print_read_to_bed12(key, reads):
    """ Merge the reads by blocks and print a single read in the BED12 format on stdout.
    It assumes that the reads are on the same TSS and contains
    barcode and fingerprint information in the read's name.

    Args:
        key: A tuple that contain the chromosome, barcode and fingerprint information.

        reads: A list of reads (in a list) from the same TSS, that have similar barcode and fingerprint.

    >>> reads = []
    >>> reads.append(['chrX', '100', '200', 'BC:AAA;FP:0012', '12', '+', '100', '110', '255,0,0', '2', '20,25', '0,75'])
    >>> reads.append(['chrX', '100', '300', 'BC:AAA;FP:0012', '12', '+', '100', '110', '255,0,0', '3', '20,25', '0,175'])
    >>> print_read_to_bed12(('chrX', 'AAA', '0012'), reads) #doctest: +NORMALIZE_WHITESPACE
    chrX    100	300	BC:AAA;FP:0012	2	+	100	120	255,0,0	3	20,25,25	0,75,175
    """
    block_sizes, block_starts = bed12.merge_overlapping_blocks(reads)
        
    #bed12
    first_read = sorted(reads, key = bed12.get_start)[0]
    chrom, barcode, fingerprint = key
    start = bed12.get_start(first_read)
    end = start + block_starts[-1] + block_sizes[-1]
    name = "BC:{0};FP:{1}".format(barcode, fingerprint)
    score = len(reads)
    
    strand = bed12.get_strand(first_read)
    
    if strand == '+':
        thick_start = start
        thick_end = start + block_sizes[0]
    else:
        thick_start = end - block_sizes[-1]
        thick_end = end
        
    color = "255,0,0"
    block_count = len(block_sizes)
    block_sizes = ','.join(map(str, block_sizes))
    block_starts = ','.join(map(str, block_starts))
    
    output = [chrom, start, end, name, score, strand, thick_start, thick_end,
              color, block_count, block_sizes, block_starts]
    
    output_str = map(str, output)
    print '\t'.join(output_str)


def main():

    #PARSER TODO: move this code somewhere else
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-d", "--directory", help="absolute path of the folder containing the bed files")
    group.add_argument("-f", "--file", help="a bed file")
    parser.add_argument("-o", help='name of the output file. Only works if the script is called with the -f option, \
                                    ignored otherwise.')

    args = parser.parse_args()

    if args.directory:
        path, folder, files = os.walk(args.directory).next()
    elif args.file:
        path = ''
        files = [args.file]
    #ENDPARSER

    #create a temporary directory
    tmp_dir = tempfile.mkdtemp()

    plus_strand_tmp_file = open(os.path.join(tmp_dir, '+'), 'w')
    minus_strand_tmp_file = open(os.path.join(tmp_dir, '-'), 'w')
    plus_and_minus_sorted_path = os.path.join(tmp_dir, '+-s')

    #creates two temporary bed files containing either reads on the plus or minus strand
    for bed_file in files:

        with open(os.path.join(path, bed_file)) as bed_file:
            reader = csv.reader(bed_file, delimiter='\t')

            for read in reader:
                strand = bed12.get_strand(read)
                if strand == '+':
                    plus_strand_tmp_file.write('\t'.join(read) + '\n')
                elif strand == '-':
                    minus_strand_tmp_file.write('\t'.join(read) + '\n')


    #close the files
    plus_strand_tmp_file.close()
    minus_strand_tmp_file.close()

    #call unix sort on the file containing reads on the plus strand by tss
    with open(os.path.join(tmp_dir, '+sorted'), "w") as outfile:
        subprocess.call(["sort", '-k2,2n', os.path.join(tmp_dir, '+')], stdout=outfile)

    #call unix sort on the file containing reads on the minus strand by tss
    with open(os.path.join(tmp_dir, '-sorted'), "w") as outfile:
        subprocess.call(["sort", '-k3,3n', os.path.join(tmp_dir, '-')], stdout=outfile)

    #concatenate the files sorted by tss
    with open(plus_and_minus_sorted_path, "w") as outfile:
        subprocess.call(['cat', os.path.join(tmp_dir, '+sorted'), os.path.join(tmp_dir, '-sorted')], stdout=outfile)

    with open(plus_and_minus_sorted_path) as bedfile:
        reader = csv.reader(bedfile, delimiter='\t')
        reads = (line for line in reader)

        #for each reads on the same tss
        for tss, reads in itertools.groupby(reads, bed12.get_tss):
            d = defaultdict(list)

            #group the reads by chr, barcode and fingerprint
            for read in reads:
                key = (bed12.get_chrom(read), get_barcode(read), get_fingerprint(read))
                d[key].append(read)

            #merge and print the reads that have similar tss, barcode and fingerprint
            for key, reads in d.iteritems():
                print_read_to_bed12(key, reads)

    shutil.rmtree(tmp_dir)


if __name__ == '__main__':
    main()

#TODO: combine this script with dedup_fingerprint