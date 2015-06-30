"""
.. module:: cage_scan_clustering
   :platform: Unix
   :synopsis: cluster cage scan reads based on arbitrary distance. It assumes that the reads are sorted by chrom, TSS and strand.
   :version: 1.0

.. moduleauthor:: Mickael Mendez <mendez.mickael@gmail.com>

"""

import csv
import argparse


import bed12

#TODO: rewrite tests
def print_read_to_bed12(reads):
    """ Merge the reads by blocks and print a single read in the BED12 format on stdout.
    It assumes that the reads are on the same TSS and contains
    fingerprint information in the read's name.

    Args:
        reads: A list of reads

    """
    block_sizes, block_starts = bed12.merge_overlapping_blocks(reads)
        
    #bed12
    first_read = sorted(reads, key=bed12.get_start)[0]
    chrom = bed12.get_chrom(first_read)
    start = bed12.get_start(first_read)
    end = start + block_starts[-1] + block_sizes[-1]

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

    name = map(str, [chrom, start, end, strand])
    name = ":".join(name)
    
    output = [chrom, start, end, name, score, strand, thick_start, thick_end,
              color, block_count, block_sizes, block_starts]
    
    output_str = map(str, output)
    print '\t'.join(output_str)


def overlapping_reads(reads, distance):
    """returns all the overlapping reads within a given distance"""

    reads_list = []
    cur_tss = 0
    cur_chrom = ''

    for read in reads:

        if not cur_tss:
            cur_tss = bed12.get_tss(read)
            reads_list.append(read)
            cur_chrom = bed12.get_chrom(read)
            continue


        tss = bed12.get_tss(read)
        chrom = bed12.get_chrom(read)

        #if not overlap
        if (tss - cur_tss > distance) or (chrom != cur_chrom):
            yield reads_list
            reads_list = [read]
            cur_tss = tss
            cur_chrom = chrom
        else:
            reads_list.append(read)
            cur_tss = tss

    yield reads_list


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('bed_file', help='input file')
    parser.add_argument('-t', '--tag_distance', default=20, type=int, help='cluster all the cage tags at distance d')

    args = parser.parse_args()

    with open(args.bed_file) as bedfile:

        reader = csv.reader(bedfile, delimiter='\t')
        reads = (line for line in reader)

        #for each reads on the same tss
        for read_list in overlapping_reads(reads, args.tag_distance):
            print_read_to_bed12(read_list)


if __name__ == '__main__':
    main()

#TODO: combine this script with fingerprint.py