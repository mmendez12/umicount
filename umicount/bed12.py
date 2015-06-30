"""
.. module:: bed12
   :platform: Unix
   :synopsis: Defines a set a generic function to parse and process bed12 files.

.. moduleauthor:: Mickael Mendez <mendez.mickael@gmail.com>

"""

__author__ = 'mickael'


import operator
import itertools

def get_chrom(read):
    """Get chromosome from a bed12 line.

    Args:
        read: A list of twelve elements where each element refers to a field in the BED format.

    Returns:
        The chromosome name

    >>> read = ['chrX', '100', '200', 'toto', '12', '+', '100', '110', '255,0,0', '2', '21,25', '0,75']
    >>> get_chrom(read)
    'chrX'
    """
    return read[0]


def get_start(read):
    """Get start position from a bed12 line.

    Args:
        read: A list of twelve elements where each element refers to a field in the BED format.

    Returns:
        An integer representing the start position of the read.

    >>> read = ['chrX', '100', '200', 'toto', '12', '+', '100', '110', '255,0,0', '2', '21,25', '0,75']
    >>> get_start(read)
    100
    """
    return int(read[1])


def get_end(read):
    """Get end position from a bed12 line.

    Args:
        read: A list of twelve elements where each element refers to a field in the BED format.

    Returns:
        An integer representing the end position of the read.

    >>> read = ['chrX', '100', '200', 'toto', '12', '+', '100', '110', '255,0,0', '2', '21,25', '0,75']
    >>> get_end(read)
    200
    """
    return int(read[2])


def get_strand(read):
    """Get strand from a bed12 line.

    Args:
        read: A list of twelve elements where each element refers to a field in the BED format.

    Returns:
        A single char representing the strand of a read

    >>> read = ['chrX', '100', '200', 'toto', '12', '+', '100', '110', '255,0,0', '2', '21,25', '0,75']
    >>> get_strand(read)
    '+'
    """
    return read[5]

def get_tss(read):
    """Get Transcription Start Site (TSS) from a bed12 line.

    Args:
        read: A list of twelve elements where each element refers to a field in the BED format.

    Returns:
        The start position as an integer if the read is on the plus strand.
        The end position as an integer if the read is on the minus strand.

    >>> read = ['chrX', '100', '200', 'toto', '12', '+', '100', '110', '255,0,0', '2', '21,25', '0,75']
    >>> get_tss(read)
    100
    >>> read = ['chrX', '100', '200', 'toto', '12', '-', '100', '110', '255,0,0', '2', '21,25', '0,75']
    >>> get_tss(read)
    200
    """
    strand = get_strand(read)

    if strand == '+':
        return get_start(read)
    else:
        return get_end(read)


def blocks_to_absolute_start_end(read):
    """Calculate the absolute start and end of the blocks from a bed12 line.

    Args:
        read: A list of twelve elements where each element refers to a field in the BED format.

    Returns:
        A list of tuple where each tuple contains the absolute start and end coordinates of a block.

    >>> read = ['chrX', '100', '200', 'toto', '12', '+', '100', '110', '255,0,0', '2', '21,25', '0,75']
    >>> blocks_to_absolute_start_end(read)
    [(100, 121), (175, 200)]
    """
    read_start = get_start(read)

    block_starts = [read_start + int(start) for start in read[11].split(',') if start]
    block_sizes = [int(size) for size in read[10].split(',') if size]

    block_starts_sizes = zip(block_starts, block_sizes)

    return [(bstart, bstart + bsize) for bstart, bsize in block_starts_sizes]


def merge_overlapping_blocks(reads):
    """Merge blocks if they overlap.

    Args:
        reads: A list of read in the BED12 format.

    Returns:
        Two lists where the first list contains the blocks sizes and the second the blocks starts.
        Values in the lists are integer.

    >>> reads = []
    >>> reads.append(['chrX', '100', '200', 'toto', '12', '+', '100', '110', '255,0,0', '2', '20,25', '0,75'])
    >>> reads.append(['chrX', '100', '200', 'toto', '12', '+', '100', '110', '255,0,0', '3', '10,10,25', '0,15,75'])
    >>> merge_overlapping_blocks(reads)
    ([25, 25], [0, 75])
    """

    blocks_list = [blocks_to_absolute_start_end(read) for read in reads]

    #flatten
    blocks = itertools.chain.from_iterable(blocks_list)

    final_blocks = []

    blocks = sorted(blocks, key = operator.itemgetter(0, 1))
    known_block_start, known_block_end = blocks[0]

    for block_start, block_end in blocks[1:]:
        if block_start <= known_block_end:
            known_block_end = max(block_end, known_block_end)
        else:
            final_blocks.append((known_block_start, known_block_end))
            known_block_start, known_block_end = (block_start, block_end)

    final_blocks.append((known_block_start, known_block_end))

    absolute_block_start = final_blocks[0][0]

    bsizes = [end - start for start, end in final_blocks]
    bstarts = [start - absolute_block_start for start, end in final_blocks]

    return bsizes, bstarts
