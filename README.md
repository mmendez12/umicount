#Introduction

umicount is a collection of Python scripts which allows to remove and count PCR duplicates from paired-end libraries
prepared with unique molecular identifiers (UMIs).

The main difference between existing approaches (rmdup or MarkDuplicates) is that it uses UMI and Transcription Start Site information to remove duplicates rather than the reads size.

It was mainly developed for single-cell CAGE and single-cell nanoCAGE protocols where a tagmentation step is performed
between two PCRs (see example below).

Note that umicount is designed to work on large sequencing data and only requires a small amount of memory. So far it was
tested on MiSeq and HiSeq data without problems.

For more information about umicount, you can look at the doc [here](http://umicount.readthedocs.org/en/latest/).
