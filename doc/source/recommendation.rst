.. _recommendation-section:

Recommendations
===============

.. _tagdust2: http://tagdust.sourceforge.net/
.. _doc: https://github.com/mmendez12/bedtools2/blob/pairedbamtobed12-pullRequest2/docs/content/tools/pairedbamtobed12.rst
.. _pairedbamtobed12: https://github.com/mmendez12/bedtools2/tree/pairedbamtobed12-pullRequest2
.. _BWA: http://bio-bwa.sourceforge.net/


umicount works with BED12 files. Here we quickly describe the softwares we tested to process and align the raw reads,
and generate BED12 files.

At first we use `tagdust2`_ to extract the UMIs from the raw reads. It has the advantage to work with paired-end data and keeps track
UMIs in the description field of the FASTQ files.

Note that some aligners remove the description field. In our case we need it since it contains a reference to the UMIs. One
workaround is to replace all the spaces in the FASTQ files by a custom separator with sed for example.

Then we align reads to a reference Genome with `BWA`_. This will work well for short reads.
For longer reads, one can simply trim them or use a different aligner, it should not be a problem as long as the description field is not lost.

Finally, we use `pairedbamtobed12`_ to convert properly paired BAM alignments to BED12 format. (read the `doc`_ to learn more about it)