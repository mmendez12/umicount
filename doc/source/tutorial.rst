How to use umicount
====================

This section explains how to use the umicount to remove PCR duplicates from a BED12 file.

First you will have to make sure that umicount is correctly installed. See :ref:`install-section` for more information.

If you don't know how to generate a valid BED12 file for umicount, you can refer to the :ref:`recommendation-section`.

Working with a single bed file
******************************

Just specify the path to a BED12 formatted file with the *-f* argument::

    umicountFP -f ./your_bed12_file.bed


Working with multiple files
***************************

Just specify the path to the folder containing your BED12 formatted files with the *-d* argument::

    umicountFP -d ./path/to/your/bed12_formatted_files/

Print usage
***********

::

    umicountFP -h
