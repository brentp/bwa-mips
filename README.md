bwa-mips
========

Map sequence from Molecular Inversion Probes with BWA, strip arms, de-dup, ..., profit


`python ../bwamips.py $REF mips-design.txt $R1 $R2` > $SAM



With molecular inversion probes, we map reads to the genome that include the
ligation and extension arms, along with a molecular tag (AKA UMI).
This script takes:

    1) ref.fasta
    2) mips design file
    3) de-multiplxed, paired-end fastqs

and moves the UMI into the read-name, aligns the reads using bwa mem,
strips that ligation and extension arms from the alignment, adjusts the
position, cigar, sequence, and quality accounting (more or less) for
insertions, deletions, and masked sequence.

The output is in SAM format with extra tags as described
in the @CO headers.

Thanks to Evan Boyle and Jay Hesselberth for (repeated) explanations.
Any mistakes are my own.

bpederse@gmail.com

LICENSE is MIT


Installation
============

------------

`bwa-meth` depends on

 + python 2.7+ 
   - `toolshed` library. can be installed with:
      * `easy_install toolshed` or
      * `pip install toolshed`

   - if you don't have root or sudo priviledges, you can run
     `python setup.py install --user` from this directory and the bwameth.py
     executable will be at: ~/.local/bin/bwameth.py

   - if you do have root or sudo run: `[sudo] python setup.py install` from
     this directory
     
   - users unaccustomed to installing their own python packages should 
     download anaconda: https://store.continuum.io/cshop/anaconda/ and
     then install the toolshed module with pip as described above.
     
 + samtools command on the `$PATH` (https://github.com/samtools/samtools)
     
 + bwa mem from: https://github.com/lh3/bwa

Example
=======

see the example/ directory.

Report
======

`bwamips.py` outputs a report it looks like this for the example dataset

    ================
    Alignment Report
    ================

    Bases in genome: 171,115,067
    Bases in target region: 85,128 (0.0497% of genome)


    MIPs found: 13,636
    Off-target reads: 8,882
    Unmapped reads: 8,053

    Observed / expected enrichment where expected is based on size
    of target region relative to size of genome.

    FOLD ENRICHMENT
    ===============
    low - high: 1618.52 - 3085.97

    Low estimate uses unmapped reads that same as off-target.

    % READS ON TARGET
    =================
    80.52%

