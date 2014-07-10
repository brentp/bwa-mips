bwa-mips
========

Map sequence from Molecular Inversion Probes with BWA, strip arms, de-dup, ..., profit


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

The output is in SAM format to "output-dir" with extra tags as described
in the @CO headers.

Thanks to Evan Boyle and Jay Hesselberth for (repeated) explanations.
Any mistakes are my own.

bpederse@gmail.com

LICENSE is MIT
