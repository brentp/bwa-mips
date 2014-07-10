Example
=======
This example will run through using `bwamips.py` on the example
fastqs in this directory and chromosome 6 from UCSC.

TLDR
----

    python ../bwamips.py $REF mips-design.txt $R1 $R2 > res.sam

And see the `run.sh` in this directory

Data
====


We have an example mips design file with mips for chromosome 6.

Design File
-----------

```Shell
$ head mips-design.txt
>index	chr	ext_probe_start	ext_probe_stop	lig_probe_start	lig_probe_stop
1	chr6	7541629	7541648	7541756	7541780
2	chr6	7541900	7541921	7542029	7542051
3	chr6	7541535	7541556	7541664	7541686
4	chr6	7541930	7541950	7542058	7542081
5	chr6	7541800	7541820	7541928	7541951
6	chr6	7541676	7541697	7541805	7541827
7	chr6	7542029	7542050	7542158	7542180
8	chr6	7590782	7590801	7590909	7590933
9	chr6	7541750	7541771	7541879	7541901
```

Your design file *must* have at least those columns.

**Note** that you'll need to make sure that your design file and the
reference to which you align are from the same genome assembly!.


Prepare Reference Fasta
-----------------------

For this example, we will download only chr6 since we know that's
where our MIPs are designed for but in practice one should always
align to the full genome.

The following block downloads and indexes the fasta file

```Shell

mkdir -p ref/
cd ref
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr6.fa.gz
gunzip ${REF}.gz
bwa index $REF
samtools faidx $REF
cd ..
```

Run bwamips
-----------


```Shell

R1=sample-202-20_S30_L001_R1_001.fastq.gz
R2=sample-202-20_S30_L001_R2_001.fastq.gz

mkdir -p results/
# now map the reads. (note that python2.7+ and python3 are supported)
python3 ../bwamips.py ref/$REF mips-design.txt $R1 $R2 --threads 16 \
                   | samtools view -bS - | samtools sort - results/sample

samtools flagstat results/sample.bam

```
