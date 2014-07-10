set -ex

REF=chr6.fa

########################
# preparation
########################

# download chr6 since we know that's where our MIPs are designed for
if [[ ! -e ref/$REF ]]; then
    mkdir -p ref/
    cd ref
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr6.fa.gz
    gunzip ${REF}.gz
    bwa index $REF
    samtools faidx $REF
    cd ..
fi

R1=sample-202-20_S30_L001_R1_001.fastq.gz
R2=sample-202-20_S30_L001_R2_001.fastq.gz

rm -f results/*
# now map the reads.
python3 ../bwamips.py ref/$REF mips-design.txt $R1 $R2 --threads 16 \
    | samtools view -bS - | samtools sort - results/sample

samtools flagstat results/sample.bam
