set -exo pipefail

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

if ! stat -t picard-tools-*
then
    curl -JLO "http://sourceforge.net/projects/picard/files/latest/download?source=file"
    unzip picard-tools-*.zip
fi
rm -f picard-tools-*.zip

R1=sample-202-20_S30_L001_R1_001.fastq.gz
R2=sample-202-20_S30_L001_R2_001.fastq.gz
#R1=/proj/Schwartz/brentp/2013/muc5b-resequencing-pilot/MIPs/NJ249-301-66_S66_L001_R1_001.fastq.gz
#R2=/proj/Schwartz/brentp/2013/muc5b-resequencing-pilot/MIPs/NJ249-301-66_S66_L001_R2_001.fastq.gz
#R1=a_R1.fq
#R2=a_R2.fq

mkdir -p results
rm -f results/*
# now map the reads.
python3 ../bwamips.py ref/$REF mips-design.txt $R1 $R2 --threads 16 \
    --picard-dir picard-tools-* \
    > results/sample.bam #| samtools view -bS - | samtools sort - results/sample

#samtools flagstat results/sample.bam
