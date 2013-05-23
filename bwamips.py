"""
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

Thanks to Evan Boyle and Jay HesselBerth for (repeated) explainations
Any mistakes are my own.

bpederse@gmail.com

LICENSE is BSD
"""

import sys
import os
import os.path as op
import string
from itertools import izip, islice, groupby, takewhile, chain
from operator import itemgetter, attrgetter
from collections import Counter
import subprocess
from math import copysign

import argparse
from toolshed import reader, nopen
import multiprocessing


NUM_CORES=multiprocessing.cpu_count()
BARCODE_LENGTH = 5
MAX_READ_LENGTH = 152
VERSION = 0.1

class Bam(object):
    __slots__ = 'read flag chrom pos mapq cigar chrom_mate pos_mate tlen \
            seq qual other'.split()
    def __init__(self, args):
        for a, v in zip(self.__slots__[:11], args):
            setattr(self, a, v)
        self.other = args[11:]
        self.flag = int(self.flag)
        self.pos = int(self.pos)
        self.tlen = int(float(self.tlen))
        try:
            self.mapq = int(self.mapq)
        except ValueError:
            pass

    def __str__(self):
        return "\t".join(str(getattr(self, s)) for s in self.__slots__[:11]) \
                + "\t" + "\t".join(self.other)

    def is_first_read(self):
        return bool(self.flag & 0x40)

    def is_second_read(self):
        return bool(self.flag & 0x80)

    def is_plus_read(self):
        return not (self.flag & 0x10)

    def is_minus_read(self):
        return bool(self.flag & 0x10)

    def is_mapped(self):
        return not (self.flag & 0x4)

    def cigs(self):
       if self.cigar == "*":
           yield (0, None)
           raise StopIteration
       cig_iter = groupby(self.cigar, lambda c: c.isdigit())
       for g, n in cig_iter:
           yield int("".join(n)), "".join(cig_iter.next()[1])

    @classmethod
    def cig_len(self, cigs):
        return sum(c[0] for c in cigs if c[1] in ("M", "D", "N", "EQ", "X", "P"))

    def right_end(self):
        # http://www.biostars.org/p/1680/#1682
        length = sum(c[0] for c in self.cigs() if c[1] in ("M", "D", "N", "EQ", "X", "P"))
        return self.start + length - 1

    @property
    def start(self):
        return self.pos

    @property
    def len(self):
        return abs(self.tlen)


    def trim(self, left_trim, right_trim):
        if "M" == "".join(c for c in self.cigar if not c.isdigit()):
            #assert len(self.seq) + 1 == self.tlen, (len(self.seq), self.tlen,
            #        self.seq, self.cigar)
            self.seq = self.seq[left_trim:-right_trim]
            self.qual = self.qual[left_trim:-right_trim]
            self.tlen = copysign(len(self.seq), self.tlen)
            self.cigar = "%iM" % self.len
            self.pos += left_trim

        else:
            cigs = [list(x) for x in self.cigs()]
            if cigs[0][1] == "M" and cigs[0][0] >= left_trim and \
               cigs[-1][1] == "M" and cigs[-1][0] >= right_trim:
                self.seq = self.seq[left_trim:-right_trim]
                self.qual = self.qual[left_trim:-right_trim]
                tlen = self.len - left_trim - right_trim
                self.tlen = copysign(tlen, self.tlen)
                cigs[0][0] -= left_trim
                cigs[-1][0] -= right_trim
                self.cigar = "".join("%i%s" % tuple(t) for t in cigs)
                self.pos += left_trim
            elif sum(c[0] for c in cigs if c[1] == "S") \
                        > len(self.seq) - (left_trim + right_trim):
                    self.mapq = 0
                    self.flag |= 0x200
            else:
                self.trim_left(left_trim)
                self.trim_right(right_trim)

    def trim_left(self, n_bp):
        cigs = [list(x) for x in self.cigs()]
        n = 0

        new_cigs = cigs[:]
        for i, (clen, cop) in enumerate(cigs):
            if not cop in ("M", "I", "S", "EQ", "X"): continue
            n += clen
            if n >= n_bp:
                new_cigs[i] = (n - n_bp, cop)
                new_cigs = new_cigs[i:]
                self.seq = self.seq[clen - (n - n_bp):]
                self.qual = self.qual[clen - (n - n_bp):]
                break
            else:
                self.seq = self.seq[clen:]
                self.qual = self.qual[clen:]

        self.pos += n_bp
        self.cigar = "".join("%i%s" % tuple(t) for t in new_cigs)

    def trim_right(self, n_bp):
        cigs = [list(x) for x in self.cigs()][::-1]
        new_cigs = cigs[:]
        n = 0
        for i, (clen, cop) in enumerate(cigs):
            if not cop in ("M", "I", "S", "EQ", "X"): continue
            if n >= n_bp:
                new_cigs[i] = (n - n_bp, cop)
                new_cigs = new_cigs[i:]
                self.seq = self.seq[::-1][clen - (n - n_bp):][::-1]
                self.qual = self.qual[::-1][clen - (n - n_bp):][::-1]
                break
            else:
                self.seq = self.seq[::-1][clen:]
                self.qual = self.qual[::-1][clen:]

        self.cigar = "".join("%i%s" % tuple(t) for t in new_cigs[::-1])

def move_umi(fqs, output_dir):
    fq1, fq2 = fqs
    f1_out = nopen("{output_dir}/{base}.fastq.gz".format(output_dir=output_dir,
                base=op.basename(fq1).rsplit(".", 2)[0]), 'w')
    f2_out = nopen("{output_dir}/{base}.fastq.gz".format(output_dir=output_dir,
                base=op.basename(fq2).rsplit(".", 2)[0]), 'w')

    for r1, r2 in izip(fqiter(fq1), fqiter(fq2)):
        assert len(r1) == 4
        assert len(r2) == 4
        umi = r2[1][:5]
        r2[1] = r2[1][5:]
        r2[3] = r2[3][5:]

        r1[0] = r1[0].rsplit(" ", 1)[0].rsplit(":", 1)[0] + ":" + umi
        r2[0] = r2[0].rsplit(" ", 1)[0].rsplit(":", 1)[0] + ":" + umi
        assert not any("\n" in l for l in r1)
        assert not any("\n" in l for l in r2)
        print >>f1_out, "\n".join(r1)
        print >>f2_out, "\n".join(r2)

    f1_out.close(); f2_out.close()
    return f1_out.name, f2_out.name

def get_umi(read_name):
    return read_name.rsplit(":", 1)[1]

    cigs = [(int(x[0]), x[1]) for x in zip(diter, diter)]
    return cigs

def read_mips(mips_file):
    print >>sys.stderr, "reading", mips_file
    m = {'ext_probe_start':{}, 'lig_probe_start':{},
         'ext_probe_stop':{}, 'lig_probe_stop':{}}
    ss = m.keys()
    for d in reader(mips_file):
        for key in ss:
            d[key] = int(d[key])
        m['ext_probe_start'][(d['chr'], int(d['ext_probe_start']))] = d
        m['lig_probe_start'][(d['chr'], int(d['lig_probe_start']))] = d
        m['lig_probe_stop'][(d['chr'], int(d['lig_probe_stop']))] = d
        m['ext_probe_stop'][(d['chr'], int(d['ext_probe_stop']))] = d
    return m

def dearm_bam(bam, mips_file):
    """
    targetting arms:
     lig is 5' of mip
     ext is 3' of mip
               ___-------___
        ___----             -----___
       |                            |
      |                              |
     |  ext                     lig   |
     ------                    --------
     e    e                    l      l
     s    s                    s      s
     t    t                    t      t
     a    o                    a      o
     r    p                    r      p
     t                         t

    for + strand reads: expected start of capture is ext_probe_start
                      : expected stop  of capture is lig_probe_stop

    for - strand reads: expected start of capture is lig_probe_start!!
                      : expected stop  of capture is ext_probe_stop!!

    """
    mips = read_mips(mips_file)

    bam_iter = reader("|samtools view -H {bam}".format(bam=bam), header=False)
    for toks in bam_iter:
        yield "\t".join(toks)

    bam_iter = reader("|samtools view {bam} -L region.bed".format(bam=bam), header=Bam)
    n = 0
    counts = [0, 0, 0, 0]
    for aln in bam_iter:

        oseq, oqual, ocigar, opos = aln.seq, aln.qual, aln.cigar, aln.pos

        if not aln.is_mapped():
            yield str(aln)
            continue
        assert aln.is_first_read() != aln.is_second_read()
        assert aln.is_plus_read() != aln.is_minus_read()
        first = right = second = False

        if aln.is_first_read() and aln.is_minus_read():
            lookup_pair = aln.chrom, aln.right_end()
            lookup_field = "lig_probe_stop"
            idx = 0

        elif aln.is_second_read() and aln.is_plus_read():
            lookup_pair = aln.chrom, aln.pos
            lookup_field = "ext_probe_start"
            idx = 1

        elif aln.is_first_read() and aln.is_plus_read():
            lookup_pair = aln.chrom, aln.pos
            lookup_field = "lig_probe_start"
            idx = 2

        if aln.is_second_read() and aln.is_minus_read():
            lookup_field = "ext_probe_stop"
            lookup_pair = aln.chrom, aln.right_end()
            idx = 3

        aln.other.append("XO:Z:%s" % oseq)
        # allow some wiggle room
        for offset in [0, -1, 1]:
            mips[lookup_field]
            try:
                pair = lookup_pair[0], lookup_pair[1] + offset
                mip = mips[lookup_field][pair]
                n += 1
                counts[idx] += 1
                #print mip['probe_strand'], ['+', '-'][int(aln.is_minus_read())], aln.is_first_read()
                break
            except KeyError:
                aln.flag |= 0x200 # not passing QC
                yield str(aln)
                continue
        else:
            aln.flag |= 0x200 # not passing QC
            yield str(aln)
            continue

        aln.other.append("BC:Z:%s" % aln.read.rsplit(":", 1)[1])
        aln.other.append("OC:Z:%s" % aln.cigar)
        aln.other.append("OP:i:%i" % aln.pos)

        aln.other.append("XI:i:%s" % mip['>index'])
        lig_len = mip['lig_probe_stop'] - mip['lig_probe_start']
        ext_len = mip['ext_probe_stop'] - mip['ext_probe_start']

        if idx in (2, 3):
            if mip['lig_probe_start'] < mip['lig_probe_stop'] < mip['ext_probe_start'] < mip['ext_probe_stop']:
                aln.trim(lig_len, ext_len)
            else:
                aln.flag |= 0x200 # not passing QC
                pass
        elif idx in (0, 1):
            if mip['ext_probe_start'] < mip['ext_probe_stop'] < mip['lig_probe_start'] < mip['lig_probe_stop']:
                aln.trim(ext_len, lig_len)
            else:
                aln.flag |= 0x200 # not passing QC
                pass

        if len(aln.seq) < 20: # weird cigar like 120S28M
            aln.seq, aln.qual, aln.cigar, aln.pos = oseq, oqual, ocigar, opos
            aln.flag |= 0x200 # not passing QC
        yield str(aln)

    print >>sys.stderr, "found %i MIPs" % n
    print >>sys.stderr, "counts", counts

from string import maketrans
def revcomp(seq, tbl=maketrans('ACGT', 'TGCA')):
    return seq.translate(tbl)


def bwamips(fastqs, ref_fasta, mips, output_dir, num_cores=NUM_CORES,
        umi_fn=move_umi):

    def base(f):
        return op.basename(f).replace('.gz', '').replace('.fastq',
                '').replace('.fq', '')
    #"""
    fq1, fq2 = move_umi(fastqs, output_dir)

    name = "".join([x[0] for x in takewhile(lambda (a, b): a == b, zip(base(fq1), base(fq2)))])
    if name.endswith("_R"): name = name[:-2]

    bam = bwa_mem((fq1, fq2), name, ref_fasta, output_dir, num_cores)
    """
    name = "test"
    bam = "tmp/t.bam"
    """
    sam_out = open('{output_dir}/{name}.sam'.format(**locals()), 'w')
    bam3 = dedup_sam(dearm_bam(bam, mips), get_umi, sam_out, mips)
    sam_out.close()

def dedup_sam(sam_iter, umi_fn, out=sys.stdout, mips_file=''):
    # first print the header
    line = None
    for line in sam_iter:
        if not line.startswith("@"):
            break
        print >>out, line

    args = " ".join(sys.argv)
    print >>out, '@PG\tID:bwamips\tPN:bwamips.py\tCL:%s\tVN:%s' \
                % (args, VERSION)
    print >>out, '@CO\tXI:i tag indicates the >index of the mip from %s' % mips_file
    print >>out, '@CO\tXO:Z tag indicates the original, mapped sequence'

    # group to reads at same position.
    counts = Counter()
    # add back the last line
    #sam_iter_ = chain([Bam(line.split("\t"))], (Bam(x.strip().split("\t")) for x in sam_iter))
    sam_iter_ = (Bam(x.strip().split("\t")) for x in sam_iter)
    for cpos, reads in groupby(sam_iter_, lambda r: (r.chrom, r.pos)):
        ureads = sorted((umi_fn(r.read), r) for r in reads)

        # group to reads with the same umi at that position
        for umi, umi_group in groupby(ureads, key=itemgetter(0)):
            rgroup = [r[1] for r in umi_group]
            # print the read with the highest quality
            counts.update([len(rgroup)])
            for i, aln in enumerate(sorted(rgroup,
                key=attrgetter('mapq'), reverse=True)):
                if i > 0:
                    aln.flag |= 0x400 # PCR or optical duplicate
                    pass
                print >>out, str(aln)
    print >>sys.stderr, counts.most_common(20)

def bwa_mem(fastqs, name, ref_fasta, output_dir, num_cores):
    fq1, fq2 = fastqs
    cmd = ("set -o pipefail && " # stolen from bcbio
           "bwa mem -M -t {num_cores} -R {rg} -v 1 {ref_fasta} "
           "{fq1} {fq2} "
           "| samtools view -b -S -u - "
           "| samtools sort -m 2G "
           "- {output_dir}/{name}")

    rg = "'@RG\\tID:%s\\tSM:%s'" % (name, name)
    bam = "{output_dir}/{name}.bam".format(**locals())

    try:
        subprocess.check_call(cmd.format(**locals()), shell=True)
    except Exception, e:
        try:
            os.unlink(bam)
        except OSError:
            pass
        raise e
    return bam

def fqiter(fq, n=4):
    with nopen(fq) as fh:
        while True:
            rec = [x.strip('\r\n') for x in islice(fh, n)]
            if not rec: raise StopIteration
            assert all(rec) and len(rec) == 4
            yield rec
def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--output-dir', help='output directory')
    p.add_argument('ref_fasta', help='reference fasta already index by bwa 0.7.4+')
    p.add_argument('mips', help='mips file')
    p.add_argument('fastqs', nargs=2, help='fastqs', metavar=('R1', 'R2'))
    args = p.parse_args()

    for f in (args.ref_fasta, args.mips, args.fastqs[1], args.fastqs[1],
              args.ref_fasta + ".sa"):
        if not op.exists(f):
            print >>sys.stderr, "%s missing" % f
            sys.exit(not p.print_help())

    args.output_dir = op.abspath(op.expanduser(args.output_dir))
    if not op.exists(args.output_dir):
        os.makedirs(args.output_dir)

    bwamips(args.fastqs, args.ref_fasta, args.mips, args.output_dir)

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
