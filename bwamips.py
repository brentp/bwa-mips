"""
With molecular inversion probes, we map reads to the genome that include the
ligation and extension arms, along with a molecular tag (AKA UMI).
This script takes:
    1) ref.fasta
    2) mips design file (likely from MIPgen)
    3) de-multiplxed, paired-end fastqs

and moves the UMI into the read-name, aligns the reads using bwa mem,
strips that ligation and extension arms from the alignment, adjusts the
position, cigar, sequence, and quality accounting (more or less) for
insertions, deletions, and masked sequence.

The output is in BAM format to stdout with extra tags as described
in the @CO headers.

Thanks to Evan Boyle and Jay Hesselberth for (repeated) explanations
Any mistakes are my own.

bpederse@gmail.com

LICENSE is MIT
"""
from __future__ import print_function, division

REPORT = """
================
Alignment Report
================

Bases in genome: {genome_size:,}
Bases in target region: {target_size:,} ({target_ratio_pct:.4f}% of genome)


MIPs found: {mip:,}
Off-target reads: {off-target:,}
Unmapped reads: {unmapped:,}

Observed / expected enrichment where expected is based on size
of target region relative to size of genome.

FOLD ENRICHMENT
===============
low - high: {lo_enrich:.2f} - {hi_enrich:.2f}

Low estimate uses unmapped reads as well as off-target.

% READS ON TARGET
=================
{on_target:.2f}%

"""

import sys
import os
import tempfile
import atexit
import os.path as op
import errno

from itertools import islice, groupby, takewhile, chain
try:
    from itertools import izip
except ImportError: #py3
    izip = zip

from operator import attrgetter
from collections import Counter
import subprocess
from math import copysign
from subprocess import Popen, PIPE

import argparse
import toolshed as ts

__version__ = "0.1.5"

MAX_READ_LENGTH = 152

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

    def __repr__(self):
        return "Bam({chr}:{start}:{read}".format(chr=self.chrom,
                start=self.pos, read=self.read)

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
           yield int("".join(n)), "".join(next(cig_iter)[1])

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
        if self.cigar.startswith("0M"):
            self.cigar = self.cigar[2:]
        # 0M, but not 10M
        if self.cigar.endswith("0M") and not self.cigar[-3].isdigit():
            self.cigar = self.cigar[:-2]

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

def get_base_name(fq1, fq2, _again=True):
    def base(f):
        return op.basename(f).replace('.gz', '')\
                .replace('.fastq', '').replace('.fq', '')

    name = "".join([x[0] for x in takewhile(lambda ab: ab[0] == ab[1], \
                                           zip(base(fq1), base(fq2)))])
    name = name.rstrip("_")
    if name.endswith("_R"): name = name[:-2]
    if len(name) < 3 and _again:
        sys.stderr.write("flipping: %s\n" % name)
        return get_base_name(base(fq1)[::-1], base(fq2)[::-1], False)
    name = (name if _again else name[::-1]).lstrip("._-")
    return name

def move_tag(fq1, fq2, umi_len):

    for r1, r2 in izip(fqiter(fq1), fqiter(fq2)):
        assert len(r1) == 4
        assert len(r2) == 4
        umi = r2[1][:umi_len]
        r2[1] = r2[1][umi_len:]
        r2[3] = r2[3][umi_len:]

        # use bwa's comment stuff to send this to the alignment.
        r1[0] = r1[0].split(" ", 1)[0] + " BC:Z:" + umi
        r2[0] = r2[0].split(" ", 1)[0] + " BC:Z:" + umi
        assert not any("\n" in l for l in r1)
        assert not any("\n" in l for l in r2)
        print("\n".join(r1))
        print("\n".join(r2))

def get_umi(read):
    return next(o for o in read.other if o.startswith("BC:Z:"))

def read_mips(mips_file):
    sys.stderr.write("reading %s\n" % mips_file)
    m = {'ext_probe_start':{}, 'lig_probe_start':{},
         'ext_probe_stop':{}, 'lig_probe_stop':{}}
    ss = m.keys()
    for d in ts.reader(mips_file):
        for key in ss:
            d[key] = int(d[key])
        m['ext_probe_start'][(d['chr'], int(d['ext_probe_start']))] = d
        m['lig_probe_start'][(d['chr'], int(d['lig_probe_start']))] = d
        m['lig_probe_stop'][(d['chr'], int(d['lig_probe_stop']))] = d
        m['ext_probe_stop'][(d['chr'], int(d['ext_probe_stop']))] = d
    return m

def target_size_from_mips(mips, pad=0):
    tmp = open(mktemp(suffix=".bed"), "w")
    for k in mips: # ext/lig_probe_start/stop
        for (chrom, pos), d in mips[k].items():
            posns = [int(d[p]) for p in "ext_probe_start ext_probe_stop lig_probe_start lig_probe_stop".split()]
            tmp.write("%s\t%i\t%i\n" % \
                       (chrom, max(0, min(posns) - pad), max(posns) + pad))
    tmp.close()
    size = 0
    for toks in ts.reader("|tail -n+2 %s | sort -k1,1 -k2,2n | bedtools merge -i stdin" % tmp.name,
            header=False):
        size += int(toks[2]) - int(toks[1])
    return size

def dearm_sam(sam_gz, mips_file):
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
    target_size = target_size_from_mips(mips)

    sam_iter = ts.reader(sam_gz, header=False)

    # calculcate genome size from bam header
    genome_size = 0
    for toks in sam_iter:
        if not toks[0].startswith("@"): break
        if toks[0].startswith("@SQ"):
            genome_size += next(int(x.split(":")[1]) for x in toks if x.startswith("LN:"))
        yield "\t".join(toks)

    target_ratio_pct = 100.0 * target_size / genome_size

    sam_iter2 = chain([Bam(toks)], (Bam(t) for t in sam_iter))
    n = 0
    counts = [0, 0, 0, 0]
    stats = {'unmapped': 0, 'mip': 0, 'off-target': 0}
    for aln in sam_iter2:
        if not aln.is_mapped():
            yield str(aln)
            stats['unmapped'] += 1
            continue

        oseq, oqual, ocigar, opos = aln.seq, aln.qual, aln.cigar, aln.pos
        assert aln.is_first_read() != aln.is_second_read()
        assert aln.is_plus_read() != aln.is_minus_read()

        idx = None
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

        elif aln.is_second_read() and aln.is_minus_read():
            lookup_pair = aln.chrom, aln.right_end()
            lookup_field = "ext_probe_stop"
            idx = 3

        # allow some wiggle room
        mip = None
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
                # check next position.
                continue
        else:
            aln.flag |= 0x200 # not passing QC
            stats['off-target'] += 1
            yield str(aln)
            continue

        aln.other.extend([
            "OP:i:%i" % aln.pos,
            "XI:i:%s" % mip.get('>index', mip.get('>mip_key')),
            "XO:Z:%s" % oseq,
            "OC:Z:%s" % aln.cigar,
            ])

        lig_len = mip['lig_probe_stop'] - mip['lig_probe_start']
        ext_len = mip['ext_probe_stop'] - mip['ext_probe_start']

        if idx in (2, 3):
            if mip['lig_probe_start'] < mip['lig_probe_stop'] < mip['ext_probe_start'] < mip['ext_probe_stop']:
                aln.trim(lig_len, ext_len)
            else:
                aln.flag |= 0x200 # not passing QC
        elif idx in (0, 1):
            if mip['ext_probe_start'] < mip['ext_probe_stop'] < mip['lig_probe_start'] < mip['lig_probe_stop']:
                aln.trim(ext_len, lig_len)
            else:
                aln.flag |= 0x200 # not passing QC

        if len(aln.seq) < 20: # weird cigar like 120S28M
            aln.seq, aln.qual, aln.cigar, aln.pos = oseq, oqual, ocigar, opos
            aln.flag |= 0x200 # not passing QC
        yield str(aln)
        stats['mip'] += 1

    hi_enrich = (stats['mip'] / stats['off-target']) / (target_size / genome_size)
    lo_enrich = (stats['mip'] / (stats['unmapped'] + stats['off-target'])) / (target_size / genome_size)

    on_target = 100.0 * (stats['mip'] / (stats['mip'] + stats['unmapped'] + stats['off-target']))

    info = locals()
    info.update(stats)
    print(REPORT.format(**info), file=sys.stderr)

def rm(f):
    try: os.unlink(f)
    except: pass

def mktemp(*args, **kwargs):
    f = tempfile.mktemp(*args, **kwargs)
    atexit.register(rm, f)
    return f

def bwamips(fastqs, ref_fasta, mips, num_cores, umi_length, picard):

    tmp_sam_name = mktemp(suffix=".sam.gz")
    name = get_base_name(*fastqs)
    sam_gz = bwa_mem(fastqs, name, ref_fasta, tmp_sam_name, num_cores, umi_length)
    sam_out = sys.stdout
    if op.exists("{picard}/FixMateInformation.jar".format(picard=picard)):
        jar = "{picard}/FixMateInformation.jar"
    else:
        jar = "{picard}/picard.jar FixMateInformation"

    out = Popen("java -jar -Xmx2G {jar} \
            SO=coordinate I=/dev/stdin O=/dev/stdout".format(jar=jar.format(picard=picard)),
            stderr=sys.stderr,
            stdout=sys.stdout, stdin=PIPE, shell=True)
    if sys.version_info[0] > 2:
        import io
        out.stdin = io.TextIOWrapper(out.stdin)

    sam_out = out.stdin
    #
    dedup_sam(dearm_sam(sam_gz, mips), get_umi, sam_out, mips)
    sam_out.close()
    out.wait()


def dedup_sam(sam_iter, get_umi_fn, out=sys.stdout, mips_file=''):
    line = None
    j = 0
    for line in sam_iter:
        if not line.startswith("@"):
            break
        out.write(line + "\n")

    args = " ".join(sys.argv)
    out.write('@PG\tID:bwamips\tPN:bwamips.py\tCL:%s\tVN:%s\n' \
                % (args, __version__))
    out.write('@CO\tXI:i tag indicates the >index or >mip_key of the mip from %s\n' %
            mips_file)
    out.write('@CO\tXO:Z tag indicates the original, mapped sequence\n')

    # group to reads at same position.
    counts = Counter()
    # add back the last line
    sorted_iter = sorted([Bam(x.strip().split("\t")) for x in sam_iter] +
                         [Bam(line.strip().split("\t"))],
                        key=lambda r: (r.chrom, r.pos, get_umi_fn(r)))
    q = 0
    for cpos, reads in groupby(sorted_iter, lambda r: (r.chrom, r.pos, get_umi_fn(r))):
        reads = list(reads)

        if cpos[0] == "*":
            q += len(reads)
            j += len(reads)
            for r in reads:
                out.write(str(r) + '\n')
            continue

        reads1 = [r for r in reads if r.is_first_read()]
        reads2 = [r for r in reads if r.is_second_read()]
        assert len(reads1) + len(reads2) == len(reads)

        # group to reads with the same umi at that position
        for readset in (reads1, reads2):
            """
            best, others = rgroup[0], rgroup[1:]
            # TODO adjust base-quality of best if the others do no match.
            if best.cigar != '*' and set([o.cigar for o in others if o.cigar != '*']) != set([best.cigar]):
                print cpos, best.cigar, [o.cigar for o in others], best.is_first_read(), best.is_plus_read()
            """
            if len(readset) == 0: continue
            # print the read with the highest quality
            counts.update([len(readset)])
            j += len(readset)
            for ir, aln in enumerate(sorted(readset, key=attrgetter('mapq'),
                reverse=True)):
                if ir > 0:
                    aln.flag |= 0x400 # PCR or optical duplicate
                out.write(str(aln) + "\n")
    out.flush()
    sys.stderr.write(str(counts.most_common(20)) + "\n")
    sys.stderr.write("wrote %i reads\n" % j)

def bwa_mem(fastqs, name, ref_fasta, tmp_sam_name, num_cores, umi_length):
    fq1, fq2 = fastqs

    # use bwa's streaming stuff and interleaved fq so we dont write temporary
    # fqs with barcode removed.
    fn = __file__
    fqbc = "'<python {fn} detag {fq1} {fq2} {umi_length}'".format(**locals())

    cmd = ("bwa mem -p -C -M -t {num_cores} -R {rg} -v 1 {ref_fasta} "
           "{fqbc} "
           "| grep -v '^@PG' "
           "| gzip -c - > {tmp_sam_name}")

    rg = "'@RG\\tID:%s\\tSM:%s\\tPL:illumina'" % (name, name)
    sys.stderr.write(cmd.format(**locals()) + "\n")

    try:
        subprocess.check_call(cmd.format(**locals()), shell=True)
    except Exception as e:
        rm(tmp_sam_name)
        raise e
    return tmp_sam_name

def fqiter(fq, n=4):
    with ts.nopen(fq) as fh:
        fqclean = (x.strip("\r\n") for x in fh if x.strip())
        while True:
            rec = [x for x in islice(fqclean, n)]
            if not rec: raise StopIteration
            assert all(rec) and len(rec) == 4
            yield rec

def main():

    # this subcommand is called internally by the main, user-visible script.
    if len(sys.argv) > 1 and sys.argv[1] == "detag":
        try:
            move_tag(sys.argv[2], sys.argv[3], int(sys.argv[4]))
        except IOError as e:
            if e.errno != errno.EPIPE:
                raise
        sys.exit(0)

    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--threads', help='number of threads for bwa mem',
            default=2, type=int)
    p.add_argument('--umi-length', help='length of umi', default=5, type=int)
    p.add_argument('--picard-dir', help='path to picard directory',
            required=True)
    p.add_argument('ref_fasta', help='reference fasta already index by bwa 0.7.4+')
    p.add_argument('mips', help='mips design file. e.g. from MIPgen')
    p.add_argument('fastqs', nargs=2, metavar=('FASTQ'))
    args = p.parse_args()

    for f in (args.ref_fasta, args.mips, args.fastqs[0], args.fastqs[1],
              args.ref_fasta + ".sa"):
        if not op.exists(f):
            sys.stderr.write("%s missing\n" % f)
            sys.exit(not p.print_help())

    bwamips(args.fastqs, args.ref_fasta, args.mips, args.threads,
            args.umi_length, args.picard_dir)

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
