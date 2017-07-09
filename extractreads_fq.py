import sys
import os
import gzip
import argparse
import commands
from timeit import default_timer as timer

parser = argparse.ArgumentParser(description='Gather potential reads for a region')
parser.add_argument('-r', help='Path of the reference fasta file', type=str)
parser.add_argument('--bed', help='Regions to extract reads for', type=str)
parser.add_argument('--r1', help='Path of the gzipped fastq file(s) of first reads. Comma separated.', type=str)
parser.add_argument('--r2', help='Path of the gzipped fastq file(s) of second reads. Comma separated.', type=str)
parser.add_argument('-o', help='Path for working folder', default=os.getcwd(), type=str)
parser.add_argument('--uniq', help='Max. number of allowed alternative mapping locations for seeds. Requires BWA.', type=int, default=-1)
parser.add_argument('--ss', help='Seed (k-mer) size', default=26, type=int)
parser.add_argument('--overhang', help='overhang at left and right flanks', default=1000, type=int)

args = parser.parse_args()


def reverse_comp(s):
    r = ''
    l = len(s)
    for i in range(l - 1, -1, -1):
        if s[i] == 'A' or s[i] == 'a':
            r += 'T'
        elif s[i] == 'T' or s[i] == 't':
            r += 'A'
        elif s[i] == 'G' or s[i] == 'g':
            r += 'C'
        elif s[i] == 'C' or s[i] == 'c':
            r += 'G'
        else:
            r += 'N'
    return r


def eliminate_repeats(seeds):
    fname = "svmap_seeds.fq.gz"
    fname = os.path.join(args.o, fname)
    myfq = gzip.open(fname, 'wb')
    qs = 'H' * args.ss
    for i in seeds:
        title = "@seed_" + seeds[i]
        read = title + '\n' + i + '\n+\n' + qs + '\n'
        myfq.write(read)
    sys.stderr.write(str(len(seeds)) + " seeds are in the list\n")
    myfq.close()
    cmd = "bwa aln -t 4 " + args.r + " " + fname + " > " + fname + ".sai"
    result = commands.getstatusoutput(cmd)
    cmd = "bwa samse -n " + str(args.uniq) + " " + args.r + " " + fname + ".sai " + fname + " > " + fname + ".sam"
    result = commands.getstatusoutput(cmd)
    selectedseeds = {}
    if result[0] == 0:
        sys.stderr.write("BWA run successful\n")
        mysam = open((fname + ".sam"), 'r')
        line = mysam.readline()
        while line:
            if line[0] == "@":
                line = mysam.readline()
                continue
            is_repeat = False
            is_acceptable = False
            is_unique = False
            line = line.split('\t')
            for c in line:
                if c == 'XT:A:R':
                    is_repeat = True
                if c == 'XT:A:U':
                    is_unique = True
                if len(c) > 4 and c[:4] == 'XA:Z':
                    is_acceptable = True
            if is_repeat and is_acceptable:
                selectedseeds[line[9]] = seeds[line[9]]
                selectedseeds[reverse_comp(line[9])] = line[0][5:]
            elif is_unique:
                selectedseeds[line[9]] = seeds[line[9]]
                selectedseeds[reverse_comp(line[9])] = line[0][5:]
            line = mysam.readline()
        mysam.close()
        return selectedseeds
    else:
        sys.stderr.write("Failed running BWA. Continuing with all of the extracted seeds\n")
        return seeds

def find_reads(seeds,files):
    found = 0
    count = 0
    start_time = timer()
    r1list = args.r1.split(',')
    r2list = args.r2.split(',')
    for i in range(len(r1list)):
        myfq1 = gzip.open(r1list[i],'rb')
        myfq2 = gzip.open(r2list[i],'rb')
        l1 = myfq1.readline()
        while l1:
            t1 = l1
            r1 = myfq1.readline()
            l1 = myfq1.readline()
            q1 = myfq1.readline()
            t2 = myfq2.readline()
            r2 = myfq2.readline()
            l2 = myfq2.readline()
            q2 = myfq2.readline()
            tmpdict = {}
            count += 1
            if count % 100000 == 0:
                current_time = timer()
                sys.stderr.write(str(count) + " pairs processed in total " +
                                 str(found) + " matching pairs found " +
                                 str(current_time - start_time) + " seconds/100k pairs\n")
                start_time = current_time
            p = 0
            while p < len(r1) - args.ss:
                if r1[p:p+args.ss] in seeds and seeds[r1[p:p + args.ss]] not in tmpdict:
                    tmp1 = t1 + r1 + "+\n" + q1
                    tmp2 = t2 + r2 + "+\n" + q2
                    files[seeds[r1[p:p + args.ss]]][0].write(tmp1)
                    files[seeds[r1[p:p + args.ss]]][1].write(tmp2)
                    tmpdict[seeds[r1[p:p + args.ss]]] = True
                    found += 1
                p += 1
            p = 0
            while p < len(r2) - args.ss:
                if r2[p:p+args.ss] in seeds and seeds[r2[p:p + args.ss]] not in tmpdict:
                    tmp1 = t1 + r1 + "+\n" + q1
                    tmp2 = t2 + r2 + "+\n" + q2
                    files[seeds[r2[p:p + args.ss]]][0].write(tmp1)
                    files[seeds[r2[p:p + args.ss]]][1].write(tmp2)
                    tmpdict[seeds[r2[p:p + args.ss]]] = True
                    found += 1
                p += 1
            l1 = myfq1.readline()
        myfq1.close()
        myfq2.close()

def main():
    bed = open(args.bed, 'r')
    header = bed.readline()
    line = bed.readline()
    seeds = {}
    read_files = {}
    start_time = timer()
    while len(line) > 3:
        array = line.rstrip('\n').split('\t')
        title = '_'.join(array)
        try:
            os.mkdir(os.path.join(args.o, array[4]))
        except:
            sys.stderr.write("type folder exists\n")
        try:
            os.mkdir(os.path.join(args.o, array[4], array[5]))
        except:
            sys.stderr.write("gt folder exists\n")
        try:
            os.mkdir(os.path.join(args.o, array[4], array[5],title))
        except:
            sys.stderr.write("title folder exists\n")
        read_files[title] = (gzip.open(os.path.join(args.o, array[4], array[5], title, "read1.fq.gz"),'wb'),gzip.open(os.path.join(args.o, array[4], array[5], title, "read2.fq.gz"),'wb'))
        sys.stderr.write(title + "\n")
        if int(array[1]) > args.overhang:
            start = int(array[1]) - args.overhang
        else:
            start = 1
        reg = array[0] + ":" + str(start) + "-" + str(int(array[2]) + args.overhang)
        cmd = "samtools faidx " + args.r + " " + reg
        result = commands.getstatusoutput(cmd)
        if result[0] == 0:
            seq = result[1].split('\n')
            seq = ''.join(seq[1:])
            for i in range(len(seq) - args.ss):
                a = seq[i:i + args.ss]
                seeds[a] = title
                seeds[reverse_comp(a)] = title
        else:
            sys.stderr.write("Error: Could not extract " + reg + " sequence with samtools")
            continue
        line = bed.readline()
    end_time = timer()
    print len(seeds)," seeds extracted ", end_time - start_time, " seconds elapsed"
    if args.uniq != -1:
        seeds = eliminate_repeats(seeds)
    sys.stderr.write(str(len(seeds)) + " seeds remained in the list\n")
    find_reads(seeds,read_files)



if __name__ == '__main__':
    main()
