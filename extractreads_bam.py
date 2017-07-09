import sys
import os
import argparse
import subprocess
from timeit import default_timer as timer

parser = argparse.ArgumentParser(description='Gather potential reads for a region')
parser.add_argument('-r', help='Path of the reference fasta file', type=str)
parser.add_argument('--bed', help='Regions where the reads are going to be extracted from', type=str)
parser.add_argument('--bam', help='Bam file path', type=str)
parser.add_argument('-o', help='Path for working folder', default=os.getcwd(), type=str)
parser.add_argument('--overhang', help='overhang at left and right flanks', default=1000, type=int)

args = parser.parse_args()

def main():
    bed = open(args.bed, 'r')
    header = bed.readline()
    line = bed.readline()
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
            os.mkdir(os.path.join(args.o, array[4], array[5], title))
        except:
            sys.stderr.write("title folder exists\n")
        if int(array[1]) > args.overhang:
            start_pos = int(array[1]) - args.overhang
        else:
            start_pos = 1
        reg = array[0] + ":" + str(start_pos) + "-" + str(int(array[2]) + args.overhang)
        bam2fqPath = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), 'bam2fastq.py')
        sys.stderr.write("Extracting reads for " + title + "\n")
        start_time = timer()
        ps = subprocess.Popen(('samtools', 'view', args.bam, reg), stdout=subprocess.PIPE)
        output = subprocess.check_output(('python', bam2fqPath, os.path.join(args.o, array[4], array[5], title, "read1.fq.gz"),
                                          os.path.join(args.o, array[4], array[5], title, "read2.fq.gz"),
                                          os.path.join(args.o, array[4], array[5], title, "unpaired.fq.gz")), stdin=ps.stdout)
        ps.wait()
        end_time = timer()
        sys.stderr.write(str(end_time - start_time) + " seconds elapsed\n")
        line = bed.readline()


if __name__ == "__main__":
	main()