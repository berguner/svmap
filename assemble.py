import sys
import shutil
import subprocess
import argparse
import os
import commands
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='Assemble reads extracted from each region')
parser.add_argument('--bed', help='.bed file of SV regions', type=str)
parser.add_argument('--dipspades', help='Path of dipspades.py', type=str)
parser.add_argument('--velvetg', help='Path of velvetg executable', type=str)
parser.add_argument('--velveth', help='Path of velveth executable', type=str)
parser.add_argument('--abyss', help='Path of abyss-pe executable', type=str, default=1)
parser.add_argument('-t', help='# of threads', type=int)
parser.add_argument('-o', help='Path for working folder', default=os.getcwd(), type=str)

args = parser.parse_args()

def assemble(line):
    array = line.rstrip('\n').split('\t')
    title = '_'.join(array)
    try:
        folder = os.path.join(args.o, array[4], array[5], title)
        sys.stderr.write(folder + "\n")
        if os.path.exists(os.path.join(folder, "unpaired.fq.gz")) and args.dipspades:
            sys.stderr.write("Running dipspades\n")
            ps = subprocess.Popen(('python', args.dipspades, '-1', os.path.join(folder, 'read1.fq.gz'), '-2',
                                   os.path.join(folder, 'read2.fq.gz'),
                                   '-s', os.path.join(folder, 'unpaired.fq.gz'), '--expect-rearrangements', '--cov-cutoff', '3',
                                   '-o', os.path.join(folder, 'dipspades')),
                                  stdout=open(os.path.join(folder, 'dipspades.log'), 'w'),
                                  stderr=open(os.path.join(folder, 'dipspades.err'), 'w'))
            result = ps.wait()
            if result != 0: sys.stderr.write("Couldn't run dipspades\n")
            else:
                os.rename(os.path.join(folder, 'dipspades/dipspades/consensus_contigs.fasta'),os.path.join(folder, 'dipspades.fasta'))
                os.rename(os.path.join(folder, 'dipspades/spades/scaffolds.fasta'),os.path.join(folder, 'spades.fasta'))
                shutil.rmtree(os.path.join(folder, 'dipspades'))
        elif args.dipspades:
            sys.stderr.write("Running dipspades\n")
            ps = subprocess.Popen(('python', args.dipspades, '-1', os.path.join(folder, 'read1.fq.gz'), '-2',
                                   os.path.join(folder, 'read2.fq.gz'),
                                   '--expect-rearrangements', '--cov-cutoff', '3', '-o', os.path.join(folder, 'dipspades')),
                                  stdout=open(os.path.join(folder, 'dipspades.log'), 'w'),
                                  stderr=open(os.path.join(folder, 'dipspades.err'), 'w'))
            result = ps.wait()
            if result != 0: sys.stderr.write("Couldn't run dipspades\n")
            else:
                os.rename(os.path.join(folder, 'dipspades/dipspades/consensus_contigs.fasta'),os.path.join(folder, 'dipspades.fasta'))
                os.rename(os.path.join(folder, 'dipspades/spades/scaffolds.fasta'),os.path.join(folder, 'spades.fasta'))
                shutil.rmtree(os.path.join(folder, 'dipspades'))
        if args.velvetg and args.velveth:
            sys.stderr.write("Running velveth\n")
            ps = subprocess.Popen(
                (args.velveth, os.path.join(folder, 'velvet'), '31', '-shortPaired', '-separate', '-fastq.gz',
                 os.path.join(folder, 'read1.fq.gz'), os.path.join(folder, 'read2.fq.gz')),
                stdout=open(os.path.join(folder, 'velveth.log'), 'w'),
                stderr=open(os.path.join(folder, 'velveth.err'), 'w'))
            result = ps.wait()
            if result != 0: sys.stderr.write("Couldn't run velveth\n")
            sys.stderr.write("Running velvetg\n")
            ps = subprocess.Popen(
                (args.velvetg, os.path.join(folder, 'velvet'), '-scaffolding', 'yes', '-very_clean', 'yes'),
                stdout=open(os.path.join(folder, 'velveth.log'), 'w'),
                stderr=open(os.path.join(folder, 'velveth.err'), 'w'))
            result = ps.wait()
            if result != 0: sys.stderr.write("Couldn't run velvetg\n")
            else:
                os.rename(os.path.join(folder, 'velvet/contigs.fa'), os.path.join(folder, 'velvet.fasta'))
        if os.path.exists(os.path.join(folder, "unpaired.fq.gz")) and args.abyss:
            if os.path.exists(os.path.join(folder, 'abyss')):
                shutil.rmtree(os.path.join(folder, 'abyss'))
                os.mkdir(os.path.join(folder, 'abyss'))
            else:
                os.mkdir(os.path.join(folder, 'abyss'))
            sys.stderr.write("Running abyss\n")
            cmd = args.abyss + ' name=' + title + ' k=64 in=\'../read1.fq.gz ../read2.fq.gz\' se=\'../unpaired.fq.gz\' --directory=' + os.path.join(
                folder, 'abyss')
            logfile = open(os.path.join(folder, 'abyss.log'), 'w')
            result = commands.getstatusoutput(cmd)
            logfile.write(result[1])
            if result[0] != 0: sys.stderr.write("Couldn't run abyss\n")
            else:
                shutil.copy2(os.path.join(folder, 'abyss/' + title + '-scaffolds.fa'), os.path.join(folder, 'abyss.fasta'))
                shutil.rmtree(os.path.join(folder, 'abyss'))
        elif args.abyss:
            if os.path.exists(os.path.join(folder, 'abyss')):
                shutil.rmtree(os.path.join(folder, 'abyss'))
                os.mkdir(os.path.join(folder, 'abyss'))
            else:
                os.mkdir(os.path.join(folder, 'abyss'))
            sys.stderr.write("Running abyss\n")
            cmd = args.abyss + ' name=' + title + ' k=64 in=\'../read1.fq.gz ../read2.fq.gz\' --directory=' + os.path.join(
                folder, 'abyss')
            logfile = open(os.path.join(folder, 'abyss.log'), 'w')
            result = commands.getstatusoutput(cmd)
            logfile.write(result[1])
            if result[0] != 0: sys.stderr.write("Couldn't run abyss\n")
            else:
                shutil.copy2(os.path.join(folder, 'abyss/' + title + '-scaffolds.fa'), os.path.join(folder, 'abyss.fasta'))
                shutil.rmtree(os.path.join(folder, 'abyss'))
    except:
        sys.stderr.write("failed this region " + title +"\n")
def main():
    bed = open(args.bed, 'r')
    header = bed.readline()
    lines = []
    line = bed.readline()
    while len(line) > 3:
        lines.append(line)
        #assemble(line)
        line = bed.readline()
    pool = Pool(args.t)
    pool.map(assemble,lines)

if __name__ == '__main__':
    main()