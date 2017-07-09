import sys
import subprocess
import argparse
import os
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='blast contigs from each assembly')
parser.add_argument('--bed', help='bed file of SV regions', type=str)
parser.add_argument('-r', help='Path of the reference fasta file(for NCBI blast)', type=str)
parser.add_argument('-t', help='# of threads', type=int)
parser.add_argument('-o', help='Path for working folder', default=os.getcwd(), type=str)

args = parser.parse_args()

def blast(line):
    array = line.rstrip('\n').split('\t')
    title = '_'.join(array)
    folder = os.path.join(args.o, array[4], array[5], title)
    if os.path.exists(os.path.join(folder, "dipspades.fasta")):
        ps = subprocess.Popen(('blastn', '-db', args.r, '-query', os.path.join(folder, "dipspades.fasta"), '-outfmt', '5', '-out', os.path.join(folder, "dipspades.blastout"),
             '-max_hsps', '10', '-max_target_seqs', '10', '-reward', '1', '-penalty', '-5'))
        result = ps.wait()
        if result == 0:
            sys.stderr.write("successfuly blasted dipspades for " + title + "\n")
        else:
            sys.stderr.write("blast failed for " + title + "\n")
    if os.path.exists(os.path.join(folder, "spades.fasta")):
        ps = subprocess.Popen(('blastn', '-db', args.r, '-query', os.path.join(folder, "spades.fasta"), '-outfmt', '5', '-out', os.path.join(folder, "spades.blastout"),
             '-max_hsps', '10', '-max_target_seqs', '10', '-reward', '1', '-penalty', '-5'))
        result = ps.wait()
        if result == 0:
            sys.stderr.write("successfuly blasted spades for " + title + "\n")
        else:
            sys.stderr.write("blast failed for " + title + "\n")
    if os.path.exists(os.path.join(folder, "velvet.fasta")):
        ps = subprocess.Popen(('blastn', '-db', args.r, '-query', os.path.join(folder, "velvet.fasta"), '-outfmt', '5', '-out', os.path.join(folder, "velvet.blastout"),
             '-max_hsps', '10', '-max_target_seqs', '10', '-reward', '1', '-penalty', '-5'))
        result = ps.wait()
        if result == 0:
            sys.stderr.write("successfuly blasted velvet for " + title + "\n")
        else:
            sys.stderr.write("blast failed for " + title + "\n")
    if os.path.exists(os.path.join(folder, "abyss.fasta")):
        ps = subprocess.Popen(('blastn', '-db', args.r, '-query', os.path.join(folder, "abyss.fasta"), '-outfmt', '5', '-out', os.path.join(folder, "abyss.blastout"),
             '-max_hsps', '10', '-max_target_seqs', '10', '-reward', '1', '-penalty', '-5'))
        result = ps.wait()
        if result == 0:
            sys.stderr.write("successfuly blasted abyss for " + title + "\n")
        else:
            sys.stderr.write("blast failed for " + title + "\n")

def main():
    bed = open(args.bed, 'r')
    header = bed.readline()
    lines = []
    line = bed.readline()
    while len(line) > 3:
        lines.append(line)
        line = bed.readline()
    pool = Pool(args.t)
    pool.map(blast,lines)

if __name__ == '__main__':
    main()