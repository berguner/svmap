from Bio.Blast import NCBIXML as nx
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import argparse

parser = argparse.ArgumentParser(description='Parse & plot blast results')
parser.add_argument('--bed', help='bed file of SV regions', type=str)
parser.add_argument('--show_plot', help='Set to 1 if you want to see interactive plot after generation of eash plot', default=0, type=int)
parser.add_argument('-o', help='Path for working folder', default=os.getcwd(), type=str)
parser.add_argument('--overhang', help='overhang at left and right flanks', default=1000, type=int)

args = parser.parse_args()

minContigSize = 150
minCoverage = 3
minHspLength = 26
minBitScore = 1.7


def ParseQueryTitle(title,assembler_name):
    tmp = {}
    if assembler_name == 'spades' or assembler_name == 'velvet':
        title = title.split('_')
        tmp["contigName"] = title[1]
        tmp["coverage"] = float(title[5])
        tmp["length"] = int(title[3])
    elif assembler_name == 'dipspades':
        title = title.split('_')
        tmp["contigName"] = title[0]
        tmp["coverage"] = 1.0
        tmp["length"] = int(title[2])
    elif assembler_name == 'abyss':
        title = title.split(' ')
        tmp["contigName"] = title[0]
        tmp["coverage"] = int(title[2])/int(title[1])
        tmp["length"] = int(title[1])
    return tmp


def parseblastout(line,assembler_name):
    array = line.rstrip('\n').split('\t')
    title = '_'.join(array)
    folder = os.path.join(args.o, array[4], array[5], title)
    file = os.path.join(folder, assembler_name + ".blastout")
    if not os.path.exists(file):
        sys.stderr.write(file + " does not exist\n")
        return
    bfile = open(file,'r')
    blast_records = nx.parse(bfile)
    blast_records = list(blast_records)
    sys.stderr.write("Blast output file is parsed successfully!\nThere are " + str(len(blast_records)) + " Blast records in the file\n")

    region = (array[0],int(array[1])-args.overhang,int(array[2])+args.overhang)
    ymax = 0
    xLabels = []
    yLabels = []
    contigList = []
    colorCodes = "bgrcmyk"
    colorCount = 0
    for b in blast_records:
        colorCount += 1
        cn = colorCodes[colorCount % 7]
        if len(b.alignments) > 0 and b.query_length >= minContigSize:
            contigInfo = ParseQueryTitle(b.query,assembler_name)
            contigBits = [0.0] * b.query_length
            if b.query_length > ymax: ymax = b.query_length + 20
            contigHsps = []
            sys.stderr.write("There are " + str(len(b.alignments)) + " alignments for " + b.query + "\n")
            for a in b.alignments:
                chromosome = a.title.split(' ')[1]
                sys.stderr.write("There are " + str(len(a.hsps)) + " hsps for chromosome " + chromosome + "\n")
                for hsp in a.hsps:
                    tmpContigBits = contigBits[hsp.query_start-1:hsp.query_start+hsp.align_length-1]
                    avgTmpContigBits = (sum(tmpContigBits)/len(tmpContigBits))
                    hspBits = [hsp.bits / hsp.align_length] * hsp.align_length
                    avgHspBits = (sum(hspBits)/len(hspBits))
                    if avgHspBits > minBitScore and avgTmpContigBits < avgHspBits-0.1:
                        contigHsps.append((chromosome,avgHspBits,hsp))
                        for i in range(hsp.query_start-1,hsp.query_start+hsp.align_length-1):
                            if i >= 0 and i < len(contigBits):
                                contigBits[i] = avgHspBits
            contigScore = (sum(contigBits)/b.query_length) * contigInfo["coverage"]
            contigHsps = sorted(contigHsps, key=lambda x: x[2].query_start)
            contigList.append((contigHsps,contigScore,contigInfo))
            for index in range(len(contigHsps)):
                chromosome, avgHspBits, hsp = contigHsps[index]
                hspPlotLabel, ctrl = '', False
                if len(contigHsps) > 1: hspPlotLabel = contigInfo['contigName']+'_'+str(index); ctrl = True
                else: hspPlotLabel = contigInfo['contigName']
                # if hsp.bits / hsp.align_length > minBitScore and hsp.align_length >= minHspLength:
                if avgHspBits > minBitScore and chromosome == region[0] and hsp.sbjct_start >= region[1] and (hsp.sbjct_start + (hsp.align_length * hsp.frame[1])) <= region[2]:
                    sys.stderr.write("plotting " + str(hsp.query_start) +" "+ str(hsp.sbjct_start) +" "+ str(hsp.query_start + (hsp.align_length * hsp.frame[0]))+" "+str(hsp.sbjct_start + (hsp.align_length * hsp.frame[1]))+"\n")
                    plt.plot([hsp.sbjct_start,hsp.sbjct_start + (hsp.align_length * hsp.frame[1])],[hsp.query_start,hsp.query_start + (hsp.align_length * hsp.frame[0])],'k-',lw=3,c=cn)
                    if ctrl: xLabels.append(hsp.sbjct_start); xLabels.append(hsp.sbjct_start + (hsp.align_length * hsp.frame[1])); yLabels.append(hsp.query_start); yLabels.append(hsp.query_start + (hsp.align_length * hsp.frame[0]));
                    plt.annotate(hspPlotLabel,xy=((hsp.sbjct_start + (hsp.sbjct_start + (hsp.align_length * hsp.frame[1])))/2,(hsp.query_start+(hsp.query_start + (hsp.align_length * hsp.frame[0])))/2),
                                 xytext=((hsp.sbjct_start + (hsp.sbjct_start + (hsp.align_length * hsp.frame[1])))/2,(hsp.query_start+(hsp.query_start + (hsp.align_length * hsp.frame[0])))/2))
        else:
            sys.stderr.write("Length of " + b.query + " is less than MinContigSize: " + str(minContigSize) + "\t Skipping...\n")
    plt.xticks(rotation=45)
    plt.xlabel('chr' + region[0])
    plt.ylim(0,ymax+50)
    plt.xlim(region[1]-50,region[2]+50)
    ax = plt.gca()
    ax.tick_params(pad=25)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.set_xticks(xLabels,minor=False)
    ax.set_yticks(yLabels,minor=False)
    ax.yaxis.grid(True,which='major')
    ax.xaxis.grid(True,which='major')
    pngname = os.path.join(args.o, array[4], array[5], title + "_" + assembler_name + ".png")
    plt.savefig(pngname,dpi=300,bbox_inches='tight')
    if args.show_plot: plt.show()
    plt.close()

    contigList = sorted(contigList,reverse=True ,key=lambda x:x[1])
    outname = os.path.join(args.o, array[4], array[5],title + "_" + assembler_name + "_svmap.out")
    outfile = open(outname,'w')
    for c in contigList:
        outfile.write("##contig_ID:" + c[2]["contigName"] +
                      " length:"+ str(c[2]["length"]) +
                      " coverage:"+str(c[2]["coverage"]) +
                      " contig_score:"+str(c[1]) +
                      " number_of_hits:"+str(len(c[0])) + "\n")
        #if len(c[0]) > 1:
        previousHsp = None
        previousChr = None
        for i in range(len(c[0])):
            hsp = c[0][i][2]
            if previousHsp:
                prevQueryEnd = (previousHsp.query_start + (previousHsp.align_length * previousHsp.frame[0]) - 1)
                if prevQueryEnd > hsp.query_start: prevQueryEnd = hsp.query_start - 1       #left align breakpoints
                events = []
                if c[0][i][0] != previousChr:
                    events.append("inter_chromosomal_translocation")
                else:
                    if previousHsp.frame[1] != hsp.frame[1]:
                        events.append("inversion")
                    if hsp.query_start - prevQueryEnd > 1:
                        events.append("insertion:Query:" + str(prevQueryEnd) +"-"+ str(hsp.query_start))
                    if abs(hsp.sbjct_start - (previousHsp.sbjct_start + (previousHsp.align_length * previousHsp.frame[1]))):
                        events.append("deletion:" + c[0][i][0] + ":" +
                                      str(min(hsp.sbjct_start , (previousHsp.sbjct_start + (previousHsp.align_length * previousHsp.frame[1])))) + "-" +
                                      str(max(hsp.sbjct_start , (previousHsp.sbjct_start + (previousHsp.align_length * previousHsp.frame[1])))))
                outfile.write("\tEVENTS: " + ','.join(events) + "\n")
            strand = ''
            if hsp.frame[1] < 0: strand= 'Rev'
            else: strand = 'Fwd'
            gaps = 0
            mismatches = 0
            for n in hsp.query:
                if n == '-': gaps += 1
            for n in hsp.sbjct:
                if n == '-': gaps += 1
            for m in hsp.match:
                if m == ' ': mismatches += 1
            outfile.write("\t#hit:" +c[2]["contigName"] +"_"+ str(i) + " strand:" + strand +
                          " contig_pos:" + str(hsp.query_start) +"-"+ str(hsp.query_start + (hsp.align_length * hsp.frame[0]) - 1) +
                          " ref_pos:"+ c[0][i][0] + ":" + str(hsp.sbjct_start) + "-" + str(hsp.sbjct_start + (hsp.align_length * hsp.frame[1]) - 1) +
                          " aligned_length:" + str(hsp.align_length) +
                          " gaps:"+ str(gaps) + " mismatches:" + str(mismatches) + "\n")
            outfile.write("\tquery: " + str(hsp.query) + "\n\tmatch: " + str(hsp.match) + "\n\tsbjct: " + str(hsp.sbjct) + "\n")
            previousHsp = hsp
            previousChr = c[0][i][0]



def main():
    bed = open(args.bed, 'r')
    #header = bed.readline()
    line = bed.readline()
    while len(line) > 3:
        sys.stderr.write("Parsing dipspades\n")
        parseblastout(line, "dipspades")
        sys.stderr.write("Parsing spades\n")
        parseblastout(line, "spades")
        sys.stderr.write("Parsing abyss\n")
        parseblastout(line, "abyss")
        sys.stderr.write("Parsing velvet\n")
        parseblastout(line, "velvet")
        line = bed.readline()

if __name__ == '__main__':
    main()