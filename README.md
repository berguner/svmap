# SVMap
SVMap is a tool for fine mapping the genomic structural variations of interest.

## Python requirements
- python2.7
- matplotlib
- biopython
## External dependencies
- BWA aligner (installed or in the $PATH)
- samtools (installed or in the $PATH)
- NCBI Blast (installed or in the $PATH)
- SPAdes assembler
- velvet assembler
- ABySS assembler
- Reference genome fasta file

## Quick Start
1. Make sure indexes of the genome for BWA and Blastn are in the same directory with the reference fasta
2. Prepare a tab delimited .bed file listing interested regions
    "CHR START	END	Info1	Info2	Info3"
3. Extract the reads by running "extractreads_fq.py"
4. If you want faster results use "extractreads_bam.py". Note that this will not extract unmapped reads.
5. Assemble the extracted reads with "assemble.py"
6. Blast the assembled contigs with "blast.py"
7. Parse and plot the results by "svmap.py"
