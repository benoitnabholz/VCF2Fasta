#! /usr/bin/env python
from __future__ import print_function
import sys
import argparse
import textwrap
import numpy as np
import random

def read_fasta(fp):
	name, seq = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith(">"):
			if name: yield (name, ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))

# Function to find information in position in FORMAT field 
def findPosInFORMATStr(infoStr, str2ID):
	RETURNPOS = 0
	findIt = "no"
	info = infoStr.split(":")
	for i in info:
		if i == str2ID:
			findIt ="yes"
			break
		RETURNPOS = RETURNPOS+1
	if findIt == "yes":
		return(RETURNPOS)
	else:
		return(-999)
		
def minAllFreqCompute(ad):
	REF_cov=float(ad[0])
	ALT_cov=float(ad[1])
	min_all_freq = ALT_cov/(REF_cov+ALT_cov)
	if min_all_freq > 0.5:
		min_all_freq = 1-min_all_freq
	return(min_all_freq)
	
	
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
description=textwrap.dedent('''\
The program convert VCF in Fasta file. 
It has been designed and tested using Freebayes https://github.com/freebayes/freebayes

If you don't use the option --report-monomorphic then freebayes only output putative SNP or
structural variant. 

In this case, I use a the following strategy:

1) compute the coverage :
sambamba depth base -L $contig -t 2 Ind.bam  >Ind.cov.txt

2) Extract the postion with a wrong coverage either too high or too low :
~/bin/extractBadCovPos.py Ind.cov.txt 9 170 >Ind_bad.cov

3) convert the VCF in fasta file assuming that the positions that are not in the VCF and have a 
"correct" coverage are monophormic. Otherwise, the site is reported as unknown (i.e, with the 
"N" genotype).

The program 1) create a sequence with the reference sequence, 2) add the SNP with good quality and 
mask the one with low quality and 3) mask site with a "bad" coverage.


Author: Benoit Nabholz
Github: https://github.com/benoitnabholz/VCF2Fasta
''')
)



parser.add_argument('-q', "--quality_threshold",  default=10.0, type=float)
parser.add_argument('-m', '--min_cov', type=int, help="Minimum coverage to be genotyped")
parser.add_argument('-M', '--max_cov', type=int, help="Maximum coverage to be genotyped")
parser.add_argument('-R', '--min_num', type=int, help="Minimum number of reads for the alternatice variant allele to be genotyped") # Variant : min. number of reads to be retained 
parser.add_argument('-f', '--min_freq', type=float, default=0.2, help="Minimum frequence of minor allele (expected = 0.5 for one diploid individual)") 

parser.add_argument('--mask_N', dest='mask_N',  default=False, action="store_true", help="Consider N in reference genome as unknown site for all individuals")

parser.add_argument('-v', '--vcf_file')
parser.add_argument('-c', '--cov_file', help="Coverage (Depth) file")
parser.add_argument('-r', '--ref_file', help="Reference genome (fasta)")



args = parser.parse_args()

fp = open(str(args.ref_file))
vcf = open(str(args.vcf_file))
cov_f = open(str(args.cov_file))

minCoverageThreshold = int(args.min_cov)
maxCoverageThreshold = int(args.max_cov)
minNumberReadsThreshold = int(args.min_num)
minAllFreqThreshold = float(args.min_freq)
QUALThres = float(args.quality_threshold)
print("Quality threshold is set to "+str(QUALThres))
mask_N=args.mask_N

# Store fasta file in dictionary
RefGen = {}
for name, seq in read_fasta(fp):
	name = name.lstrip(">")
	RefGen[name] = seq
fp.close()
print("Ref genome loaded")

# Store cov in dict
Cov = {}
pos = 0
seqName = ""
for line in cov_f:
	line = line.rstrip()
	arline = line.split("\t")
	pos=int(arline[1])
	
	if arline[0] not in Cov:
		p = [pos]
		Cov[arline[0]] = p
	else:
		Cov[arline[0]].append(pos)
cov_f.close()
print("Cov file loaded")

# Read VCF
nameSeq = []
nameInd = []
sequenceDNA = ""
countOutFasta = 0
referenceSeq = ""
genotypedSeq = []
seqName = ""


countSNP = 0
countSite = 0
NumberOfSeq = ""
sites = []

TotCov = 0

pos = -1
oldpos = 0
logfile = open("logfile_VCF2Fasta.txt", "w")

count = 0

for line in vcf:

	# count number of individuals and store names in array
	if line[0:6] == "#CHROM" and len(nameInd) == 0:
		arrline = line.split()
		
		# individuals are named using the popphyl format 
		for i in arrline[9:]:
			indseq1 = i+"|Allele1"
			indseq2 = i+"|Allele2"
			nameSeq.append(indseq1)
			nameSeq.append(indseq2)
			nameInd.append(i)
		NumberOfSeq = len(nameSeq)
		print("There is "+str(NumberOfSeq/2)+" individuals")


	if line[0] == "#":
		continue
	
	arrline = line.split()
	chromosome = arrline[0]
	pos = int(arrline[1])
	REF = str(arrline[3])
	ALT = str(arrline[4])
	site = []
	
	## Create the matrix with reference sequence
	if seqName == "":
		seqName = chromosome
		referenceSeq = list(RefGen[seqName])
		genotypedSeq = np.empty((0,len(referenceSeq)))

		for i in range(0, len(nameSeq)):
			seq = np.array(referenceSeq, dtype=object)
			genotypedSeq = np.append(genotypedSeq, [seq], axis=0)


	if chromosome != seqName and seqName != "":

		# put 'N' in low or high cov sites
		if seqName in Cov:
			for bad_pos in Cov[seqName]:
				for i in range(0, len(nameSeq)):
					seqNumber = i*2 # convert individual number in seq number
					genotypedSeq[i][bad_pos-1]="N"

		# print sequences
		fasta = open(seqName+".fst", "w")
		for i in range(0, len(nameSeq)):
			print(">"+nameSeq[i], file=fasta)
			print("".join(genotypedSeq[i]), file=fasta)

		## set chromosome name
		seqName = chromosome

		## Create the matrix with reference sequence
		referenceSeq = list(RefGen[seqName])
		genotypedSeq = np.empty((0,len(referenceSeq)))
		for i in range(0, len(nameSeq)):
			seq = np.array(referenceSeq, dtype=object)
			genotypedSeq = np.append(genotypedSeq, [seq], axis=0)


	if len(REF) > 1 or len(ALT) > 3 or ALT == "*" or (len(ALT) == 3 and "," not in ALT) :  # exlcude indels
		
		for i in range(0, len(nameSeq)):
			genotypedSeq[i][pos-1]="N"

		continue

	if referenceSeq[pos-1] == "N" and mask_N :  # exclude position with "N" in the reference
		for i in range(0, len(nameSeq)):
			genotypedSeq[i][pos-1]="N"
		continue

	covPos = findPosInFORMATStr(arrline[8], "DP") # find coverage position
	gtPos  = findPosInFORMATStr(arrline[8], "GT")

	if  QUALThres > float(arrline[5]): # exclude lowQ SNP 

		if len(ALT) == 3:  # exlcude three alleles cases
			for i in range(0, len(nameSeq)):
				genotypedSeq[i][pos-1]="N"
		else:
			for i in range(0, len(nameInd)):
				seqNumber = i*2 # convert individual number in seq number

				ind = arrline[i+9]
				if ind == "." : # no reads for this individual
					genotypedSeq[seqNumber][pos-1]="N"
					genotypedSeq[seqNumber+1][pos-1]="N"
					continue
						
				arrInd = ind.split(":")
				if arrInd[covPos] == ".":
					genotypedSeq[seqNumber][pos-1]="N"
					genotypedSeq[seqNumber+1][pos-1]="N"
					continue

				cov = int(arrInd[covPos])
				GT = str(arrInd[gtPos]) # genotype information
					
				if cov < minCoverageThreshold or cov > maxCoverageThreshold:
					genotypedSeq[seqNumber][pos-1]="N"
					genotypedSeq[seqNumber+1][pos-1]="N"
					continue
							
				if GT == "0/1":
					genotypedSeq[seqNumber][pos-1]="N"
					genotypedSeq[seqNumber+1][pos-1]="N"
				if GT == "0/0":
					genotypedSeq[seqNumber][pos-1]=REF
					genotypedSeq[seqNumber+1][pos-1]=REF
				if GT == "1/1":
					genotypedSeq[seqNumber][pos-1]="N"
					genotypedSeq[seqNumber+1][pos-1]="N"
	else:

		covPos = findPosInFORMATStr(arrline[8], "DP") # find coverage position
		gtPos  = findPosInFORMATStr(arrline[8], "GT")
					
		ALT1 = "N"
		ALT2 = "N"
		if len(ALT) == 3:
			ALT1 = ALT.split(",")[0]
			ALT2 = ALT.split(",")[1]
			ALT = ALT.split(",")[0] # Security if I forgot something below

					
		## Option possible??
		##INFO=<ID=AO,Number=A,Type=Integer,Description="Count of full observations of this alternate haplotype.">
		##INFO=<ID=SAF,Number=A,Type=Integer,Description="Number of alternate observations on the forward strand">
		##INFO=<ID=SAR,Number=A,Type=Integer,Description="Number of alternate observations on the reverse strand">
		#print(line)
					
		for i in range(0, len(nameInd)):

			seqNumber = i*2 # convert individual number in seq number

			ind = arrline[i+9]
			if ind == ".": # no reads for this individual
				genotypedSeq[seqNumber][pos-1]="N"
				genotypedSeq[seqNumber+1][pos-1]="N"
				continue

			arrInd = ind.split(":")
			if arrInd[covPos] == ".":
				genotypedSeq[seqNumber][pos-1]="N"
				genotypedSeq[seqNumber+1][pos-1]="N"
				continue

			cov = int(arrInd[covPos])
			GT = str(arrInd[gtPos]) # genotype information

			if cov < minCoverageThreshold or cov > maxCoverageThreshold:
				genotypedSeq[seqNumber][pos-1]="N"
				genotypedSeq[seqNumber+1][pos-1]="N"
				continue
	
			if GT == "0/0":
				genotypedSeq[seqNumber][pos-1]=REF
				genotypedSeq[seqNumber+1][pos-1]=REF
							
			if GT == "0/1":
				adPos = findPosInFORMATStr(arrline[8], "AD") #D=AD,Number=R,Type=Integer,Description="Number of observation for each allele 
				ad = arrInd[adPos].split(",")
				
				min_all_freq = minAllFreqCompute(ad)
					
				if int(ad[0]) > minNumberReadsThreshold and int(ad[1]) > minNumberReadsThreshold and min_all_freq >= minAllFreqThreshold: # must be supported by at least minNumberReadsThreshold reads and coverage higher than min allele freq
					genotypedSeq[seqNumber][pos-1]=REF
					genotypedSeq[seqNumber+1][pos-1]=ALT
				else:
					genotypedSeq[seqNumber][pos-1]="N"
					genotypedSeq[seqNumber+1][pos-1]="N"
				
			if GT == "1/1":
				if len(ALT) == 3:
					genotypedSeq[seqNumber][pos-1]=ALT1
					genotypedSeq[seqNumber+1][pos-1]=ALT1
				else:
					genotypedSeq[seqNumber][pos-1]=ALT
					genotypedSeq[seqNumber+1][pos-1]=ALT
							
			if GT == "0/2":

				adPos = findPosInFORMATStr(arrline[8], "AD") #D=AD,Number=R,Type=Integer,Description="Number of observation for each allele 
				ad = arrInd[adPos].split(",")
				min_all_freq = minAllFreqCompute(ad)
			
				if int(ad[0]) > minNumberReadsThreshold and int(ad[2]) > minNumberReadsThreshold and min_all_freq >= minAllFreqThreshold: # must be supported by at least XX reads
					genotypedSeq[seqNumber][pos-1]=REF
					genotypedSeq[seqNumber+1][pos-1]=ALT2
				else:
					genotypedSeq[seqNumber][pos-1]="N"
					genotypedSeq[seqNumber+1][pos-1]="N"

			if GT == "2/2":
				genotypedSeq[seqNumber][pos-1]=ALT2
				genotypedSeq[seqNumber+1][pos-1]=ALT2


			if GT == "1/2":
				adPos = findPosInFORMATStr(arrline[8], "AD") #D=AD,Number=R,Type=Integer,Description="Number of observation for each allele 
				ad = arrInd[adPos].split(",")
				min_all_freq = minAllFreqCompute(ad)
				
				if int(ad[1]) > minNumberReadsThreshold and int(ad[2]) > minNumberReadsThreshold and min_all_freq >= minAllFreqThreshold: # must be supported by at least XX reads
					genotypedSeq[seqNumber][pos-1]=ALT1
					genotypedSeq[seqNumber+1][pos-1]=ALT2
				else:
					genotypedSeq[seqNumber][pos-1]="N"
					genotypedSeq[seqNumber+1][pos-1]="N"

# put 'N' in low or high cov sites
if seqName in Cov:
	for bad_pos in Cov[seqName]:
		for i in range(0, len(nameSeq)):
			seqNumber = i*2 # convert individual number in seq number
			genotypedSeq[i][bad_pos-1]="N"


fasta = open(seqName+".fst", "w")
for i in range(0, len(nameSeq)):
	print(">"+nameSeq[i], file=fasta)
	print("".join(genotypedSeq[i]), file=fasta)			
