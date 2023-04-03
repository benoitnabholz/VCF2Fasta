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
	
def minAllFreqCompute(InfoField):
	info = InfoField.split(";")
	for i in info:
		if "AC=" in i:
			COV=i.split("=")[1].split(",")
			All1=float(COV[0])
			All2=float(COV[1])
			break
	min_all_freq = All1/(All1+All2)
	if min_all_freq > 0.5:
		min_all_freq = 1-min_all_freq
	return(min_all_freq)

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
description=textwrap.dedent('''\
The program convert VCF in Fasta file. 
VCF obtained using longshot https://github.com/pjedge/longshot

Author: Benoit Nabholz
Github: https://github.com/benoitnabholz/VCF2Fasta
''')
)



parser.add_argument('-q', "--quality_threshold",  default=10.0, type=float)
parser.add_argument('-v', '--vcf_file')
parser.add_argument('-c', '--cov_file', help="Coverage (Depth) file")
parser.add_argument('-r', '--ref_file', help="Reference genome (fasta)")
parser.add_argument('-f', '--min_freq', type=float, default=0.2, help="Minimum frequence of minor allele (expected = 0.5 for one diploid individual)") 
parser.add_argument('-d', '--dn', default=1 ,type=int, help="Keep (--dn 0) or exclude (--dn 1, default) the dn SNP")
parser.add_argument('--mask_N', dest='mask_N',  default=False, action="store_true", help="Consider \"N\" in reference genome as unknown site for all individuals")


args = parser.parse_args()
# Quality threshold QUAL field
QUALThres = float(args.quality_threshold)
fp = open(str(args.ref_file))
vcf = open(str(args.vcf_file))
cov_f = open(str(args.cov_file))
minAllFreqThreshold = float(args.min_freq)
dn = int(args.dn)
mask_N=args.mask_N


print("Quality threshold is set to "+str(QUALThres))

# Store fasta file in dictionary
RefGen = {}
RefGenPrinted = {}
for name, seq in read_fasta(fp):
	name = name.lstrip(">")
	RefGen[name] = seq
	RefGenPrinted[name] = 0
fp.close()
print("Ref genome loaded")

# Store cov in dict
Cov = {}
pos = 0
scaf = ""
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
pos=0
scaf=""
All1=""
All2=""
for line in vcf:
	if line[0] == "#":
		continue

	arline = line.split("\t")
	if scaf == "":
		scaf = arline[0]
		All1 = list(RefGen[scaf])
		All2 = list(RefGen[scaf])
		refgen = list(RefGen[scaf])

	if scaf != arline[0]:
		# put 'N' in low or high cov sites
		if scaf in Cov:
			for bad_pos in Cov[scaf]:
				All1[bad_pos-1] = "N"
				All2[bad_pos-1] = "N"
			
		# Print seq
		fasta = open(scaf+".fst", "w")
		print(">"+scaf+"|All1"+"\n"+"".join(All1), file=fasta)
		print(">"+scaf+"|All2"+"\n"+"".join(All2), file=fasta)

		# new seq
		scaf = arline[0]
		All1 = list(RefGen[scaf])
		All2 = list(RefGen[scaf])
		refgen = list(RefGen[scaf])
		RefGenPrinted[scaf] = 1

	Ref = arline[3]
	Alt = arline[4]
	pos = int(arline[1])
	InfoField = arline[7]
	
	if len(Alt) > 1:
		print(line)

	# if "N" in ref genome (e.g., masked for repeat)
	if refgen[pos-1] == "N" and mask_N :
		All1[pos-1] = "N"
		All2[pos-1] = "N"
		continue
		
	if All1[pos-1] != Ref:
		print("POBLEM\nRef allele does not correspond to "+name+" sequence\t"+line)
		All1[pos-1] = "N"
		All2[pos-1] = "N"
		sys.exit()


	# Quality
	if float(arline[5]) < QUALThres:
		#print(line)
		All1[pos-1] = "N"
		All2[pos-1] = "N"
		continue

	# Exclude dn position
	if arline[6] == "dn" and dn == 1:
		All1[pos-1] = "N"
		All2[pos-1] = "N"
		continue

	if arline[6] == "dn;dp":
		All1[pos-1] = "N"
		All2[pos-1] = "N"
		continue

	if arline[6] == "dp":
		All1[pos-1] = "N"
		All2[pos-1] = "N"
		continue
		
	if  minAllFreqThreshold > minAllFreqCompute(InfoField):
		All1[pos-1] = "N"
		All2[pos-1] = "N"
		continue
		
	GT = arline[9].split(":")[0]
	if GT == "1|0":
		All1[pos-1] = Alt
	if GT == "0|1":
		All2[pos-1] = Alt
	if GT == "1|1":
		All1[pos-1] = Alt
		All2[pos-2] = Alt

	if GT == "0/1":
		if random.uniform(0, 1) >= 0.5:
			All1[pos-1] = Alt
		else:
			All2[pos-1] = Alt

	if GT == "1/1":
		All1[pos-1] = Alt
		All2[pos-1] = Alt

RefGenPrinted[scaf] = 1
# put 'N' in low or high cov sites
if scaf in Cov:
	for bad_pos in Cov[scaf]:
		All1[bad_pos-1] = "N"
		All2[bad_pos-1] = "N"
	 
# Print seq
fasta = open(scaf+".fst", "w")
print(">"+scaf+"|All1"+"\n"+"".join(All1), file=fasta)
print(">"+scaf+"|All2"+"\n"+"".join(All2), file=fasta)


# take care of the contig that were not present in the VCF (e.g., contig without SNP)
for scaf in RefGenPrinted:
	if RefGenPrinted[scaf] == 0:
		# new seq
		scaf = arline[0]
		All1 = list(RefGen[scaf])
		All2 = list(RefGen[scaf])
		refgen = list(RefGen[scaf])
		
		# put 'N' in low or high cov sites
		if scaf in Cov:
			for bad_pos in Cov[scaf]:
				All1[bad_pos-1] = "N"
				All2[bad_pos-1] = "N"
				
	 	# Print seq
		fasta = open(scaf+".fst", "w")
		print(">"+scaf+"|All1"+"\n"+"".join(All1), file=fasta)
		print(">"+scaf+"|All2"+"\n"+"".join(All2), file=fasta)

	


