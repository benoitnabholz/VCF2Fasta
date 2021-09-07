#! /usr/bin/env python
### Version du Janvier 2021
### Benoit Nabholz benoit.nabholz@umontpellier.fr

from __future__ import print_function
import sys
import numpy as np
import optparse
import argparse
import glob
import textwrap

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
	


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
description=textwrap.dedent('''\
The program convert VCF in Fasta file. 
It has been designed and tested using Freebayes https://github.com/freebayes/freebayes
Author: Benoit Nabholz
Github: 
''')
)

parser.add_argument('-q', "--quality_threshold",  default=10.0, type=float)
parser.add_argument('-m', '--min_cov', type=int, help="Minimum coverage to be genotyped")
parser.add_argument('-M', '--max_cov', type=int, help="Maximum coverage to be genotyped")
parser.add_argument('-R', '--min_num', type=int, help="Minimum number of reads for the alternatice variant allele to be genotyped") # Variant : min. number of reads to be retained 
parser.add_argument('-f', '--vcf_file')
parser.add_argument('--Print_all_positions', default=False, action='store_true', help="Print all positions corresponding to the reference scaffolds (i.e., from start to end)")

args = parser.parse_args()

minCoverageThreshold = int(args.min_cov)
maxCoverageThreshold = int(args.max_cov)
minNumberReadsThreshold = int(args.min_num)

# Quality threshold QUAL field
QUALThres = float(args.quality_threshold)
print("Quality threshold is set to "+str(QUALThres))

namefile = str(args.vcf_file)
fp = open(namefile)

############## MAIN ############
nameSeq = []
nameInd = []
sequenceDNA = ""
countOutFasta = 0

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

## Dict for Scaffolds size
Scaffold_size = {}

for line in fp:
	line = line.rstrip()
	
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

	# Store all scaffolds size in dictionary
	if line[0:8] == "##contig":
		arline=line.split(",")
		Scaf_ID = arline[0]
		Scaf_ID = Scaf_ID.replace("##contig=<ID=", "")
		Scaf_size = arline[1].rstrip(">")
		Scaf_size = Scaf_size.replace("length=","")
		Scaffold_size[Scaf_ID] = int(Scaf_size)
		continue

	if line[0] != "#":
		
		arrline = line.split()
		chromosome = arrline[0]
		pos = int(arrline[1])
		REF = str(arrline[3])
		ALT = str(arrline[4])
		site = []
		#oldsite = []
		

		if (pos % 100000) == 0:
			print("Positions "+str(pos)+" of "+chromosome)

		
		if chromosome != seqName and seqName != "": # print seq when new chromosome found

			#countOutFasta=glob.glob(seqName+'*.fst')
			#countOutFasta=len(countOutFasta)+1
			Scaf_size = Scaffold_size[seqName]
			if len(sites) < Scaf_size:
				numberSite = (int(Scaf_size)+1) - len(sites)
				for j in range(1,numberSite):
					site = []
					for i in range(0, len(nameSeq)):
						site.append("N")
					sites.append(site)
					site = []
			fasta = open(seqName+".fst", "w")
			sites = np.ravel(sites)
			sequenceDNA = np.reshape(sites, newshape=(int(len(sites)/NumberOfSeq), NumberOfSeq))
			for i, seqName in enumerate(nameSeq):
					print(">"+seqName, file=fasta)
					print(''.join(map(str, sequenceDNA[0:,i])), file=fasta)
					
			seqName = chromosome
			sites = []
			oldpos=0
			
			
		if seqName == "":
			seqName = chromosome
		
		######## Loop if current site is not exactly one position after oldpos site put "N" ########################## 
		## Usefull for dealing with indels
		if oldpos+1 < pos and pos != -1 :
			if args.Print_all_positions and oldpos == 0: # if first position of VCF not 1.
				numberSite = int(pos-oldpos)
				logfile.write("Adding position at: "+str(pos)+"\n"+"Number of position added: "+str(numberSite))
				for j in range(1,numberSite):
					site = []
					for i in range(0, len(nameSeq)):
						site.append("N")
					sites.append(site)
					site = []
			elif oldpos != 0:
				numberSite = int(pos-oldpos)
				logfile.write("Adding position at: "+str(pos)+"\n"+"Number of position added: "+str(numberSite))
				for j in range(1,numberSite):
					site = []
					for i in range(0, len(nameSeq)):
						site.append("N")
					sites.append(site)
					site = []
		################################################################################################################


		if len(REF) > 1 or len(ALT) > 3 or ALT == "*" or (len(ALT) == 3 and "," not in ALT) :  # exlcude indels
			for i in range(0, len(nameSeq)):
				site.append("N")
			
		else:
			covPos = findPosInFORMATStr(arrline[8], "DP") # find coverage position
			gtPos  = findPosInFORMATStr(arrline[8], "GT")
			
			if ALT == ".": # if no ALT allele ( = monorphic site)
				for i in range(0, len(nameInd)):
					ind = arrline[i+9]
					if ind == "." : # no reads for this individual
						site.append("N")
						site.append("N")
						continue

					arrInd = ind.split(":")
					if arrInd[covPos] == ".":
						site.append("N")
						site.append("N")
						continue

					cov = int(arrInd[covPos])
									
					if cov > minCoverageThreshold and cov < maxCoverageThreshold:
					#	print("OK" + "\t" +str(ind)+ "\t" +str(cov)+ "\t" +line)
						site.append(REF)
						site.append(REF)
					else:
					#	print("NO" + "\t" +str(ind)+ "\t" +str(cov)+ "\t" +line)
						site.append("N")
						site.append("N")
			else:
				if  QUALThres >  float(arrline[5]): # exclude lowQ SNP / only look for homozygote

					if len(ALT) == 3:  # exlcude three alleles cases
						for i in range(0, len(nameSeq)):
							site.append("N")
					else:
						for i in range(0, len(nameInd)):
							ind = arrline[i+9]
							if ind == "." : # no reads for this individual
								site.append("N")
								site.append("N")
								continue
						
							arrInd = ind.split(":")
							if arrInd[covPos] == ".":
						        site.append("N")
						        site.append("N")
						        continue
							cov = int(arrInd[covPos])
							GT = str(arrInd[gtPos]) # genotype information
					
							if cov < minCoverageThreshold or cov > maxCoverageThreshold:
								#print("POLY NO" + "\t" +str(ind)+ "\t" +str(cov)+ "\t" +line)
								site.append("N")
								site.append("N")
								continue
							
							if GT == "0/1":
								site.append("N")
								site.append("N")
							if GT == "0/0":
								site.append(REF)
								site.append(REF)
							if GT == "1/1":
								site.append("N")
								site.append("N")

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
						ind = arrline[i+9]
						if ind == ".": # no reads for this individual
							site.append("N")
							site.append("N")
							continue

						arrInd = ind.split(":")
						if arrInd[covPos] == ".":
						        site.append("N")
						        site.append("N")
						        continue

						cov = int(arrInd[covPos])
						GT = str(arrInd[gtPos]) # genotype information
						
						if cov < minCoverageThreshold or cov > maxCoverageThreshold:
							site.append("N")
							site.append("N")
							continue
						
						if GT == "0/0":
							site.append(REF)
							site.append(REF)
							
						if GT == "0/1":
							adPos = findPosInFORMATStr(arrline[8], "AD") #D=AD,Number=R,Type=Integer,Description="Number of observation for each allele 
							ad = arrInd[adPos].split(",")
							if int(ad[0]) > minNumberReadsThreshold and int(ad[1]) > minNumberReadsThreshold: # must be supported by at least XX reads
								site.append(REF)
								site.append(ALT)
							else:
								site.append("N")
								site.append("N")
					
						if GT == "1/1":
							if len(ALT) == 3:
								site.append(ALT1)
								site.append(ALT1)
							else:
								site.append(ALT)
								site.append(ALT)
							
						if GT == "0/2":
							adPos = findPosInFORMATStr(arrline[8], "AD") #D=AD,Number=R,Type=Integer,Description="Number of observation for each allele 
							ad = arrInd[adPos].split(",")
							if int(ad[0]) > minNumberReadsThreshold and int(ad[2]) > minNumberReadsThreshold: # must be supported by at least XX reads
								site.append(REF)
								site.append(ALT2)
							else:
								site.append("N")
								site.append("N")
								
						if GT == "2/2":
							site.append(ALT2)
							site.append(ALT2)
							
						if GT == "1/2":
							adPos = findPosInFORMATStr(arrline[8], "AD") #D=AD,Number=R,Type=Integer,Description="Number of observation for each allele 
							ad = arrInd[adPos].split(",")
							if int(ad[1]) > minNumberReadsThreshold and int(ad[2]) > minNumberReadsThreshold: # must be supported by at least XX reads
								site.append(ALT1)
								site.append(ALT2)
							else:
								site.append("N")
								site.append("N")
							
						
		sites.append(site)
		## Test if the new site is of correct length [ WARNING : This test take a huge amount of time ]
		#print(line)
		#print(str(pos))


		#if not (float(len(site))/float(NumberOfSeq)).is_integer():
		#	print(line)
		#	print(site)
		#	print(chromosome)
		#	print(str(pos))
		#	sys.exit("NOOOO :-)")
		
		oldpos = pos
		
		
if args.Print_all_positions:
	Scaf_size = Scaffold_size[seqName]
	if len(sites) < Scaf_size:
		print("Add site at the end.")
		numberSite = (int(Scaf_size)+1) - len(sites)
		for j in range(1,numberSite):
			site = []
			for i in range(0, len(nameSeq)):
				site.append("N")
			sites.append(site)
			site = []

fasta = open(seqName+".fst", "w")
sites = np.ravel(sites)
sequenceDNA = np.reshape(sites, newshape=(int(len(sites)/NumberOfSeq), NumberOfSeq))
for i, seqName in enumerate(nameSeq):
		print(">"+seqName, file=fasta)
		print(''.join(map(str, sequenceDNA[0:,i])), file=fasta)


#print "Number of SNP "+str(count)

fp.close()

