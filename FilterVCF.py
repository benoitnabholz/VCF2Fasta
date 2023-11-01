#! /usr/bin/env python3
from __future__ import print_function
import sys
import argparse
import textwrap
import numpy as np
import random

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
Clean VCF based on criteria not available in BCFTOOLS.
Author: Benoit Nabholz
''')
)


parser.add_argument('-R', '--min_num', type=int, help="Minimum number of reads for the alternatice variant allele to be genotyped") # Variant : min. number of reads to be retained 
parser.add_argument('-f', '--min_freq', type=float, default=0.2, help="Minimum frequence of minor allele (expected = 0.5 for one diploid individual)") 
parser.add_argument('-s', '--min_samples', type=int, help="Number of samples in the VCF") 
parser.add_argument('-v', '--vcf_file')



args = parser.parse_args()

vcf = open(str(args.vcf_file))
minNumberReadsThreshold = int(args.min_num)
minAllFreqThreshold = float(args.min_freq)
num_sample = int(args.min_samples)

# Read VCF
nameSeq = []
nameInd = []

countSNP = 0
countSite = 0
NumberOfSeq = ""
sites = []


pos = -1
oldpos = 0
logfile = open("logfile_VCF2Fasta.txt", "w")

count = 0

for line in vcf:

	line = line.rstrip()
	
	if line[0] == "#":
		print(line)
		continue
	
	arrline = line.split()
	chromosome = arrline[0]
	
	
	covPos = findPosInFORMATStr(arrline[8], "DP") # find coverage position
	gtPos  = findPosInFORMATStr(arrline[8], "GT")

	printSite = "YES"
	
	for i in range(0, num_sample):

		ind = arrline[i+9]
		if ind == "." : # no reads for this individual
			continue
						
		arrInd = ind.split(":")
		if arrInd[covPos] == ".":
			continue

		cov = int(arrInd[covPos])
		GT = str(arrInd[gtPos]) # genotype information
					
		if GT == "0/1":
			adPos = findPosInFORMATStr(arrline[8], "AD") #D=AD,Number=R,Type=Integer,Description="Number of observation for each allele 
			ad = arrInd[adPos].split(",")
				
			min_all_freq = minAllFreqCompute(ad)
					
			if int(ad[0]) > minNumberReadsThreshold and int(ad[1]) > minNumberReadsThreshold and min_all_freq >= minAllFreqThreshold: # must be supported by at least minNumberReadsThreshold reads and coverage higher than min allele freq
				#printSite = "YES"
				continue
			else:
				#print(line+" "+str(min_all_freq))
				printSite = "NO"
	
	
	if printSite == "YES":
		print(line)
		
		
