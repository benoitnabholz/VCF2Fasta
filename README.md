# These programs convert VCF in Fasta file. 

It has been designed and tested using [Freebayes](https://github.com/freebayes/freebayes)

### VCF2Fasta.py : works with the option "--report-monomorphic" of freebayes

### VCF2Fasta_no_mono.py : If you don't use the option "--report-monomorphic" of freebayes (*much faster*)

Both program report all sites *including monomorphic sites*.

For using **VCF2Fasta_no_mono.py** you need :

1) to compute the coverage only on good quality bases and mapped reads:
```
samtools depth -q 30 -Q 50 Ind.bam  >Ind.cov.txt
```

2) Extract the postion with a "bad" coverage (either too high or too low, with the threshold >=10x and <170x) :
```
extractBadCovPos.py Ind.cov.txt 9 170 >Ind_bad.cov
```

Alternatively, you can use `awk` straigth out of samtools:
```
min_cov=10
max_cov=170

samtools depth -q 30 -Q 30 Ind.bam | awk -v mincov=${min_cov} -v maxcov=${max_cov} '($3> maxcov || $3<mincov){ print $0}' >Ind_bad.cov
```

3) Then, I convert the VCF in fasta assuming that the positions that are no in the VCF and have a 
"correct" coverage are monophormic. Otherwinse, the site is a "N". 

The program 1) create a matrix with the reference sequence, 2) add the SNP with good quality and 
mask the one with low quality and 3) mask site with a "bad" coverage.

*Author:* Benoit Nabholz

--------

### usage: VCF2Fasta_no_mono_try.py [-h] [-q QUALITY_THRESHOLD] [-m MIN_COV] [-M MAX_COV] [-R MIN_NUM] [-f MIN_FREQ] [--mask_N] [-v VCF_FILE] [-c COV_FILE] [-r REF_FILE]

  -h, --help            show this help message and exit
  
  -q QUALITY_THRESHOLD, --quality_threshold QUALITY_THRESHOLD
  
  -m MIN_COV, --min_cov MIN_COV
                        Minimum coverage to be genotyped
                        
  -M MAX_COV, --max_cov MAX_COV
                        Maximum coverage to be genotyped
                        
  -R MIN_NUM, --min_num MIN_NUM
                        Minimum number of reads for the alternatice variant allele to be genotyped
                        
  -f MIN_FREQ, --min_freq MIN_FREQ
                        Minimum frequence of minor allel (expected = 0.5 for one diploid individual)
                        
  --mask_N              Consider N in reference genome as unknown site for all individuals
  
  -v VCF_FILE, --vcf_file VCF_FILE
  
  -c COV_FILE, --cov_file COV_FILE
                        Coverage (Depth) file
                        
  -r REF_FILE, --ref_file REF_FILE
                        Reference genome (fasta)

--------

### usage: VCF2Fasta.py [-h] [-q QUALITY_THRESHOLD] [-m MIN_COV] [-M MAX_COV] [-R MIN_NUM] [-f VCF_FILE] [--Print_all_positions]


  -q QUALITY_THRESHOLD, --quality_threshold QUALITY_THRESHOLD

  
  -m MIN_COV, --min_cov MIN_COV
                        Minimum coverage to be genotyped


  -M MAX_COV, --max_cov MAX_COV
                        Maximum coverage to be genotyped


  -R MIN_NUM, --min_num MIN_NUM
                        Minimum number of reads for the alternatice variant allele to be genotyped


  -f VCF_FILE, --vcf_file VCF_FILE

  -c COV_FILE, --cov_file COV_FILE
                        Coverage (Depth) file

  -r REF_FILE, --ref_file REF_FILE
                        Reference genome (fasta)
-----

## Example file
A tool dataset is present in the directory data. This dataset have been created using the data of [Delmore et al. 2020](https://elifesciences.org/articles/54462).


To use the program:
``` 
python3 ./VCF2Fasta.py -q 10 -m 10 -M 100 -R 2 -f data/try.vcf
```

The output should be
```
CM020465.1.fst
```
, containing the 30 sequences (two per individuals).


