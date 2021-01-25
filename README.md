# The program convert VCF in Fasta file. 

It has been designed and tested using [Freebayes](https://github.com/freebayes/freebayes)

*Author:* Benoit Nabholz

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


  --Print_all_positions
                        Print all positions corresponding to the reference scaffolds (i.e., from start to end)

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


