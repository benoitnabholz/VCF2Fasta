# These programs convert VCF in Fasta file. 

It has been designed and tested using [Freebayes](https://github.com/freebayes/freebayes)

### VCF2Fasta_no_mono.py : If you don't use the option "--report-monomorphic" of freebayes (*much faster*)

### VCF2FastaLongShot.py : for the [LongShot](https://github.com/pjedge/longshot) software. To genotype long-reads from ONT technology

Both program report all sites *including monomorphic sites*.

They work one individual at the time so extract individual if you have a multi-samples VCF.

For using **VCF2Fasta_no_mono.py** you need :

1) Calculate sequencing depth and identify callable sites:
Calculate the sequencing depth and retain only positions with a coverage within a predefined range considered reliable (here MINCOV and MAXCOV).
Convert the positions from 1-based coordinates (samtools output) to 0-based coordinates to ensure compatibility with BEDtools.
```
samtools depth -a -q 30 -Q 20 -@ 20 bamfile_from_mapping.bam | awk '{for(i=3;i<=NF;i++) t+=$i; print $1, $2, t; t=0}' |  awk -v MiC=$MINCOV -v MaC=$MAXCOV '$3 < MaC && $3 > MiC' |  awk 'OFS="\t"{print $1, $2-1}' >good_sites.txt
```

2) Merge consecutive callable positions into BED intervals.
Convert the list of individual genomic positions into BED intervals by merging consecutive positions belonging to the same chromosome.
```
awk '
BEGIN { OFS="\t" }
NR==1 { chr=$1; start=$2; end=$2; next }
{
    if ($1==chr && $2==end+1) {
        end=$2
    } else {
        print chr, start-1, end
        chr=$1; start=$2; end=$2
    }
}
END { print chr, start-1, end }
' good_sites.txt > good_sites.bed
```

Then **VCF2Fasta_no_mono.py** will 1) create a matrix with the reference sequence, 2) add the SNP with good quality and 
mask the one with low quality and 3) mask site with a "bad" coverage.

** *Warning:* The VCF should not contain haplotype calls or complex alleles. ** Please, use `freebayes ... | vcfallelicprimitives -kg >calls.vcf` or `vt decompose_blocksub` (see https://github.com/freebayes/freebayes#normalizing-variant-representation )

 

*Author:* Benoit Nabholz

--------
### Script to genotype lots of individuals at once using FreeBayes

```
REFGEN=PATH/TO/REF_GENOME
export PATH=$PATH:$HOME/softwares/freebayes/scripts/:$HOME/bin/:$HOME/softwares/freebayes/vcflib/scripts/:$HOME/softwares/freebayes/vcflib/bin/:/media/bigvol/benoit/softwares/speedseq/bin/
```


- compute interval with same coverage
```
rm cov.bash
for chr in $(cat list_chromosomes); do
echo "samtools depth -r $chr -@ 100 -a *ALLBAM*.bam | awk -F'\t'  '{for(i=3;i<=NF;i++) t+=\$i; print \$1\"\\t\"\$2\"\\t\"t; t=0}' | python2.7 coverage_to_regions.py $REFGEN.fai 1000 >regions_cov_$chr.bed" >>cov.bash
done
cat cov.bash | parallel

cat regions_cov_*.bed >regions_cov.bed

rm regions_cov_*.bed

```

- Run Freebayes
```
ulimit -n `ulimit -Hn`
time freebayes-parallel regions_cov.bed 120 -f $REFGEN --use-best-n-alleles 4 -L list_bamAllWGS_cleaned.txt >AllWGS_cleaned.vcf

bgzip AllWGS_cleaned.vcf
```

- correct SNP with multiple bases  (i.e, haplotype calls or complex alleles  = SNP being returned with multiple padding bases : see https://github.com/freebayes/freebayes#normalizing-variant-representation  and https://github.com/freebayes/freebayes/issues/161)
Using VT https://github.com/atks/vt
```
~/bin/vt decompose_blocksub AllWGS_cleaned.vcf.gz >AllWGS_cleaned_vt.vcf
bgzip AllWGS_cleaned_vt.vcf

~/bin/vt uniq AllWGS_cleaned_vt.vcf.gz >tmp && rm AllWGS_cleaned_vt.vcf.gz*
mv tmp AllWGS_cleaned_vt.vcf && bgzip AllWGS_cleaned_vt.vcf
tabix AllWGS_cleaned_vt.vcf.gz
```

- alternatively, you can use vcfallelicprimitives from vcflib (https://github.com/vcflib/vcflib)
```
cat AllWGS_cleaned.vcf | vcfallelicprimitives -kg >calls.vcf
```

- Or, with vcfwaves from vcflib >=1.0.9 : https://github.com/vcflib/vcflib
```
vcfwave AllWGS_cleaned.vcf.gz > AllWGS_cleaned_2.vcf
```


- filter vcf and exclude missing data
```
vcftools --gzvcf AllWGS_cleaned.vcf.gz --out AllWGS_cleaned_goodSNP_noindel_nosexchr --remove-indels --max-missing 1.0 --max-alleles 2 --minQ 200 --minDP 15 --maxDP 150  --recode --recode-INFO-all

python3 ~/bin/FilterVCF.py -R 3 -f 0.2 -s 37 -v AllWGS_cleaned_goodSNP_noindel_nosexchr.recode.vcf >AllWGS_cleaned_goodSNP_noindel_nosexchr.vcf
bgzip AllWGS_cleaned_goodSNP_noindel_nosexchr.vcf
```                   


--------

### usage: VCF2Fasta_no_mono_try.py [-h] [-q QUALITY_THRESHOLD] [-m MIN_COV] [-M MAX_COV] [-R MIN_NUM] [-f MIN_FREQ] [--mask_N] [-v VCF_FILE] [-b BED_FILE] [-r REF_FILE]

  -h, --help            show this help message and exit
  
  -q QUALITY_THRESHOLD, --quality_threshold QUALITY_THRESHOLD
  
  -m MIN_COV, --min_cov MIN_COV
                        Minimum coverage to be genotyped
                        
  -M MAX_COV, --max_cov MAX_COV
                        Maximum coverage to be genotyped
                        
  -R MIN_NUM, --min_num MIN_NUM
                        Minimum number of reads for the alternatice variant allele to be genotyped
                        
  -f MIN_FREQ, --min_freq MIN_FREQ
                        Minimum frequence of minor allele (expected = 0.5 for one diploid individual / default = 0.2)
                        
  --mask_N              Consider N in reference genome as unknown site for all individuals
  
  -v VCF_FILE, --vcf_file VCF_FILE
  
  -b BED_FILE, --bed_good_sites BED_FILE
                       Bed file with good sites
  
  -r REF_FILE, --ref_file REF_FILE
                        Reference genome (fasta)

--------
###     VCF2FastaLongShot.py [-h] [-q QUALITY_THRESHOLD] [-v VCF_FILE] [-c COV_FILE] [-r REF_FILE] [-f MIN_FREQ] [-d DN] [--mask_N]
 
options:
  -h, --help            show this help message and exit
  
  -q QUALITY_THRESHOLD, --quality_threshold QUALITY_THRESHOLD
  
  -v VCF_FILE, --vcf_file VCF_FILE
  
  -c COV_FILE, --cov_file COV_FILE
                        Coverage (Depth) file
                        
  -r REF_FILE, --ref_file REF_FILE
                        Reference genome (fasta)
                        
  -f MIN_FREQ, --min_freq MIN_FREQ
                        Minimum frequence of minor allele (expected = 0.5 for one diploid individual / default = 0.2)
                        
  -d DN, --dn DN        Keep (--dn 0) or exclude (--dn 1, default) the dn SNP
  
  --mask_N              Consider "N" in reference genome as unknown site for all individual

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


