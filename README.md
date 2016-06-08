# BAM2GENES

This code is developed for identification of genes from NGS data of pooled PCR amplicons.

1. map the fq file to ref genome (bam file)
2. sort and index the bam file

**Command**
```sh
python 1.continuity_top_three.py [bamfile] [ref.fa]
```
**Result**
```
$ python 1.continuity_top_three.py T2-1_C09_S33.sam.sorted.bam /ref/analysis/References/Creinhardtii/Creinhardtii_281_v5.0.fa
Filename First_Genename Cov1x Cov10x  Cov30x  Depth  Total.depth  Second_Genename Cov1x Cov10x  Cov30x  Depth  Total.depth Third_Genename Cov1x Cov10x  Cov30x  Depth  Total.depth
T2-1_C09_S33.sam.sorted.bam     Cre17.g704850.v5.5      0.994505494505  0.644688644689  0.0     12.4175824176   6780    Cre14.g617950.v5.5      0.793304221252  0.647743813683  0.0     8.91848617176   6127    Cre06.g280900.v5.5 0.996138996139  0.0978120978121 0.0     6.85070785071   5323
```
- Total.depth : sum(depth per base on gene region)
- Depth : Total.depth/genelength


