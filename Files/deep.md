
-) deep resequencing so high depth VCF

SS=/data/data/Software/selscan/bin/linux

$SS/selscan --xpehh --vcf Data/PEL.chr11.vcf --vcf-ref Data/CHB.chr11.vcf --out Results/PEL --map genetic_map_GRCh37_chr11.txt 

$SS/selscan --nsl --vcf Data/PEL.chr11.vcf --out pel --max-extend-nsl 50 --maf 0.04
$SS/norm --ihs --files pel.nsl.out --bins 40

$SS/selscan --xpehh --vcf Data/CLM.chr11.vcf --vcf-ref Data/CHB.chr11.vcf -out Results/CLM.xpehh.txt

echo nsl
$SS/selscan --vcf TLR2.vcf --threads 4 --nsl --out nsl

echo norm
$SS/norm --ihs --files nsl.nsl.out --bins 40 > log

selscan to compute nSL/iHS

-----------------------

plot network with popart (maybe before nSL?)

```
> Results/FADS.fs
Rscript Scripts/vcf2fasta.R Data/PEL.chr11.vcf PEL >> Results/FADS.fa
Rscript Scripts/vcf2fasta.R Data/CLM.chr11.vcf CLM >> Results/FADS.fa
Rscript Scripts/vcf2fasta.R Data/CHB.chr11.vcf CHB >> Results/FADS.fa
```







--------------

OPTIONAL

Indetifying top sites (max PBS)


some final considerations on targeted variants? final aim identifcation of causal variant? what is the functional effect? already associated to some phenotypes? gwas-catalog? ucsc? dgGap? functionality like in roadmap?

------------------------

[HOME](https://github.com/mfumagalli/Weggis)





