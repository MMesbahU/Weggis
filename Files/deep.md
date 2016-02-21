
-) deep resequencing so high depth VCF

SS=/data/data/Software/selscan/bin/linux

#$SS/selscan --xpehh --vcf Data/PEL.chr11.vcf --vcf-ref Data/CHB.chr11.vcf --out Results/PEL --map genetic_map_GRCh37_chr11.txt 

explain options in nsl

```
$SS/selscan --nsl --vcf Data/PEL.chr11.vcf --out Results/PEL --max-extend-nsl 200 --maf 0.04
$SS/norm --ihs --files Results/PEL.nsl.out --bins 20
```
The output file is called `Results/PEL.nsl.out.20bins.norm`.

Plot the results:
```
Rscript Scripts/plotnSL.R Results/PEL.nsl.out.20bins.norm Results/PEL.nsl.pdf
```

Why not high?
Investigate haplotypes.

-----------------------


```
> Results/FADS.fa
Rscript Scripts/vcf2fasta.R Data/PEL.fads.vcf PEL Results/PEL.fads.snp >> Results/FADS.fa
Rscript Scripts/vcf2fasta.R Data/CLM.fads.vcf CLM NULL >> Results/FADS.fa
Rscript Scripts/vcf2fasta.R Data/CHB.fads.vcf CHB NULL >> Results/FADS.fa
```
Have a look at the resulting file:
```
less -S Results/FADS.fa
```

Plot the haplotype network.
```
Rscript Scripts/plotNet.R Results/FADS.fa Results/PEL.fads.snp Results/FADS.pdf 2> /dev/null > Results/FADS.diff
```
Open the plot:
```
evince Results/FADS.pdf
```

Look at candidate causal variants:
```
less -S Results/FADS.diff
grep -P "V \t XXXVI" Results/FADS.diff | head -n 2 | tail -n 1 > Results/FADS.cause.diff
```

Another tool for visualising haplotypes is popart (give link).

--------------

OPTIONAL

Calculate allele frequencies in CHB, CLM and PEL for the previously identified causal variants.
Calculate PBS.

some final considerations on targeted variants? final aim identifcation of causal variant? what is the functional effect? already associated to some phenotypes? gwas-catalog? ucsc? dgGap? functionality like in roadmap?

------------------------

[HOME](https://github.com/mfumagalli/Weggis)





