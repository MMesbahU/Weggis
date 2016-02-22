
We have some suggestive results pointing towards selection acting in the FADS region based on patterns of genetic differentiation against East Asians and Europeans.
However, we employed low-depth data and our analyses were not based on called genotypes.
Therefore, we did not have the power, for instance, to identify the causal variant (if any) and patterns of haplotype distribution.

It is usual to then perform a targeted deep resequencing of our region of interest, to further refione our selection analyses and highlight putative causal variants.
Here we can even include more samples (as the experimental cost will be lower anyway).
Therefore, we are now using high-depth phased data (in VCF format) for  180 samples from CHB, PEL and another Latin American population (CLM, from Colombia) to assess whether selection signatures are shared across other populations.

-------------------------- 

As we are using phased data, we are able to perform selection tests based on haplotype diversity and/or homozygosity, as seen during the lecture.
We are going to use the software [selscan](https://github.com/szpiech/selscan), which implements several selection tests (iHS, nSL, XP-EHH and so on).



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





