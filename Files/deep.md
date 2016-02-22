
We have some suggestive results pointing towards selection acting in the FADS region based on patterns of genetic differentiation against East Asians and Europeans.
However, we employed low-depth data and our analyses were not based on called genotypes.
Therefore, we did not have the power, for instance, to identify the causal variant (if any) and patterns of haplotype distribution.

It is usual to then perform a targeted deep resequencing of our region of interest, to further refione our selection analyses and highlight putative causal variants.
Here we can even include more samples (as the experimental cost will be lower anyway).
Therefore, we are now using high-depth phased data (in VCF format) for  180 samples from CHB, PEL and another Latin American population (CLM, from Colombia) to assess whether selection signatures are shared across other populations.

-------------------------- 

As we are using phased data, we are able to perform selection tests based on haplotype diversity and/or homozygosity, as seen during the lecture.
We are going to use the software [selscan](https://github.com/szpiech/selscan), which implements several selection tests (iHS, nSL, XP-EHH and so on).
According to its manual *[these] statistics are designed to use phased genotypes to identify putative regions of recent or ongoing positive selection in genomes. They are all based on the model of a hard selective sweep, where a de novo adaptive mutation arises on a haplotype that quickly sweeps toward fixation, reducing diversity around the locus. If selection is strong enough, this occurs faster than recombination or mutation can act to break up the haplotype, and thus a signal of high haplotype homozygosity can be observed extending from an adaptive locus.*

As an illustration, we are here calculating nSL.

First, specify the path to selscan and then look at the help page:
```
SS=/data/data/Software/selscan/bin/linux
$SS/selscan --help
```
The basic usage of selscan to compute nSL from VCF file is:
```
$SS/selscan --nsl --vcf Data/PEL.chr11.vcf --out Results/PEL
```
As you can see, unlike for iHS and EHH, we are not required to provide a genetic map since values are computed in windows based on physical distances.

For nSL, the integration is not cut after a certain threshold, but we need to set a value of maximum length for the window (in number of SNPs units) for building the haplotype.
This is set with the option `--max-extend-nsl`.
It is also common to filter out variant with very low frequency.

Therefore our command line might be:
```
$SS/selscan --nsl --vcf Data/PEL.chr11.vcf --out Results/PEL --max-extend-nsl 200 --maf 0.02
```
Have a look at the output file, knowning that the header is:
`< locusID > < physicalPos > < ’1 ’ freq > <sl1 > <sl0 > < unstandardized nSL >`
```
less -S Results/PEL.nsl.out
```

Finally, these unstandardized scores are normalized in allele frequency bins across the entire genome.
This can be achieved by the extra command `norm`:
```
$SS/norm --help
...
--bins <int>: The number of frequency bins in [0,1] for score normalization.
	Default: 100
...
```
Thus, our command would be:
```
$SS/norm --ihs --files Results/PEL.nsl.out --bins 20
```
The output file is called `Results/PEL.nsl.out.20bins.norm`.
```
less -S Results/PEL.nsl.out.20bins.norm
```
In theory, we should perform the normalization on genome-wide data (or at least chromosome-wide).

We can do the same for our other populations, CHB and CLM.
```
$SS/selscan --nsl --vcf Data/CHB.chr11.vcf --out Results/CHB --max-extend-nsl 200 --maf 0.02
$SS/norm --ihs --files Results/CHB.nsl.out --bins 20
$SS/selscan --nsl --vcf Data/CLM.chr11.vcf --out Results/CLM --max-extend-nsl 200 --maf 0.02
$SS/norm --ihs --files Results/CLM.nsl.out --bins 20
```

We can plot these results:
```
Rscript Scripts/plotnSL.R Results/PEL.nsl.out.20bins.norm Results/PEL.nsl.pdf
Rscript Scripts/plotnSL.R Results/CHB.nsl.out.20bins.norm Results/CHB.nsl.pdf
Rscript Scripts/plotnSL.R Results/CLM.nsl.out.20bins.norm Results/CLM.nsl.pdf
```

What is happening here? What is wrong? Why is there not a clear patterns of high/low scores?
Note that in these VCF alleles have been "polarise" according to the reference sequence, and not the ancestral one.

To overcome this issue, either you should use absolute values of nSL and perform a normalisation on the minor allele frequency rather than the non-reference allele frequency.

-----------------------

Another powerful statistic to delect selection from hard sweeps is XP-EHH, which measures the differential decay of haplotype homozygosity.
XP-EHH are not standardised based on allele frequency bins and therefore are not sensitive to misspecification of the ancestral state.
This statistic can be again computed using selscan, provided that we have a genetic map file.

A genetic map for chromosome 11 (based on 1000 Genomes data) is provided here:
```
less -S Files/genetic_map_GRCh37_chr11.map
```
However we need to extract only the sites that correspond in our VCF file.
A simple (slow) R script to do that is here:
```

```
Now we can run XP-EHH giving the resulting file as input.

```
$$SS/selscan --xp-ehh --vcf Data/PEL.chr11.vcf --vcf-ref Data/CHB.chr11.vcf --map Files/genetic.map --out Results/PEL
```
The output file has the header:
`< locusID > < physicalPos > < geneticPos > < popA ’1 ’ freq > < ihhA > < popB ’1 ’ freq > < ihhB > < unstandardized XPEHH >`
```
$SS/norm --xp-ehh --files Results/PEL.nsl.out --bins 20
```

Have a look at the results:
```

```
Plot them:
```

```

-----------------------------


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





