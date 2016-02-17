
Brief intro on ANGSD


--------------------

1) Investigate population structure our samples: PEL,(Peruvians), TSI (Europeans), LWK (Africans), and CHB (East Asians).

One solution would be to perform a Principal Component Analysis (PCA), Multidimensional Scaling (MDS) or some clustering based on genetic distances among samples.
Then we can check whether some PEL fall within EUR/LWK or all PEL form a separate clade.

To do this, we first need to assign genotypes (or their associated probabilities).
We now see how to use ANGSD to call genotypes.
The specific option is `-doGeno`:
```
$ANGSD/angsd -doGeno
...
-doGeno 0
        1: write major and minor
        2: write the called genotype encoded as -1,0,1,2, -1=not called
        4: write the called genotype directly: eg AA,AC etc
        8: write the posterior probability of all possible genotypes
        16: write the posterior probability of called gentype
        32: write the posterior probability of called gentype as binary
        -> A combination of the above can be choosen by summing the values, EG write 0,1,2 types with majorminor as -doGeno 3
        -postCutoff=0.333333 (Only genotype to missing if below this threshold)
        -geno_minDepth=-1       (-1 indicates no cutof)
        -geno_maxDepth=-1       (-1 indicates no cutof)
        -geno_minMM=-1.000000   (minimum fraction af major-minor bases)
        -minInd=0       (only keep sites if you call genotypes from this number of individuals)

        NB When writing the posterior the -postCutoff is not used
        NB geno_minDepth requires -doCounts
        NB geno_maxDepth requires -doCounts
```

Therefore, if we set `-doGeno 2`, genotypes are coded as 0,1,2, as the number of alternate alleles.
If we want to print the major and minor alleles as well then we set `-doGeno 3`.

To calculate the posterior probability of genotypes we need to define a model.
```
angsd -doPost

...
-doPost 0       (Calculate posterior prob 3xgprob)
        1: Using frequency as prior
        2: Using uniform prior
...
```
`-doPost 1` uses the estimate per-site allele frequency as a prior for genotype proportions, assuming Hardy Weinberg Equilibrium.
When the assumption of HWE is not valid, you can use an estimate of the inbreeding coefficient, for instance calculated using [ngsF](https://github.com/fgvieira/ngsF).

However, since our data is low-depth, genotypes cannot be assigned with high confidence and therefore we want to use **genotype posterior probabilities** instead, using options `-doGeno 8 -doPost 1`.

For more details and examples on genotype calling with ANGSD see [here](https://github.com/mfumagalli/WoodsHole/genocall.md)

--------------------------------

Furthermore, we want to restrict this analysis on a set of putative polymorphic sites (SNPs), as non-variable sites (across all samples) will not carry information regarding population structure or differentiation.

This search can be firstly translated in the estimation of the allele frequency for each site.
In other words, at each site we want to to estimate (or count) how many copies of different alleles (two in case of biallelic SNPs) we observe in our sample (across all sequenced individuals).

ANGSD has an option to estimate **allele frequencies**:

```
$ANGSD/angsd -doMaf
...
-doMaf  0 (Calculate persite frequencies '.mafs.gz')
        1: Frequency (fixed major and minor)
        2: Frequency (fixed major unknown minor)
        4: Frequency from genotype probabilities
        8: AlleleCounts based method (known major minor)
        NB. Filedumping is supressed if value is negative
-doPost 0       (Calculate posterior prob 3xgprob)
        1: Using frequency as prior
        2: Using uniform prior
Filters:
        -minMaf         -1.000000       (Remove sites with MAF below)
        -SNP_pval       1.000000        (Remove sites with a pvalue larger)
        -rmTriallelic   0.000000        (Remove sites with a pvalue lower)
Extras:
        -ref    (null)  (Filename for fasta reference)
        -anc    (null)  (Filename for fasta ancestral)
        -eps    0.001000 [Only used for -doMaf &8]
        -beagleProb     0 (Dump beagle style postprobs)
        -indFname       (null) (file containing individual inbreedcoeficients)
NB These frequency estimators requires major/minor -doMajorMinor
```

Therefore, the estimation of allele frequencies requires the specification of how to assign the major and minor alleles (if biallelic).
```
$ANGSD/angsd -doMajorMinor
...
        -doMajorMinor   0
        1: Infer major and minor from GL
        2: Infer major and minor from allele counts
        3: use major and minor from a file (requires -sites file.txt)
        4: Use reference allele as major (requires -ref)
        5: Use ancestral allele as major (requires -anc)
        -rmTrans: remove transitions 0
        -skipTriallelic 0
```

Finally, you need to specify which genotype likelihood model to use.
```
$ANGSD/angsd -GL
...
        -GL=0:
        1: SAMtools
        2: GATK
        3: SOAPsnp
        4: SYK
        5: phys
        -trim           0               (zero means no trimming)
        -tmpdir         angsd_tmpdir/   (used by SOAPsnp)
        -errors         (null)          (used by SYK)
        -minInd         0               (0 indicates no filtering)

Filedumping:
        -doGlf  0
        1: binary glf (10 log likes)    .glf.gz
        2: beagle likelihood file       .beagle.gz
        3: binary 3 times likelihood    .glf.gz
        4: text version (10 log likes)  .glf.gz
```

We may be interested in looking at allele frequencies only for sites that are actually variable in our sample.
Therefore we want to perform a **SNP calling**.
There are two main ways to call SNPs using ANGSD with these options:
```
        -minMaf         0.000000        (Remove sites with MAF below)
        -SNP_pval       1.000000        (Remove sites with a pvalue larger)
```
Therefore we can consider assigning as SNPs sites whose estimated allele frequency is above a certain threhsold (e.g. the frequency of a singleton) or whose probability of being variable is above a specified value.

--------------------------------

Back to our example, firstly, we need to compute genotype posterior probabilities for all samples with ANGSD only on putative polymorphic sites.

ANGSD can also perform some basic filtering of the data.
These filters are based on:

* quality and depth, see [here](http://www.popgen.dk/angsd/index.php/Filters)
* SNP quality, see [here](http://popgen.dk/angsd/index.php/SnpFilters)
* sites, see [here](http://popgen.dk/angsd/index.php/Sites)

For more details and examples on how to filtering data with ANGSD see [here](https://github.com/mfumagalli/WoodsHole/filtering.md)

As an illustration,...

```
# Assuming HWE, without filtering based on probabilities, with SNP calling
$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 40 -setMinDepth 40 -setMaxDepth 400 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-3 \
        -doGeno 8 -doPost 1 &> /dev/null
```

Record how many sites we retrieve.
```
N_SITES=`zcat Results/ALL.mafs.gz | tail -n+2 | wc -l`
echo $N_SITES
```

For plotting purposes, we now create a file with labels indicating the population of interest for each sample.
```
Rscript -e 'cat(paste(rep(c("LWK","TSI","CHB","PEL"),each=20), rep(1:20, 4), sep="_"), sep="\n", file="pops.label")'
cat pops.label
```

With [ngsDist](https://github.com/fgvieira/ngsDist) we can compute pairwise genetic distances without relying on individual genotype calls.
```
$NGSDIST/ngsDist -verbose 1 -geno Results/ALL.geno.gz -probs -n_ind 80 -n_sites $N_SITES -labels pops.label -o Results/ALL.dist -n_threads 4
less -S Results/ALL.dist
```

From these distances, we can perform a MDS analysis and investigate the population genetic structure of our samples.

N_SAMPLES=80
tail -n +3 Data/ALL.dist | head -n $N_SAMPLES | Rscript --vanilla --slave Scripts/get_MDS.R --no_header --data_symm -n 4 -m "mds" -o Data/ALL.mds


We can visualise the pairwise genetic distances in form of a tree (in Newick format).
```
$FASTME -D 1 -i Results/ALL.dist -o Results/ALL.tree -m b -n b &> /dev/null
cat Results/ALL.tree
```

We can some R packages to plot the resulting tree.
```
Rscript -e 'library(ape); library(phangorn); pdf(file="Results/ALL.tree.pdf"); plot(read.tree("Results/ALL.tree"), cex=0.5); dev.off();' &> /dev/null
evince Results/ALL.tree.pdf
```

Therefore, PEL samples appear very close to EUR.
The next step would be to compute admixture proportions in order to select a subset of putative Native American (unadmixed) samples.

------------------------
2) Compute admixture proportions across samples.

We use ngsAdmix, which again works on genotype probabilities and not on individual calls.
This is suitable for low-depth data.

ngsAdmix requires genotype likelihoods in BEAGLE format as input.
We can compute these quantities with ANGSD with `-doGlf 2`.
```
angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-3 \
        -doGlf 2 &> /dev/null
```

We assume 3 ancestral populations (Europeans, Africans and Native Americans) making up the genetic diversity of our samples.
Therefore we compute admixture proportions with 3 ancestral components.
```
K=3
NGSadmix -likes Results/ALL.beagle.gz -K $K -outfiles Results/ALL.admix.K$K -P 4 -minMaf 0.02
```

Combine samples IDs with admixture proportions and inspect the results.
```
paste ALL.bamlist Results/ALL.admix.K$K.qopt > Results/ALL.admix.K$K.txt
less -S Results/ALL.admix.K$K.txt
```

From these quantities we can extract how many samples (and which ones) have a high proportion of Native American ancestry (e.g. >0.90).
We can also plot the individual ancestral proportions for PEL samples.
We need to specify the order of ancestral populations in the admixture file `Results/ALL.admix.K$K.txt`.
In my case these are PEL TSI LWK
```
Rscript Scripts/getUnadmixed.R 0.90 PEL TSI LWK
```
Inspect the results.
```
less -S Results/PEL_unadm.BAMs.txt
evince Results/ALL.admix.PEL.pdf
```

Now we have a subset of putative Native American samples.
We can calculate allele frequencies for only these samples.

------------------------

[HOME](https://github.com/mfumagalli/Weggis)




