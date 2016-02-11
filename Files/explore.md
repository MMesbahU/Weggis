

1) Investigate population structure (or clustering) of PEL samples (Peruvians), EUR (Europeans) and LWK (Africans).

One solution would be to perform a PCA, MDS or some clustering based on genetic distances among samples.
Then we can check whether some PEL fall within EUR/LWK or all PEL form a separate clade.

First, compute genotype posterior probabilities for all samples.

```
# Assuming HWE, without filtering based on probabilities, with SNP calling
angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-3 \
        -doGeno 8 -doPost 1 &> /dev/null
```

Record how many sites we retrieve.
```
N_SITES=`zcat Results/ALL.mafs.gz | tail -n+2 | wc -l`
echo $N_SITES
```

Create a file with labels indicating the population of interest for each sample.
```
Rscript -e 'cat(paste(rep(c("LWK","TSI","PEL"),each=20), rep(1:20, 3), sep="_"), sep="\n", file="pops.label")'
cat pops.label
```

With [ngsDist](https://github.com/fgvieira/ngsDist) we can computer pairwise genetic distances without relying on individual genotype calls.
```
ngsDist -verbose 1 -geno Results/ALL.geno.gz -probs -n_ind 60 -n_sites $N_SITES -labels pops.label -o Results/ALL.dist -n_threads 4
less -S Results/ALL.dist
```

We can visualise the pairwise genetic distances in form of a tree.
```
fastme -D 1 -i Results/ALL.dist -o Results/ALL.tree -m b -n b &> /dev/null
cat Results/ALL.tree
```

Plot the tree.
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




