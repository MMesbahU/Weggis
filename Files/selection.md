

### Allele frequency differentiation

The joint-SFS can be considered as a summary statistics for the level of genetic differentiation between populations.
We have also seen how it can be used as prior information when computing FST (and related metrics) without relying on genotype calling.

Here we see how to compute the 2D-SFS, FST, PBS, and other summary statistics from low-depth data using ANGSD.
Our final goal is to detect signatures of selection in our data.

-------------------------------

One of the most important aspect of data analysis for population genetics is the estimate of the Site Frequency Spectrum (SFS). 
SFS records the proportions of sites at different allele frequencies. It can be folded or unfolded, and the latter case implies the use of an outgroup species to define the ancestral state. 
SFS is informative on the demography of the population or on selective events (when estimated at a local scale).

We use ANGSD to estimate SFS using on example dataset, using the methods described [here](http://www.ncbi.nlm.nih.gov/pubmed/22911679).
Details on the implementation can be found [here](http://popgen.dk/angsd/index.php/SFS_Estimation).
Briefly, from sequencing data one computes genotype likelihoods (as previously described). 
From these quantities ANGSD computes posterior probabilities of Sample Allele Frequency (SAF), for each site. 
Finally, an estimate of the SFS is computed.

These steps can be accomplished in ANGSD using `-doSaf 1/2` options and the program `realSFS`.

```
$ANGSD/angsd -doSaf
...
-doSaf		0
	1: perform multisample GL estimation
	2: use an inbreeding version
	3: calculate genotype probabilities
	4: Assume genotype posteriors as input (still beta) 
	-doThetas		0 (calculate thetas)
	-underFlowProtect	0
	-fold			0 (deprecated)
	-anc			(null) (ancestral fasta)
	-noTrans		0 (remove transitions)
	-pest			(null) (prior SFS)
	-isHap			0 (is haploid beta!)
NB:
	  If -pest is supplied in addition to -doSaf then the output will then be posterior probability of the sample allelefrequency for each site
```

The SFS is typically computed for each population separately.
We need to slightly modify the filtering options as now each population has 20 samples. 
So now we set `-minInd 10 -setMinDepth 10 -setMaxDepth 100`.
Also, we want to estimate the unfolded SFS and we use a putative ancestral sequence to polarise our alleles (to ancestral and derived states).

```
for POP in LWK TSI CHB PEL
do
        echo $POP
        $ANGSD/angsd -P 4 -b $POP.bamlist -ref $REF -anc $ANC -out Results/$POP \
		-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
		-minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 10 -setMaxDepth 200 -doCounts 1 \
		-GL 1 -doSaf 1 &> /dev/null
done
```

Have a look at the output file.
```
$ANGSD/misc/realSFS print Results/PEL.saf.idx | less -S
```
These values represent the sample allele frequency likelihoods at each site, as seen during the lecture.

The next step would be to use these likelihoods and estimate the overall SFS.
This is achieved by the program `realSFS`.

```
$ANGSD/misc/realSFS
	-> ---./realSFS------
	-> EXAMPLES FOR ESTIMATING THE (MULTI) SFS:

	-> Estimate the SFS for entire genome??
	-> ./realSFS afile.saf.idx 

	-> 1) Estimate the SFS for entire chromosome 22 ??
	-> ./realSFS afile.saf.idx -r chr22 

	-> 2) Estimate the 2d-SFS for entire chromosome 22 ??
	-> ./realSFS afile1.saf.idx  afile2.saf.idx -r chr22 

	-> 3) Estimate the SFS for the first 500megabases (this will span multiple chromosomes) ??
	-> ./realSFS afile.saf.idx -nSites 500000000 

	-> 4) Estimate the SFS around a gene ??
	-> ./realSFS afile.saf.idx -r chr2:135000000-140000000 

	-> Other options [-P nthreads -tole tolerence_for_breaking_EM -maxIter max_nr_iterations -bootstrap number_of_replications]

	-> See realSFS print for possible print options
	-> Use realSFS print_header for printing the header

	->------------------
	-> NB: Output is now counts of sites instead of log probs!!
	-> NB: You can print data with ./realSFS print afile.saf.idx !!
	-> NB: Higher order SFS can be estimated by simply supplying multiple .saf.idx files!!
	-> NB: Program uses accelerated EM, to use standard EM supply -m 0
```

This command will estimate the SFS for each population:
```
for POP in LWK TSI CHB PEL
do
        echo $POP
        $ANGSD/misc/realSFS Results/$POP.saf.idx 2> /dev/null > Results/$POP.sfs
done
```

Have a look at the output file:
```
cat Results/PEL.sfs
```
How many values do you expect?
```
awk -F' ' '{print NF; exit}' Results/PEL.sfs 
```

OPTIONAL

Let us plot the SFS for each pop using this simple R script.
```
Rscript Scripts/plotSFS.R Results/LWK.sfs Results/TSI.sfs Results/CHB.sfs Results/PEL.sfs
evince Results/LWK_TSI_CHB_PEL.pdf
```
Do they behave like expected? 
Which population has more SNPs?
Which population has a higher proportion of common (not rare) variants?
Why is it bumpy?

It is sometimes convenient to generate bootstrapped replicates of the SFS, by sampling with replacements genomic segments.
This could be used for instance to get confidence intervals when using the SFS for demographic inferences (as you may seen in the next days).
This can be achieved in ANGSD using:
```
$ANGSD/misc/realSFS Results/PEL.saf.idx -bootstrap 10  2> /dev/null > Results/PEL.boots.sfs
cat Results/PEL.boots.sfs
```
This command may take some time.
The output file has one line for each boostrapped replicate.

More examples on how to estimate the SFS with ANGSD can be found [here](https://github.com/mfumagalli/WoodsHole/blob/master/Files/sfs.md).

---------------------------------------

OPTIONAL

It is very useful to estimate a **multi-dimensional SFS**, for instance the joint SFS between 2 populations (2D).
This can be used for making inferences on their divergence process (time, migration rate and so on).
However, here we are interested in estimating the 2D-SFS as prior information for our FST/PBS.

An important issue when doing this is to be sure that we are comparing the exactly same corresponding sites between populations.
ANGSD does that automatically and considers only a set of overlapping sites.
The 2D-SFS between all populations and PEL are computed with:
```
for POP in LWK TSI CHB
do
        echo $POP
	$ANGSD/misc/realSFS -P 4 Results/$POP.saf.idx Results/PEL.saf.idx 2> /dev/null > Results/$POP.PEL.sfs
done
```

The output file is a flatten matrix, where each value is the count of sites with the corresponding joint frequency ordered as [0,0] [0,1] and so on.
```
less -S Results/LWK.PEL.sfs
```
You can plot it, but you need to define how many samples you have per population.
```
Rscript Scripts/plot2DSFS.R Results/LWK.PEL.sfs 20 20
evince Results/LWK.PEL.sfs.pdf
```

You can even estimate SFS with higher order of magnitude.
This command may take some time.
```
$ANGSD/misc/realSFS -P 4 Results/LWK.saf.idx Results/TSI.saf.idx Results/PEL.saf.idx 2> /dev/null > Results/LWK.TSI.PEL.sfs
```

------------------------------------

Here we are going to calculate allele frequency differentiation using the PBS (population branch statistic) metric.
Again, we can achieve this by avoid genotype calling using ANGSD.
From the sample allele frequencies likelihoods (.saf files) we can estimate PBS using the following pipeline.

Note that here we use the previously calculated SFS as prior information.
Also, PEL is our target population, while CHB and TSI are reference populations.
Therefore, we need to compute also the 2D-SFS between TSI and CHB and PEL (if not already done):
```
for POP in TSI CHB
do
        $ANGSD/misc/realSFS -P 4 Results/$POP.saf.idx Results/PEL.saf.idx 2> /dev/null > Results/$POP.PEL.sfs
done
$ANGSD/misc/realSFS -P 4 Results/TSI.saf.idx Results/CHB.saf.idx 2> /dev/null > Results/TSI.CHB.sfs
```

PBS can be calculated (in windows) using the following commands:
```
$ANGSD/misc/realSFS fst index Results/TSI.saf.idx Results/CHB.saf.idx Results/PEL.saf.idx -sfs Results/TSI.CHB.sfs -sfs Results/TSI.PEL.sfs -sfs Results/CHB.PEL.sfs -fstout Results/PEL.pbs
# the nex command will perform a sliding-window analysis
$ANGSD/misc/realSFS fst stats2 Results/PEL.pbs.idx -win 50000 -step 10000 > Results/PEL.pbs.txt
```

Let us have a look at the output file:
./realSFS fst print


SUMMARY STATS in PEL

$ANGSD/angsd -b $BAM/PEL_unadm.txt -anc $ANC -minInd 10 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minQ 20 -minMapQ 20 -out Data/PEL -GL 1 -doMajorMinor 5 -doThetas 1 -pest Data/PEL.ml -doSaf 1
$ANGSD/misc/thetaStat make_bed Data/PEL.thetas.gz
$ANGSD/misc/thetaStat do_stat Data/PEL.thetas.gz -nChr 1 -win 20000 -step 2000 -outnames Data/PEL.thetas

plot


OPTIONAL

FST - this can be in the section when I do simulations

$ANGSD/misc/realSFS fst index Data/TSI.saf.idx Data/CHB.saf.idx -sfs Data/TSI.CHB.ml -fstout Data/TSI.CHB

$ANGSD/misc/realSFS fst stats2 Data/TSI.CHB.fst.idx -win 20000 -step 2000 > Data/TSI.CHB.fst.txt

./realSFS fst print 

cut -f 2- Data/LWK.PEL.fst.txt | cut -c 2- > Data/LWK.PEL.fst.slidwindtxt

# Rscript Scripts/plotFST.R

also plot genes below as in locuszoom?

..........................................................

3) Estimate allele frequencies for SNPs in FADS genes of interest

maybe move this after simulations? as it is more on haplotype based analysis?

In ANGSD we can restrict our analyses on a subset of positions of interest using the `-sites` option.
The positions we are looking at are the one found under selection in Inuit, shown [here](https://github.com/mfumagalli/WoodsHole/blob/master/Files/snps_inuit.png):
- 11 61627960 <br>
- 11 61631510 <br>
- 11 61632310 <br>
- 11 61641717 <br>
- 11 61624414 <br>
- 11 61597212 <br>

The file with these positions need to be formatted as (chromosome positions).
```
> snps.txt
echo 11 61627960 >> snps.txt
echo 11 61631510 >> snps.txt
echo 11 61632310 >> snps.txt
echo 11 61641717 >> snps.txt
echo 11 61624414 >> snps.txt
echo 11 61597212 >> snps.txt
```
Inspect the file.
```
cat snps.txt
```

We need to index this file in order for ANGSD to process it.
```
angsd sites index snps.txt
```
We are interested in calculating the derived allele frequencies, so are using the ancestral sequence to polarise the alleles.
Create new lists of BAM files.
```
head -n 20 ALL.bamlist > LWK.sub.bamlist
tail -n 20 ALL.bamlist > TSI.sub.bamlist
cp Results/PEL_unadm.BAMs.txt PEL.sub.bamlist
```

Run ANGSD to compute allele frequencies.
Here we change the filtering (more relaxed) since we are interested in outputting all sites.
```
for POP in LWK TSI PEL
do
        echo $POP
        angsd -P 4 -b $POP.sub.bamlist -ref $REF -anc $ANC -out Results/$POP \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 1 -setMinDepth 10 -setMaxDepth 500 -doCounts 1 \
                -GL 1 -doMajorMinor 5 -doMaf 1 -skipTriallelic 1 \
                -sites snps.txt &> /dev/null
done
```

Inspect the results.
```
zcat Results/LWK.mafs.gz Results/TSI.mafs.gz Results/PEL.mafs.gz
```

Do you see any allele frequency differentiation?

----------------------------

4) Sliding windows scan for PBS


------------------------

[HOME](https://github.com/mfumagalli/Weggis)




