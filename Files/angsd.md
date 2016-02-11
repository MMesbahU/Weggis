
# Brief tutorial on ANGSD

input files \\
how to build up a command line \\
some examples on geno/snp calling \\
intro on summary stats \\

MAPPED DATA > FILTERING

In this section, we will show how to perform a basic filtering of sites, after the reads have been mapped or aligned.

For most of the examples, we will use the program [ANGSD](http://popgen.dk/wiki/index.php/ANGSD) (Analysis of Next Generation Sequencing Data) developed by Thorfinn Korneliussen and Anders Albrechtsen at the University of Copenhagen.
More information about its rationale and implemented methods can be found [here](http://www.ncbi.nlm.nih.gov/pubmed/25420514).

We will use 60 BAM files of human samples (of African, European, and Native American descent), a reference genome, and putative ancestral sequence.
The human data represents a small genomic region (1MB on chromosome 11) extracted from the 1000 Genomes Project data set.

## Preparation

Please set the path for all programs and data we will be using.
As an example these are my paths.
```
ANGSD=/data/data/Software/angsd
SAMTOOLS=/data/data/Software/samtools-1.3
NGSDIST=/data/Software/ngsDist
NGSTOOLS=/data/Software/ngsTools
NGSADMIX=/data/data/Software/NGSadmix/NGSadmix
FASTME=/data/data/Software/fastme-2.1.4/src/fastme
```
However, these paths have been sym-linked to your /usr/bin so they can be called by simply typing their name, e.g. `angsd`.

If you downloaded the data using the provided script, this is what you should specify.
```
REF=Data/hs37d5.fa.gz
ANC=Data/hg19ancNoChr.fa.gz
```

#### ANGSD

Here we will use ANGSD to analyse our data. To see a full list of options type
```
angsd
```
and you should see something like
```
...
Overview of methods:
        -GL             Estimate genotype likelihoods
        -doCounts       Calculate various counts statistics
        -doAsso         Perform association study
        -doMaf          Estimate allele frequencies
        -doError        Estimate the type specific error rates
        -doAncError     Estimate the errorrate based on perfect fastas
        -doHWE          Est inbreedning per site
        -doGeno         Call genotypes
        -doFasta        Generate a fasta for a BAM file
        -doAbbababa     Perform an ABBA-BABA test
        -sites          Analyse specific sites (can force major/minor)
        -doSaf          Estimate the SFS and/or neutrality tests genotype calling
        -doHetPlas      Estimate hetplasmy by calculating a pooled haploid frequency

        Below are options that can be usefull
        -bam            Options relating to bam reading
        -doMajorMinor   Infer the major/minor using different approaches
        -ref/-anc       Read reference or ancestral genome
        -doSNPstat      Calculate various SNPstat
        many others

For information of specific options type:
        ./angsd METHODNAME eg
                ./angsd -GL
                ./angsd -doMaf
                ./angsd -doAsso etc
                ./angsd sites for information about indexing -sites files
Examples:
        Estimate MAF for bam files in 'list'
                './angsd -bam list -GL 2 -doMaf 2 -out RES -doMajorMinor 1'
```


We will see later of to perform SNP and genotype calling (and many other things) with ANGSD.

ANGSD can accept several input files, as described [here](http://popgen.dk/angsd/index.php/Input):

* BAM/CRAM
* Pileup
* Genotype likelihood/probability files
* VCF

#### Basic filtering post-mapping

Here we show how ANGSD can also perform some basic filtering of the data.
These filters are based on:

* quality and depth, see [here](http://www.popgen.dk/angsd/index.php/Filters)
* SNP quality, see [here](http://popgen.dk/angsd/index.php/SnpFilters)
* sites, see [here](http://popgen.dk/angsd/index.php/Sites)

If the input file is in BAM format, the possible options are:
```
angsd -bam
...
---------------
parseArgs_bambi.cpp: bam reader:
        -bam/-b         (null)  (list of BAM/CRAM files)
        -i              (null)  (Single BAM/CRAM file)
        -r              (null)  Supply a single region in commandline (see examples below)
        -rf             (null)  Supply multiple regions in a file (see examples below)
        -remove_bads    1       Discard 'bad' reads, (flag >=256)
        -uniqueOnly     0       Discards reads that doesnt map uniquely
        -show           0       Mimic 'samtools mpileup' also supply -ref fasta for printing reference column
        -minMapQ        0       Discard reads with mapping quality below
        -minQ           13      Discard bases with base quality below
        -trim           0       Number of based to discard at both ends of the reads
        -only_proper_pairs 1    Only use reads where the mate could be mapped
        -C              0       adjust mapQ for excessive mismatches (as SAMtools), supply -ref
        -baq            0       adjust qscores around indels (as SAMtools), supply -ref
        -checkBamHeaders 1      Exit if difference in BAM headers
        -doCheck        1       Keep going even if datafile is not suffixed with .bam/.cram
        -downSample     0.000000        Downsample to the fraction of original data
        -nReads         50      Number of reads to pop from each BAM/CRAMs
        -minChunkSize   250     Minimum size of chunk sent to analyses

Examples for region specification:
                chr:            Use entire chromosome: chr
                chr:start-      Use region from start to end of chr
                chr:-stop       Use region from beginning of chromosome: chr to stop
                chr:start-stop  Use region from start to stop from chromosome: chr
                chr:site        Use single site on chromosome: chr
```
Some basic filtering consists in removing, for instance, read with low quality and/or with multiple hits, and this can be achieved using the parameters ```-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1```.

Also, you may want to remove reads with low mapping quality and sites with low quality or covered by few reads (low depth).
Under these circumnstances, the assignment of individual genotypes and SNPs is problematics, and can lead to errors.

For more details and examples on how to filtering data with ANGSD see [here](https://github.com/mfumagalli/WoodsHole/filtering.md)

-----------------------------------------

**WORKFLOW**:

MAPPED DATA > FILTERING > SNP CALLING

In this section, we will go through some examples on how to assign variable sites from BAM files, once the data has been filtered.
We will show how to call SNPs with different methods, and we will compare their results.

### Estimating allele frequencies and calling SNPs

One of the first quantities of interest is the identification of which sites are variable, or polymorphic, in our sample.
This search can be firstly translated in the estimation of the allele frequency for each site.
In other words, at each site we want to to estimate (or count) how many copies of different alleles (two in case of biallelic SNPs) we observe in our sample (across all sequenced individuals).

ANGSD has an option to estimate **allele frequencies**:

```
angsd -doMaf
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
angsd -doMajorMinor
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

Finally, one needs to specify with genotype likelihood model to use.
```
angsd -GL
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
A description of these different implementation can be found [here](http://www.popgen.dk/angsd/index.php/Genotype_likelihoods).
The GATK model refers to the first GATK paper (as seen in the lecture), SAMtools is somehow more sophisticated (non-independence of errors), SOAPsnp requires a reference sequence for recalibration of quality scores, SYK is error-type specific.
For most applications and data, GATK and SAMtools models should give similar results and be preferred.

From these observations, our command line could be:
```
angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 -sites sites.txt\
        -GL 1 -doMajorMinor 4 -doMaf 1 -skipTriallelic 1 &> /dev/null
```
where we specify:
* -GL 1: genotype likelihood model is in SAMtools
* -doMajorMinor 4: force the major allele to be the reference (the minor is inferred)
* -doMaf 1: major and minor are fixed

Which are the output files?
```
->"Results/ALL.arg"
->"Results/ALL.mafs.gz"
```
`.args` file is a summary of all options used, while `.mafs.gz` file shows the allele frequencies computed at each site.

Have a look at this file which contains estimates of **Minor Allele Frequency** (MAF) values.
```
zcat Results/ALL.mafs.gz | head
```
and you should see
```
chromo  position        major   minor   ref     knownEM nInd
11      61005556        C       A       C       0.000001        49
11      61005584        C       A       C       0.000001        49
11      61005597        C       T       C       0.001915        49
11      61005994        A       C       A       0.000001        46
11      61005995        G       A       G       0.000001        46
11      61005996        C       A       C       0.000001        46
11      61005997        C       A       C       0.000001        46
11      61005998        T       A       T       0.000001        46
11      61005999        G       A       G       0.000000        46
```

The first and second column indicate the position of each site, then we have major and minor alleles (based on the reference sequence), the estimated allele frequency, and the number of samples with data.

------------

We may be interested in looking at allele frequencies only for sites that are actually variable in our sample.
Therefore we want to perform a **SNP calling**.
There are several ways to call SNPs using ANGSD, for instance by using these options:
```
        -minMaf         0.000000        (Remove sites with MAF below)
        -SNP_pval       1.000000        (Remove sites with a pvalue larger)
```
Therefore we can consider assigning as SNPs sites whose estimated allele frequency is above a certain threhsold (e.g. the frequency of a singleton) or whose probability of being variable is above a specified value.

As an illustration, let us call SNPs by computing:
 - genotype likelihoods using GATK method;
 - major and minor alleles from allele counts (you need to specify -doCounts 1; not recommended in many cases);
 - frequency from known major allele but unknown minor;
 - SNPs as those having MAF=>0.01.

```
angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 -sites sites.txt \
        -GL 2 -doMajorMinor 2 -doMaf 2 -skipTriallelic 1  \
        -minMaf 0.01 &> /dev/null
```

You can have a look at the results:
```
zcat Results/ALL.mafs.gz | head

chromo  position        major   minor   ref     unknownEM       nInd
11      61006040        G       A       G       0.010250        53
11      61007218        T       C       C       0.444894        58
11      61007710        C       G       C       0.303258        59
11      61007781        A       G       G       0.010488        59
11      61007786        G       T       T       0.025118        59
11      61007804        T       C       C       0.329408        55
11      61014040        C       A       C       0.051804        60
11      61014194        C       T       C       0.014751        59
11      61014749        A       T       A       0.010908        60
```

For more details and examples on SNP calling with ANGSD see [here](https://github.com/mfumagalli/WoodsHole/snpcall.md)

-----------------------

Geno call?

**WORKFLOW**:

FILTERED DATA > SNP CALLING > GENOTYPE CALLING

Here we will explore several ways to call genotypes from sequencing data.
We will use ANGSD (and SAMtools as additional material), and compare results using different options.

We now see how to use ANGSD to call genotypes.
The specific option is `-doGeno`:
```
angsd -doGeno
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
We will see later what to do when the assumption of HWE is not valid.

A typical command for genotype calling assuming HWE is:

```
angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 30 -setMinDepth 210 -setMaxDepth 700 -doCounts 1 -sites sites.txt\
        -GL 1 -doMajorMinor 1 -doMaf 2 -skipTriallelic 1 \
        -SNP_pval 1e-3 \
        -doGeno 3 -doPost 1 -postCutoff 0 &> /dev/null
```

Have a look at the output file:
```
less -S Results/ALL.geno.gz
```

For more details and examples on genotype calling with ANGSD see [here](https://github.com/mfumagalli/WoodsHole/genocall.md)

------------------------

[HOME](https://github.com/mfumagalli/Weggis)




