# Weggis

Bioinformatics for Adaptation Genomics - [Winter School 2016](http://www.adaptation.ethz.ch/education/winter-school-2016.html)

INFERENCE OF EVOLUTIONARY SIGNAL FROM SNP DATA

Part 1: Detection of signatures of selection

2nd March 2016

## Material

The data has been already downloaded and it is provided here `/gdc_home5/groups/bag2016/wednesday`.
To download and access all the material in this web page locally use [git](http://git-scm.com/).
Uncomment it before running it.
```
# git clone https://github.com/mfumagalli/Weggis.git
# cd Weggis
# git pull # to be sure you have the latest version, if so you should see "Already up-to-date."
```


## Data

As an illustration, we will use 80 BAM files of human samples (of African, European, East Asian, and Native American descent), a reference genome, and putative ancestral sequence.
We will also use VCF files for 180 individuals from the same populations.
The human data represents a small genomic region (1MB on chromosome 11) extracted from the 1000 Genomes Project data set.
More information on this project can be found [here](http://www.1000genomes.org/), including their last publication available [here](http://www.nature.com/nature/journal/v526/n7571/full/nature15393.html).

All data is publicly available.
A pipeline to retrieve such data is provided [here](https://github.com/mfumagalli/Weggis/blob/master/Files/data.sh).
You need to have 'samtools', 'bgzip' and 'Rscript' installed in your /usr/bin to run this.
Uncomment it before running it.
```
# bash Files/data.sh
```
`Data` and `Results` folder will be created automatically inside the `Weggis` folder.
Data will be saved (but not pushed to git main repository) in `Data` folder.
If you do not use the provided script, please make sure you create `Results` folder:
```
mkdir Results
```

Additional scripts are be provided in the `Scripts/` folder.

## Agenda

### Lecture

Slides are available in [pdf]() and [pptx]() format:
* The effect of selection on the genome
* Methods to detect selection signals
* The problem of assessing significance
* Summary statistics from low-depth data

### Practical

Genomic scan for selection from large-scale data set.
Case study: identification of allele frequency differentation between, with admixture assessment and quantification, from low-depth data: the case of FADS genetic variation in Native Americans

* [Rationale](https://github.com/mfumagalli/Weggis/blob/master/Files/rationale.md) and plan of action
* Brief introduction to [ANGSD](http://www.popgen.dk/angsd/index.php/Main_Page) and [population structure](https://github.com/mfumagalli/Weggis/blob/master/Files/explore.md) analysis
* Selection scan based on [genetic differentiation](https://github.com/mfumagalli/Weggis/blob/master/Files/selection.md) from low-depth data
* Assessing significance through [simulations](https://github.com/mfumagalli/Weggis/blob/master/Files/simulation.md)
* Selection test based on [haplotype diversity](https://github.com/mfumagalli/Weggis/blob/master/Files/deep.md) from high-depth data

#### Additional material

This material provides some guidelines on how to perform SNP and genotype calling and estimate population genetics metrics using [ANGSD](http://popgen.dk/wiki/index.php/ANGSD), suitable for low-depth data.
* Basic [filtering](https://github.com/mfumagalli/WoodsHole/blob/master/Files/filtering.md)
* Estimation of allele frequencies and [SNP calling](https://github.com/mfumagalli/WoodsHole/blob/master/Files/snpcall.md)
* [Genotype calling](https://github.com/mfumagalli/WoodsHole/blob/master/Files/genocall.md)
* Advanced methods to estimate [SFS](https://github.com/mfumagalli/WoodsHole/blob/master/Files/sfs.md)

## Credits

Some materials have been borrowed (and then adapted) from [Thorfinn Korneliussen](http://scholar.google.co.uk/citations?user=-YNWF4AAAAAJ&hl=en), [Anders Albrechtsen](http://popgen.dk/albrecht/web/WelcomePage.html), [Tyler Linderoth](http://scholar.google.com/citations?user=dTuxmzkAAAAJ&hl=en), [Filipe G. Vieira](http://scholar.google.com/citations?user=gvZmPNQAAAAJ&hl=en), [Dean Ousby](https://www.linkedin.com/in/deanousby), [Javier Mendoza](https://www.ucl.ac.uk/candela/candela-news/new-fellow-javiermendoza)  and possibly many more.


