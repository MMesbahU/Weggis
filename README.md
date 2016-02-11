# Weggis

Bioinformatics for Adaptation Genomics - [Winter School 2016](http://www.adaptation.ethz.ch/education/winter-school-2016.html)

INFERENCE OF EVOLUTIONARY SIGNAL FROM SNP DATA

Part 1: Detection of signatures of selection

2nd March 2016

## Material

To download and access all the material in this web page use [git](http://git-scm.com/):
```
git clone https://github.com/mfumagalli/Weggis
cd Weggis
# git pull # to be sure you have the latest version, if so you should see "Already up-to-date."
```
`Data` and `Results` folder will be created automatically inside the `Weggis` folder.

## Data

As an illustration, we will use 60 BAM files of human samples (of African, European, and Native American descent), a reference genome, and putative ancestral sequence.
The human data represents a small genomic region (1MB on chromosome 11) extracted from the 1000 Genomes Project data set.
More information on this project can be found [here](http://www.1000genomes.org/), including their last publication available [here](http://www.nature.com/nature/journal/v526/n7571/full/nature15393.html).

All data is publicly available.
A pipeline to retrieve such data is provided [here](https://github.com/mfumagalli/Weggis/blob/master/Files/data.md).
You need to have 'samtools', 'bgzip' and 'Rscript' installed in your /usr/bin to run this.
```
bash Files/data.sh
```
Data will be saved (but not pushed to git main repository) in `Data` folder.

Additional scripts are be provided in the `Scripts/` folder.

## Agenda

### Lecture

* Session 1 - The effect of selection on the genome + Methods to detect selection signals (positive, negative, balancing)

* Session 2 - The problem of assessing significance + Critical discussion of case studies

### Practical

* Session 3 - Genomic scan for selection from large-scale data set

## Credits

Some materials have been borrowed (and then adapted) from [Thorfinn Korneliussen](http://scholar.google.co.uk/citations?user=-YNWF4AAAAAJ&hl=en), [Anders Albrechtsen](http://popgen.dk/albrecht/web/WelcomePage.html), [Tyler Linderoth](http://scholar.google.com/citations?user=dTuxmzkAAAAJ&hl=en), [Filipe G. Vieira](http://scholar.google.com/citations?user=gvZmPNQAAAAJ&hl=en), [Dean Ousby](https://www.linkedin.com/in/deanousby), and possibly many more.






