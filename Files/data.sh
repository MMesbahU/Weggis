
## Pipeline to download and process the data to be used for this workshop.

# set path
# SAMTOOLS=/data/data/Software/samtools-1.3/samtools
SAMTOOLS=samtools
echo Is this your path to samtools? $SAMTOOLS

mkdir Data
mkdir Results

# get unrelated samples iDs
Rscript Scripts/getUnrelated.R

# download BAM files
bash Scripts/getBams.sh
# this creates files and folders in Data/PEL.BAMs/* and TSI and LWK and CHB

# create file with list of BAMs
ls Data/LWK.BAMs/*.bam Data/TSI.BAMs/*.bam Data/CHB.BAMs/*.bam Data/PEL.BAMs/*.bam > ALL.bamlist
ls Data/LWK.BAMs/*.bam Data/TSI.BAMs/*.bam Data/PEL.BAMs/*.bam > ALL_noCHB.bamlist
ls Data/LWK.BAMs/*.bam > LWK.bamlist
ls Data/TSI.BAMs/*.bam > TSI.bamlist
ls Data/CHB.BAMs/*.bam > CHB.bamlist
ls Data/PEL.BAMs/*.bam > PEL.bamlist

# download ancestral sequence
echo Downloading and processing ancestral sequence...
wget http://dna.ku.dk/~thorfinn/hg19ancNoChr.fa.gz &>/dev/null
gunzip hg19ancNoChr.fa.gz &>/dev/null
$SAMTOOLS faidx hg19ancNoChr.fa
mv hg19ancNoChr.* Data/.

# download reference sequence
echo Downloading and processing reference sequence...
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz &>/dev/null
gunzip hs37d5.fa.gz &>/dev/null
$SAMTOOLS faidx hs37d5.fa
mv hs37d5.* Data/.

# download VCF files
echo Downloading VCF file...
wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz &> /dev/null
wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi &> /dev/null
echo Processing VCF file...
VCFLIB=/data/data/Software/vcflib/bin
echo Is this your path to VCFlib? $VCFLIB
VCF=ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# whole region for selscan
$VCFLIB/vcffilter -f "VT = SNP" -f "AC > 10" -r 11:61000000-62000000 $VCF > ALL.chr11.tmp.vcf
grep -v MULTI_ALLELIC ALL.chr11.tmp.vcf > Data/ALL.chr11.vcf

$VCFLIB/vcfkeepsamples Data/ALL.chr11.vcf HG01565 HG01566 HG01571 HG01572 HG01577 HG01578 HG01892 HG01893 HG01917 HG01918 HG01920 HG01921 HG01923 HG01924 HG01926 HG01927 HG01932 HG01933 HG01935 HG01938 HG01939 HG01941 HG01942 HG01944 HG01945 HG01947 HG01948 HG01950 HG01951 HG01953 HG01954 HG01961 HG01965 HG01967 HG01968 HG01970 HG01971 HG01973 HG01974 HG01976 HG01977 HG01979 HG01980 HG01982 HG01991 HG01992 HG01995 HG01997 HG02002 HG02003 HG02008 HG02089 HG02090 HG02104 HG02105 HG02146 HG02147 HG02252 HG02253 HG02259 > Data/PEL.chr11.vcf
sed -i -e 's/\t\./\t0|0/g' Data/PEL.chr11.vcf

$VCFLIB/vcfkeepsamples Data/ALL.chr11.vcf NA18524 NA18525 NA18526 NA18527 NA18528 NA18529 NA18530 NA18531 NA18532 NA18533 NA18534 NA18535 NA18536 NA18537 NA18538 NA18539 NA18540 NA18541 NA18542 NA18543 NA18544 NA18545 NA18546 NA18547 NA18548 NA18549 NA18550 NA18552 NA18553 NA18555 NA18557 NA18558 NA18559 NA18560 NA18561 NA18562 NA18563 NA18564 NA18565 NA18566 NA18567 NA18569 NA18570 NA18571 NA18572 NA18573 NA18574 NA18575 NA18576 NA18577 NA18579 NA18580 NA18582 NA18583 NA18591 NA18592 NA18593 NA18595 NA18596 NA18597 > Data/CHB.chr11.vcf
sed -i -e 's/\t\./\t0|0/g' Data/CHB.chr11.vcf

$VCFLIB/vcfkeepsamples Data/ALL.chr11.vcf HG01112 HG01113 HG01119 HG01121 HG01122 HG01124 HG01125 HG01130 HG01131 HG01133 HG01134 HG01136 HG01137 HG01139 HG01140 HG01142 HG01148 HG01149 HG01250 HG01251 HG01253 HG01254 HG01256 HG01257 HG01259 HG01260 HG01269 HG01271 HG01272 HG01274 HG01275 HG01277 HG01278 HG01280 HG01281 HG01284 HG01341 HG01342 HG01344 HG01345 HG01347 HG01348 HG01350 HG01351 HG01353 HG01354 HG01356 HG01357 HG01359 HG01360 HG01362 HG01363 HG01365 HG01366 HG01369 HG01372 HG01374 HG01375 HG01377 HG01378 > Data/CLM.chr11.vcf
sed -i -e 's/\t\./\t0|0/g' Data/CLM.chr11.vcf

# fads region 11:61567097-61659006
# fads2
$VCFLIB/vcffilter -f "VT = SNP" -f "AC > 10" -r 11:61595713-61659006 $VCF > ALL.fads.tmp.vcf
grep -v MULTI_ALLELIC ALL.fads.tmp.vcf > Data/ALL.fads.vcf

$VCFLIB/vcfkeepsamples Data/ALL.fads.vcf HG01565 HG01566 HG01571 HG01572 HG01577 HG01578 HG01892 HG01893 HG01917 HG01918 HG01920 HG01921 HG01923 HG01924 HG01926 HG01927 HG01932 HG01933 HG01935 HG01938 HG01939 HG01941 HG01942 HG01944 HG01945 HG01947 HG01948 HG01950 HG01951 HG01953 HG01954 HG01961 HG01965 HG01967 HG01968 HG01970 HG01971 HG01973 HG01974 HG01976 HG01977 HG01979 HG01980 HG01982 HG01991 HG01992 HG01995 HG01997 HG02002 HG02003 HG02008 HG02089 HG02090 HG02104 HG02105 HG02146 HG02147 HG02252 HG02253 HG02259 > Data/PEL.fads.vcf
sed -i -e 's/\t\./\t0|0/g' Data/PEL.fads.vcf

$VCFLIB/vcfkeepsamples Data/ALL.fads.vcf NA18524 NA18525 NA18526 NA18527 NA18528 NA18529 NA18530 NA18531 NA18532 NA18533 NA18534 NA18535 NA18536 NA18537 NA18538 NA18539 NA18540 NA18541 NA18542 NA18543 NA18544 NA18545 NA18546 NA18547 NA18548 NA18549 NA18550 NA18552 NA18553 NA18555 NA18557 NA18558 NA18559 NA18560 NA18561 NA18562 NA18563 NA18564 NA18565 NA18566 NA18567 NA18569 NA18570 NA18571 NA18572 NA18573 NA18574 NA18575 NA18576 NA18577 NA18579 NA18580 NA18582 NA18583 NA18591 NA18592 NA18593 NA18595 NA18596 NA18597 > Data/CHB.fads.vcf
sed -i -e 's/\t\./\t0|0/g' Data/CHB.fads.vcf

$VCFLIB/vcfkeepsamples Data/ALL.fads.vcf HG01112 HG01113 HG01119 HG01121 HG01122 HG01124 HG01125 HG01130 HG01131 HG01133 HG01134 HG01136 HG01137 HG01139 HG01140 HG01142 HG01148 HG01149 HG01250 HG01251 HG01253 HG01254 HG01256 HG01257 HG01259 HG01260 HG01269 HG01271 HG01272 HG01274 HG01275 HG01277 HG01278 HG01280 HG01281 HG01284 HG01341 HG01342 HG01344 HG01345 HG01347 HG01348 HG01350 HG01351 HG01353 HG01354 HG01356 HG01357 HG01359 HG01360 HG01362 HG01363 HG01365 HG01366 HG01369 HG01372 HG01374 HG01375 HG01377 HG01378 > Data/CLM.fads.vcf
sed -i -e 's/\t\./\t0|0/g' Data/CLM.fads.vcf


# clean up
rm $VCF $VCF.tbi ALL.chr11.tmp.vcf ALL.fads.tmp.vcf


wget http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz .

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/CHB_omni_recombination_20130507.tar

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/CLM_omni_recombination_20130507.tar

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/PEL_omni_recombination_20130507.tar


echo Done!
ls -lh Data/* > Data/download.log
echo Open Data/download.log to see what files have been generated.





