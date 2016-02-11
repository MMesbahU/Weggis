
3) Estimate allele frequencies for SNPs in FADS genes of interest

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




