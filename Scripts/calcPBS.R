
# not a general way to do it, it assumes known sample sizes and their order!

# calculate FST and PBS from called genotypes from ANGSD

args=commandArgs(T)
fin=args[1]
fout=args[2]
rm(args)

source("/gdc_home5/groups/bag2016/wednesday/Scripts/popgen.R")

# fin="DUMB.geno.gz"

genos=read.table(fin, head=F, stringsAsFact=F)

haplos=list()
hap1=hap2=hap3=c()
for (i in 1:20) hap1[i]=paste(genos[,(i+22)], sep="",collapse="")
for (i in 1:20) hap2[i]=paste(genos[,(i+42)], sep="",collapse="")
for (i in 1:20) hap3[i]=paste(genos[,(i+62)], sep="",collapse="")
haplos=list(hap1,hap2,hap3)

lwin=50000
lstep=10000
pos=genos[,2]

inizi=seq(min(pos),max(pos),lstep)
fini=inizi+lwin-1
mid=inizi+(lstep/2)

pbs_pel=pbs_tsi=pbs_chb=c()

for (i in 1:length(inizi)) {

	ind=which(pos>=inizi[i] & pos<fini[i])
	if (length(ind)>0) {
	subhap1=substring(hap1, ind[1], ind[length(ind)])
	subhap2=substring(hap2, ind[1], ind[length(ind)])
	subhap3=substring(hap3, ind[1], ind[length(ind)])
	subhaplos=list(subhap1,subhap2,subhap3)

	fsts=reynolds(subhaplos)
	# 1 is TSI CHB, 2 is TSI PEL, 3 is CHB PEL
	pbs_pel=c(pbs_pel, dopbs(fsts[2],fsts[3],fsts[1]))
	pbs_tsi=c(pbs_pel, dopbs(fsts[1],fsts[2],fsts[3]))
	pbs_chb=c(pbs_pel, dopbs(fsts[1],fsts[3],fsts[2]))
	} else {
	pbs_pel=c(pbs_pel, NA)
        pbs_tsi=c(pbs_pel, NA)
        pbs_chb=c(pbs_pel, NA)
	}
}

cat("Maximum PBS value:", max(pbs_pel, na.rm=T), "\n")

pdf(file=fout)
plot(x=pos, y=pbs_pel, xlab="Chromosome 11", ylab="PBS", ty="l", lty=1, col="red")
lines(x=pos, y=pbs_chb, col="blue")
lines(x=pos, y=pbs_tsi, col="green")
legend("topright", col=c("red","blue","green"), legend=c("PEL","CHB","TSI"))
dev.off()



