
library(pegas)

args=commandArgs(T)
fin=args[1]
fout=args[2]
rm(args)

x=read.dna("fin", format="fasta")

h=haplotype(x)
h <- sort(h, what = "labels")

ind.hap<-with(
	              stack(setNames(attr(h, "index"), rownames(h))),
		              table(hap=ind, pop=rownames(x)[values])
		              )

net=haploNet(h)

pdf(file=fout)

plot(net, size=5*log2(attr(net, "freq"))+5, scale.ratio = 50, cex = 0.8, pie=ind.hap, labels=F, legend=F, show.mutation=0, th=0, lwd = (1 + round(.15*(net[,'step']))))
legend("topleft", colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=20)

dev.off()




