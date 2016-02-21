
library(ape)
library(ade4)
library(adegenet)
library(pegas)

args=commandArgs(T)
fin=args[1]
fout=args[2]
rm(args)

x=read.dna(fin, format="fasta")

h=haplotype(x)
h=subset(h,5)
h <- sort(h, what = "labels")

ind.hap<-with(
	              stack(setNames(attr(h, "index"), rownames(h))),
		              table(hap=ind, pop=rownames(x)[values])
		              )

subset.haplotype <- function(x, freqmin = 1, freqmax = Inf, ...)
{
	    oc <- oldClass(x)
    idx <- attr(x, "index")
        f <- sapply(idx, length)
        s <- f <= freqmax & f >= freqmin
	    x <- x[s, ]
	    attr(x, "index") <- idx[s]
	        class(x) <- oc
	        x
}

net=haploNet(h)

pdf(file=fout)

plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 0.8, fast=TRUE, pie=ind.hap, labels=F, legend=F, show.mutation=0, th=0, lwd = (1 + round(.10*(net[,'step']))))
legend("topleft", colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=20)

dev.off()




