
# adapted from a script by Javier Mendoza Revilla

args=commandArgs(T)
fin=args[1]
fout=args[2]
rm(args)

source("Scripts/plotGenes.R")

pdf(file=fout)
par(mfrow=c(2,1))
par(mar=c(0, 4, 4, 2) + 0.1)

pbs <- read.table(fin, header=TRUE)

plot(pbs$midPos/1e6, pbs$PBS2, type ="l", col="#d95f02",ylab = "PBS", xlab=paste("Chromosome",pbs$chr[1]), main="PBS scan in PEL" , lwd=1, xaxt="n",cex.main = 0.9, cex.axis = 0.6, cex.lab = 0.68 ) #orange
    #title(ylab="Selection Statistics", line=2.2, cex.lab=0.69)
    #abline(h=1,lty=2,col="grey",lwd=1)


#Add a second panel showing the genes in that region
par(mar=c(5, 4, 0.5, 2) + 0.1)
plotGenes(pbs$chr[1],min(pbs$midPos)/1e6,max(pbst$midPos)/1e6)





