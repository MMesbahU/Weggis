
# adapted from a script by Javier Mendoza Revilla

args=commandArgs(T)
fin=args[1]
fout=args[2]
rm(args)

source("Scripts/plotGenes.R")

pdf(file=fout)
par(mfrow=c(2,1))
par(mar=c(0, 4, 4, 2) + 0.1)

nsl <- read.table(fin, header=TRUE)

cat("Maximum XP_EHH value:", max(nsl[,8], na.rm=T), "\n")

plot(nsl[,2]/1e6, nsl[,8], type ="l", col="orange",ylab = "nSL", xlab="Chromosome", main="XP-EHH scan" , lwd=1, xaxt="n",cex.main = 0.9, cex.axis = 0.6, cex.lab = 0.68 ) #orange

#Add a second panel showing the genes in that region
par(mar=c(5, 4, 0.5, 2) + 0.1)
plotGenes(11,min(nsl[,2]/1e6),max(nsl[,2]/1e6))





