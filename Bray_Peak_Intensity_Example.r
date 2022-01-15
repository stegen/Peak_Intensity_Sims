rm(list=ls())

library('vegan')

for (peak.rich in c(100,1000)) {

bray.full.comp = numeric()

for (bray.it in 1:100) {

# generate two samples that differ by some random amount. goal is to have different levels of compositional differences across iterations
niche.breadth = sample(x = c(1,5,10,25,50,100,200),size = 1)
Samp1.niche = sample(x = seq(0,1000),size=1)
Samp2.niche = sample(x = seq(0,1000),size=1)

Samp1 = numeric()
Samp2 = numeric()
for (i in 1:peak.rich) {

  Samp1 = rbind(Samp1,c(max(0,rnorm(n = 1,mean = 100*min(i/Samp1.niche,Samp1.niche/i),sd = niche.breadth))))  
  Samp2 = rbind(Samp2,c(max(0,rnorm(n = 1,mean = 100*min(i/Samp2.niche,Samp2.niche/i),sd = niche.breadth))))  

}

# cut at 0.01 to get about 3-4 orders of magnitude range in concentrations
#Samp1[Samp1 < 0.01] = 0
#Samp2[Samp2 < 0.01] = 0

# order of mag range in real concentrations
order.of.mag = c(log10(max(Samp1)) - log10(min(Samp1[Samp1>0])),
log10(max(Samp2)) - log10(min(Samp2[Samp2>0])))


# build species by site matrix
spxsite = data.frame("Samp1" = Samp1,"Samp2" = Samp2)
double.zero = which(spxsite$Samp1 == 0 & spxsite$Samp2 == 0)
if (length(double.zero) > 0) {
  spxsite = spxsite[-which(spxsite$Samp1 == 0 & spxsite$Samp2 == 0),]
}
head(spxsite)

# based on assumption that a given peak has a given level of detectability, and that detection ability does not vary across samples
bray.comp = numeric()

for (i in 1:100) {

  spxsite.mult = spxsite
  
  peaks.detect = runif(n = nrow(spxsite),min = 0,max = 100)
  
  spxsite.mult = spxsite*peaks.detect
  
  bray.comp = rbind(bray.comp,c(vegdist(x = t(spxsite.mult),method = "bray"),
                    diversity(spxsite.mult$Samp1, index = "shannon"),
                    diversity(spxsite.mult$Samp2, index = "shannon"),
                    weighted.mean(as.numeric(row.names(spxsite.mult)),spxsite.mult$Samp1),
                    weighted.mean(as.numeric(row.names(spxsite.mult)),spxsite.mult$Samp2)
                    )
  )

}

bray.comp = as.data.frame(bray.comp)
colnames(bray.comp) = c("Bray.Peak.Var","Samp1_Shannon.Peak.Var","Samp2_Shannon.Peak.Var","Samp1_Weighted.Peak.Var","Samp2_Weighted.Peak.Var")

peak.var = bray.comp

# based on assumption that a each peak has a different level of detectability across samples and they are uncorrelated across samples
bray.comp = numeric()

for (i in 1:100) {
  
  spxsite.mult = spxsite
  
  samp1.peaks.detect = runif(n = nrow(spxsite),min = 0,max = 100)
  samp2.peaks.detect = runif(n = nrow(spxsite),min = 0,max = 100)
  
  spxsite.mult$Samp1 = spxsite$Samp1*samp1.peaks.detect
  spxsite.mult$Samp2 = spxsite$Samp2*samp2.peaks.detect
  
  bray.comp = rbind(bray.comp,c(vegdist(x = t(spxsite.mult),method = "bray"),
                                diversity(spxsite.mult$Samp1, index = "shannon"),
                                diversity(spxsite.mult$Samp2, index = "shannon"),
                                weighted.mean(as.numeric(row.names(spxsite.mult)),spxsite.mult$Samp1),
                                weighted.mean(as.numeric(row.names(spxsite.mult)),spxsite.mult$Samp2)
  )
  )
  
}

bray.comp = as.data.frame(bray.comp)
colnames(bray.comp) = c("Bray.Peak.Samp.Var","Samp1_Shannon.Peak.Samp.Var","Samp2_Shannon.Peak.Samp.Var","Samp1_Weighted.Peak.Samp.Var","Samp2_Weighted.Peak.Samp.Var")

peak.and.samp.var = bray.comp

# based on assumption that a each peak has a different level of detectability across samples and they are uncorrelated across samples, but the two samples have different ranges in detectability
bray.comp = numeric()

for (i in 1:100) {
  
  spxsite.mult = spxsite
  
  samp1.peaks.detect = runif(n = nrow(spxsite),min = 0,max = 10)
  samp2.peaks.detect = runif(n = nrow(spxsite),min = 0,max = 100)
  
  spxsite.mult$Samp1 = spxsite$Samp1*samp1.peaks.detect
  spxsite.mult$Samp2 = spxsite$Samp2*samp2.peaks.detect
  
  bray.comp = rbind(bray.comp,c(vegdist(x = t(spxsite.mult),method = "bray"),
                                diversity(spxsite.mult$Samp1, index = "shannon"),
                                diversity(spxsite.mult$Samp2, index = "shannon"),
                                weighted.mean(as.numeric(row.names(spxsite.mult)),spxsite.mult$Samp1),
                                weighted.mean(as.numeric(row.names(spxsite.mult)),spxsite.mult$Samp2)
  )
  )
  
}

bray.comp = as.data.frame(bray.comp)
colnames(bray.comp) = c("Bray.Samp.Detect.Var","Samp1_Shannon.Samp.Detect.Var","Samp2_Shannon.Samp.Detect.Var","Samp1_Weighted.Samp.Detect.Var","Samp2_Weighted.Samp.Detect.Var")

samples.diff.detect.range = bray.comp

# based on assumption that the level of detectability varies with the trait value, but does not vary across samples
bray.comp = numeric()

for (i in 1:100) {
  
  spxsite.mult = spxsite
  
  spxsite.mult = spxsite*as.numeric(row.names(spxsite))
  
  bray.comp = rbind(bray.comp,c(vegdist(x = t(spxsite.mult),method = "bray"),
                                diversity(spxsite.mult$Samp1, index = "shannon"),
                                diversity(spxsite.mult$Samp2, index = "shannon"),
                                weighted.mean(as.numeric(row.names(spxsite.mult)),spxsite.mult$Samp1),
                                weighted.mean(as.numeric(row.names(spxsite.mult)),spxsite.mult$Samp2)
  )
  )
  
}

bray.comp = as.data.frame(bray.comp)
colnames(bray.comp) = c("Bray.Trait.Corr","Samp1_Shannon.Trait.Corr","Samp2_Shannon.Trait.Corr","Samp1_Weighted.Trait.Corr","Samp2_Weighted.Trait.Corr")

trait.corr = bray.comp


# make full matrix
bray.full.comp = rbind(bray.full.comp,cbind(rep(vegdist(x = t(spxsite),method = "bray"),nrow(peak.var)),
                                            rep(diversity(spxsite$Samp1, index = "shannon"),nrow(peak.var)),
                                            rep(diversity(spxsite$Samp2, index = "shannon"),nrow(peak.var)),
                                            weighted.mean(as.numeric(row.names(spxsite)),spxsite$Samp1),
                                            weighted.mean(as.numeric(row.names(spxsite)),spxsite$Samp2),
                                            peak.var,peak.and.samp.var,samples.diff.detect.range,trait.corr))

#check niche breadth and order of magnitude variation
print(c(niche.breadth,order.of.mag))

}

bray.full.comp = as.data.frame(bray.full.comp)
colnames(bray.full.comp)[1:5] = c("Bray.True","Samp1.Shannon.True","Samp2.Shannon.True","Samp1.Weighted.True","Samp2.Weighted.True")

# make plots

pdf(paste0("Bray_Simulations_",peak.rich,".pdf"),height = 20)
par(pty="s",cex.lab=2,cex.axis=1.5,mfrow=c(4,1))

mod.to.plot = bray.full.comp$Bray.Peak.Var ~ bray.full.comp$Bray.True
plot(mod.to.plot,ylab="Simulated Bray-Curtis",xlab="True Bray-Curtis",main="Peaks Vary in Detectability",ylim=c(0,1),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Bray.Peak.Samp.Var ~ bray.full.comp$Bray.True
plot(mod.to.plot,ylab="Simulated Bray-Curtis",xlab="True Bray-Curtis",main="Peak Detectability Varies Across Peaks and Samples",ylim=c(0,1),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Bray.Samp.Detect.Var ~ bray.full.comp$Bray.True
plot(mod.to.plot,ylab="Simulated Bray-Curtis",xlab="True Bray-Curtis",main="Range in Peak Detectability Varies Across Samples",ylim=c(0,1),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Bray.Trait.Corr ~ bray.full.comp$Bray.True
plot(mod.to.plot,ylab="Simulated Bray-Curtis",xlab="True Bray-Curtis",main="Peak Detectability Correlates with Trait Value",ylim=c(0,1),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

dev.off()

# make plots

pdf(paste0("Shannon_Simulations_",peak.rich,".pdf"),height = 20,width=10)
par(pty="s",cex.lab=2,cex.axis=1.5,mfrow=c(4,2))

mod.to.plot = bray.full.comp$Samp1_Shannon.Peak.Var ~ bray.full.comp$Samp1.Shannon.True
plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main="Peaks Vary in Detectability (Samp1)",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Samp2_Shannon.Peak.Var ~ bray.full.comp$Samp2.Shannon.True
plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main="Peaks Vary in Detectability (Samp2)",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Samp1_Shannon.Peak.Samp.Var ~ bray.full.comp$Samp1.Shannon.True
plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main="Peak Detectability Varies Across Peaks and Samples (Samp1)",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Samp2_Shannon.Peak.Samp.Var ~ bray.full.comp$Samp2.Shannon.True
plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main="Peak Detectability Varies Across Peaks and Samples (Samp2)",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Samp1_Shannon.Samp.Detect.Var ~ bray.full.comp$Samp1.Shannon.True
plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main="Range in Peak Detectability Varies Across Samples (Samp1)",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Samp2_Shannon.Samp.Detect.Var ~ bray.full.comp$Samp2.Shannon.True
plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main="Range in Peak Detectability Varies Across Samples (Samp2)",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Samp1_Shannon.Trait.Corr ~ bray.full.comp$Samp1.Shannon.True
plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main="Peak Detectability Correlates with Trait Value (Samp1)",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Samp2_Shannon.Trait.Corr ~ bray.full.comp$Samp2.Shannon.True
plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main="Peak Detectability Correlates with Trait Value (Samp2)",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

dev.off()

#

pdf(paste0("Weighted_Simulations_",peak.rich,".pdf"),height = 15,width=10)
par(pty="s",cex.lab=2,cex.axis=1.5,mfrow=c(3,2))

mod.to.plot = bray.full.comp$Samp1_Weighted.Peak.Var ~ bray.full.comp$Samp1.Weighted.True
plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main="Peaks Vary in Detectability (Samp1)",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Samp2_Weighted.Peak.Var ~ bray.full.comp$Samp2.Weighted.True
plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main="Peaks Vary in Detectability (Samp2)",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Samp1_Weighted.Peak.Samp.Var ~ bray.full.comp$Samp1.Weighted.True
plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main="Peak Detectability Varies Across Peaks and Samples (Samp1)",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Samp2_Weighted.Peak.Samp.Var ~ bray.full.comp$Samp2.Weighted.True
plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main="Peak Detectability Varies Across Peaks and Samples (Samp2)",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Samp1_Weighted.Samp.Detect.Var ~ bray.full.comp$Samp1.Weighted.True
plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main="Range in Peak Detectability Varies Across Samples (Samp1)",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Samp1_Weighted.Samp.Detect.Var ~ bray.full.comp$Samp1.Weighted.True
plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main="Range in Peak Detectability Varies Across Samples (Samp2)",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Samp1_Weighted.Trait.Corr ~ bray.full.comp$Samp1.Weighted.True
plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main="Peak Detectability Correlates with Trait Value (Samp1)",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Samp1_Weighted.Trait.Corr ~ bray.full.comp$Samp1.Weighted.True
plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main="Peak Detectability Correlates with Trait Value (Samp2)",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

dev.off()

}
