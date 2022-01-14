library('vegan')

bray.full.comp = numeric()

for (bray.it in 1:100) {

# generate two samples that differ by some random amount. goal is to have different levels of compositional differences across iterations
niche.breadth = sample(x = c(5,10,25,50,100,200,400,700,1000),size = 1)
Samp1.niche = sample(x = seq(0,1000),size=1)
Samp2.niche = sample(x = seq(0,1000),size=1)

Samp1 = numeric()
Samp2 = numeric()
for (i in 1:1000) {

  Samp1 = rbind(Samp1,c(max(0,rnorm(n = 1,mean = 100*min(i/Samp1.niche,Samp1.niche/i),sd = niche.breadth))))  
  Samp2 = rbind(Samp2,c(max(0,rnorm(n = 1,mean = 100*min(i/Samp2.niche,Samp2.niche/i),sd = niche.breadth))))  

}

# built species by site matrix
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
  
  bray.comp = rbind(bray.comp,c(vegdist(x = t(spxsite.mult),method = "bray"),i))

}

bray.comp = as.data.frame(bray.comp)
colnames(bray.comp) = c("Bray_Curtis","Iteration")

bray.mean.peak.var = mean(bray.comp$Bray_Curtis)

#mod.to.plot = bray.comp$Bray_Curtis ~ bray.comp$Iteration
#plot(mod.to.plot,ylim=c(0,1), ylab="Bray-Curtis",xlab="Iteration",main="Peaks Vary in Detectability Across Peaks but not Across Samples")
#abline(h = vegdist(x = t(spxsite),method = "bray"),col=2,lwd=2)

# based on assumption that a each peak has a different level of detectability across samples and they are uncorrelated across samples
bray.comp = numeric()

for (i in 1:100) {
  
  spxsite.mult = spxsite
  
  samp1.peaks.detect = runif(n = nrow(spxsite),min = 0,max = 100)
  samp2.peaks.detect = runif(n = nrow(spxsite),min = 0,max = 100)
  
  spxsite.mult$Samp1 = spxsite$Samp1*samp1.peaks.detect
  spxsite.mult$Samp2 = spxsite$Samp2*samp2.peaks.detect
  
  bray.comp = rbind(bray.comp,c(vegdist(x = t(spxsite.mult),method = "bray"),i))
  
}

bray.comp = as.data.frame(bray.comp)
colnames(bray.comp) = c("Bray_Curtis","Iteration")

bray.mean.peak.and.samp.var = mean(bray.comp$Bray_Curtis)

#mod.to.plot = bray.comp$Bray_Curtis ~ bray.comp$Iteration
#plot(mod.to.plot,ylim=c(0,1), ylab="Bray-Curtis",xlab="Iteration",main="Peaks Vary in Detectability Across Peaks and Samples")
#abline(h = vegdist(x = t(spxsite),method = "bray"),col=2,lwd=2)

# based on assumption that a each peak has a different level of detectability across samples and they are uncorrelated across samples, but the two samples have different ranges in detectability
bray.comp = numeric()

for (i in 1:100) {
  
  spxsite.mult = spxsite
  
  samp1.peaks.detect = runif(n = nrow(spxsite),min = 0,max = 10)
  samp2.peaks.detect = runif(n = nrow(spxsite),min = 0,max = 100)
  
  spxsite.mult$Samp1 = spxsite$Samp1*samp1.peaks.detect
  spxsite.mult$Samp2 = spxsite$Samp2*samp2.peaks.detect
  
  bray.comp = rbind(bray.comp,c(vegdist(x = t(spxsite.mult),method = "bray"),i))
  
}

bray.comp = as.data.frame(bray.comp)
colnames(bray.comp) = c("Bray_Curtis","Iteration")

bray.mean.samples.diff.detect.range = mean(bray.comp$Bray_Curtis)

#mod.to.plot = bray.comp$Bray_Curtis ~ bray.comp$Iteration
#plot(mod.to.plot,ylim=c(0,1), ylab="Bray-Curtis",xlab="Iteration",main="Samples Differ in Detectability Range")
#abline(h = vegdist(x = t(spxsite),method = "bray"),col=2,lwd=2)

bray.full.comp = rbind(bray.full.comp,c(vegdist(x = t(spxsite),method = "bray"),bray.mean.peak.var,bray.mean.peak.and.samp.var,bray.mean.samples.diff.detect.range))

}

bray.full.comp = as.data.frame(bray.full.comp)
colnames(bray.full.comp) = c("Bray.True","Bray.Peak.Var","Bray.Peak.Samp.Var","Bray.Samp.Detect.Var")

# make plots

pdf("Bray_Simulations.pdf",height = 15)
par(pty="s",cex.lab=2,cex.axis=1.5,mfrow=c(3,1))

mod.to.plot = bray.full.comp$Bray.Peak.Var ~ bray.full.comp$Bray.True
plot(mod.to.plot,ylab="Simulated Bray-Curtis",xlab="True Bray-Curtis",main="Peaks Vary in Detectability",ylim=c(0,1))
abline(0,1,col=2,lwd=2)

mod.to.plot = bray.full.comp$Bray.Peak.Samp.Var ~ bray.full.comp$Bray.True
plot(mod.to.plot,ylab="Simulated Bray-Curtis",xlab="True Bray-Curtis",main="Peak Detectability Varies Across Peaks and Samples",ylim=c(0,1))
abline(0,1,col=2,lwd=2)

mod.to.plot = bray.full.comp$Bray.Samp.Detect.Var ~ bray.full.comp$Bray.True
plot(mod.to.plot,ylab="Simulated Bray-Curtis",xlab="True Bray-Curtis",main="Range in Peak Detectability Varies Across Samples",ylim=c(0,1))
abline(0,1,col=2,lwd=2)

dev.off()
