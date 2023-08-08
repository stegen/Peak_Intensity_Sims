rm(list=ls())

library('vegan')

complete.metrics = numeric()
unique.id = 1

for (peak.rich in c(100,500,1000,5000,10000)) {

  bray.full.comp = numeric()
    
for (min.detect in c(10,20,80)) {

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

# based on assumption that a each peak has a different level of detectability across samples and they are uncorrelated across samples
# systematically varying the fraction of peaks that have zero ionization (or not enough to detect)
bray.comp = numeric()

for (i in 1:100) {
  
  spxsite.mult = spxsite
  
  samp1.peaks.detect = runif(n = nrow(spxsite),min = 0,max = 100)
  samp2.peaks.detect = runif(n = nrow(spxsite),min = 0,max = 100)

  samp1.peaks.detect[samp1.peaks.detect < min.detect] = 0
  samp2.peaks.detect[samp2.peaks.detect < min.detect] = 0  
  spxsite.mult$Samp1 = spxsite$Samp1*samp1.peaks.detect
  spxsite.mult$Samp2 = spxsite$Samp2*samp2.peaks.detect

  bray.comp = rbind(bray.comp,c(unique.id,
                                vegdist(x = t(spxsite.mult),method = "bray"),
                                diversity(spxsite.mult$Samp1, index = "shannon"),
                                diversity(spxsite.mult$Samp2, index = "shannon"),
                                weighted.mean(as.numeric(row.names(spxsite.mult)),spxsite.mult$Samp1),
                                weighted.mean(as.numeric(row.names(spxsite.mult)),spxsite.mult$Samp2),
                                length(spxsite.mult$Samp1[spxsite.mult$Samp1 > 0]),
                                length(spxsite.mult$Samp2[spxsite.mult$Samp2 > 0])
  )
  )
  
}

bray.comp = as.data.frame(bray.comp)
colnames(bray.comp) = c("Unique.ID","Bray.Peak.Samp.Var","Samp1_Shannon.Peak.Samp.Var","Samp2_Shannon.Peak.Samp.Var","Samp1_Weighted.Peak.Samp.Var","Samp2_Weighted.Peak.Samp.Var","Samp1_Richness","Samp2_Richness")

peak.and.samp.var = bray.comp

unique.id = unique.id + 1

#####

# make full matrix
bray.full.comp = rbind(bray.full.comp,cbind(rep(min.detect,nrow(peak.and.samp.var)),rep(vegdist(x = t(spxsite),method = "bray"),nrow(peak.and.samp.var)),
                                            rep(diversity(spxsite$Samp1, index = "shannon"),nrow(peak.and.samp.var)),
                                            rep(diversity(spxsite$Samp2, index = "shannon"),nrow(peak.and.samp.var)),
                                            weighted.mean(as.numeric(row.names(spxsite)),spxsite$Samp1),
                                            weighted.mean(as.numeric(row.names(spxsite)),spxsite$Samp2),
                                            length(spxsite$Samp1[spxsite$Samp1 > 0]),
                                            length(spxsite$Samp2[spxsite$Samp2 > 0]),
                                            peak.and.samp.var))

#check niche breadth and order of magnitude variation
print(c(niche.breadth,order.of.mag,bray.it,date(),min.detect,peak.rich,length(spxsite.mult$Samp1[spxsite.mult$Samp1 > 0])))

}

}

bray.full.comp = as.data.frame(bray.full.comp)
colnames(bray.full.comp)[1:8] = c("Min.Detect","Bray.True","Samp1.Shannon.True","Samp2.Shannon.True","Samp1.Weighted.True","Samp2.Weighted.True","Samp1.Richness.True","Samp2.Richness.True")

complete.metrics = rbind(complete.metrics,bray.full.comp)

# make plots

pdf(paste0("Bray_Simulations_",peak.rich,"_Detect_Varies.pdf"),height = 15)
par(pty="s",cex.lab=2,cex.axis=1.5,mfrow=c(3,1))

unique.min.detect = 1
mod.to.plot = bray.full.comp$Bray.Peak.Samp.Var[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])] ~ bray.full.comp$Bray.True[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])]
plot(mod.to.plot,ylab="Simulated Bray-Curtis",xlab="True Bray-Curtis",main=paste0("Peak Detect Varies Across Peaks & Samples (Min Detect = ",unique(bray.full.comp$Min.Detect)[unique.min.detect],")"),ylim=c(0,1),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

unique.min.detect = 2
mod.to.plot = bray.full.comp$Bray.Peak.Samp.Var[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])] ~ bray.full.comp$Bray.True[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])]
plot(mod.to.plot,ylab="Simulated Bray-Curtis",xlab="True Bray-Curtis",main=paste0("Peak Detect Varies Across Peaks & Samples (Min Detect = ",unique(bray.full.comp$Min.Detect)[unique.min.detect],")"),ylim=c(0,1),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

unique.min.detect = 3
mod.to.plot = bray.full.comp$Bray.Peak.Samp.Var[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])] ~ bray.full.comp$Bray.True[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])]
plot(mod.to.plot,ylab="Simulated Bray-Curtis",xlab="True Bray-Curtis",main=paste0("Peak Detect Varies Across Peaks & Samples (Min Detect = ",unique(bray.full.comp$Min.Detect)[unique.min.detect],")"),ylim=c(0,1),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

dev.off()

#### need to edit below here
# make plots

pdf(paste0("Shannon_Simulations_",peak.rich,"_Detect_Varies.pdf"),height = 15,width=10)
par(pty="s",cex.lab=2,cex.axis=1.5,mfrow=c(3,2))

unique.min.detect = 1
mod.to.plot = bray.full.comp$Samp1_Shannon.Peak.Samp.Var[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])] ~ bray.full.comp$Samp1.Shannon.True[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])]
plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main=paste0("Peak Detect Varies Across Peaks & Samples (Samp1, Min Detect = ",unique(bray.full.comp$Min.Detect)[unique.min.detect],")"),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

unique.min.detect = 1
mod.to.plot = bray.full.comp$Samp2_Shannon.Peak.Samp.Var[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])] ~ bray.full.comp$Samp2.Shannon.True[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])]
plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main=paste0("Peak Detect Varies Across Peaks & Samples (Samp2, Min Detect = ",unique(bray.full.comp$Min.Detect)[unique.min.detect],")"),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

unique.min.detect = 2
mod.to.plot = bray.full.comp$Samp1_Shannon.Peak.Samp.Var[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])] ~ bray.full.comp$Samp1.Shannon.True[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])]
plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main=paste0("Peak Detect Varies Across Peaks & Samples (Samp1, Min Detect = ",unique(bray.full.comp$Min.Detect)[unique.min.detect],")"),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

unique.min.detect = 2
mod.to.plot = bray.full.comp$Samp2_Shannon.Peak.Samp.Var[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])] ~ bray.full.comp$Samp2.Shannon.True[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])]
plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main=paste0("Peak Detect Varies Across Peaks & Samples (Samp2, Min Detect = ",unique(bray.full.comp$Min.Detect)[unique.min.detect],")"),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

unique.min.detect = 3
mod.to.plot = bray.full.comp$Samp1_Shannon.Peak.Samp.Var[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])] ~ bray.full.comp$Samp1.Shannon.True[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])]
plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main=paste0("Peak Detect Varies Across Peaks & Samples (Samp1, Min Detect = ",unique(bray.full.comp$Min.Detect)[unique.min.detect],")"),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

unique.min.detect = 3
mod.to.plot = bray.full.comp$Samp2_Shannon.Peak.Samp.Var[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])] ~ bray.full.comp$Samp2.Shannon.True[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])]
plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main=paste0("Peak Detect Varies Across Peaks & Samples (Samp2, Min Detect = ",unique(bray.full.comp$Min.Detect)[unique.min.detect],")"),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

dev.off()

#

pdf(paste0("Weighted_Simulations_",peak.rich,"_Detect_Varies.pdf"),height = 15,width=10)
par(pty="s",cex.lab=2,cex.axis=1.5,mfrow=c(3,2))

unique.min.detect = 1
mod.to.plot = bray.full.comp$Samp1_Weighted.Peak.Samp.Var[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])] ~ bray.full.comp$Samp1.Weighted.True[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])]
plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main=paste0("Peak Detect Varies Across Peaks & Samples (Samp1, Min Detect = ",unique(bray.full.comp$Min.Detect)[unique.min.detect],")"),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

unique.min.detect = 1
mod.to.plot = bray.full.comp$Samp2_Weighted.Peak.Samp.Var[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])] ~ bray.full.comp$Samp2.Weighted.True[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])]
plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main=paste0("Peak Detect Varies Across Peaks & Samples (Samp1, Min Detect = ",unique(bray.full.comp$Min.Detect)[unique.min.detect],")"),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

unique.min.detect = 2
mod.to.plot = bray.full.comp$Samp1_Weighted.Peak.Samp.Var[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])] ~ bray.full.comp$Samp1.Weighted.True[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])]
plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main=paste0("Peak Detect Varies Across Peaks & Samples (Samp1, Min Detect = ",unique(bray.full.comp$Min.Detect)[unique.min.detect],")"),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

unique.min.detect = 2
mod.to.plot = bray.full.comp$Samp2_Weighted.Peak.Samp.Var[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])] ~ bray.full.comp$Samp2.Weighted.True[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])]
plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main=paste0("Peak Detect Varies Across Peaks & Samples (Samp1, Min Detect = ",unique(bray.full.comp$Min.Detect)[unique.min.detect],")"),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

unique.min.detect = 3
mod.to.plot = bray.full.comp$Samp1_Weighted.Peak.Samp.Var[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])] ~ bray.full.comp$Samp1.Weighted.True[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])]
plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main=paste0("Peak Detect Varies Across Peaks & Samples (Samp1, Min Detect = ",unique(bray.full.comp$Min.Detect)[unique.min.detect],")"),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

unique.min.detect = 3
mod.to.plot = bray.full.comp$Samp2_Weighted.Peak.Samp.Var[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])] ~ bray.full.comp$Samp2.Weighted.True[which(bray.full.comp$Min.Detect == unique(bray.full.comp$Min.Detect)[unique.min.detect])]
plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main=paste0("Peak Detect Varies Across Peaks & Samples (Samp1, Min Detect = ",unique(bray.full.comp$Min.Detect)[unique.min.detect],")"),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

dev.off()

}

# analyses with complete.metric data frame

complete.metrics$Shannon.diff = 100*abs(complete.metrics$Samp1_Shannon.Peak.Samp.Var - complete.metrics$Samp1.Shannon.True)/complete.metrics$Samp1.Shannon.True
complete.metrics$Bray.diff = 100*abs(complete.metrics$Bray.Peak.Samp.Var - complete.metrics$Bray.True)/complete.metrics$Bray.True
complete.metrics$Weighted.diff = 100*abs(complete.metrics$Samp1_Weighted.Peak.Samp.Var - complete.metrics$Samp1.Weighted.True)/complete.metrics$Samp1.Weighted.True
complete.metrics$Mean.obs.rich = rowMeans(complete.metrics[,c('Samp1_Richness','Samp2_Richness')])

write.csv(complete.metrics,"Raw_Metrics_Min.Detect_Variation.csv",quote = F,row.names = F)


# plot difference between true and observed Shannon
# only doing for min.detect = 10

# shannon error as function of observed richness
mod.to.plot = complete.metrics$Shannon.diff[which(complete.metrics$Min.Detect == 10)] ~ log10(complete.metrics$Samp1_Richness[which(complete.metrics$Min.Detect == 10)])
plot(mod.to.plot,xlab="Log10(Observed Peak Richness)",ylab="Shannon (% Error)",cex=0.4)

# bray error as function of observed mean richness
mod.to.plot = complete.metrics$Bray.diff[which(complete.metrics$Min.Detect == 10)] ~ log10(complete.metrics$Samp1_Richness[which(complete.metrics$Min.Detect == 10)])
plot(mod.to.plot,xlab="Log10(Observed Peak Richness)",ylab="Bray-Curtis (% Error)",cex=0.4)

# weighted trait error as function of observed richness
mod.to.plot = complete.metrics$Weighted.diff[which(complete.metrics$Min.Detect == 10)] ~ log10(complete.metrics$Samp1_Richness[which(complete.metrics$Min.Detect == 10)])
plot(mod.to.plot,xlab="Log10(Observed Peak Richness)",ylab="Weighted Trait (% Error)",cex=0.4)

# calculate coef of variance around each true value and plot as function of observed peak richness
# only doing for min.detect = 10
# only doing Samp1 for the within-sample metrics

ids.to.use = unique(complete.metrics$Unique.ID[which(complete.metrics$Min.Detect == 10)])
var.comp = numeric()

for (i in ids.to.use) {
  
  shannon.cv = sd(x = complete.metrics$Samp1_Shannon.Peak.Samp.Var[which(complete.metrics$Unique.ID == i)])/mean(x = complete.metrics$Samp1_Shannon.Peak.Samp.Var[which(complete.metrics$Unique.ID == i)])
  bray.cv = sd(x = complete.metrics$Bray.Peak.Samp.Var[which(complete.metrics$Unique.ID == i)])/mean(x = complete.metrics$Bray.Peak.Samp.Var[which(complete.metrics$Unique.ID == i)])
  weighted.cv = sd(x = complete.metrics$Samp1_Weighted.Peak.Samp.Var[which(complete.metrics$Unique.ID == i)])/mean(x = complete.metrics$Samp1_Weighted.Peak.Samp.Var[which(complete.metrics$Unique.ID == i)])
  samp1.mean.rich = mean(complete.metrics$Samp1_Richness[which(complete.metrics$Unique.ID == i)])
  two.samp.mean.rich = mean(c(complete.metrics$Samp1_Richness[which(complete.metrics$Unique.ID == i)],complete.metrics$Samp2_Richness[which(complete.metrics$Unique.ID == i)]))
  
  var.comp = rbind(var.comp,c(shannon.cv,
                              bray.cv,
                              weighted.cv,
                              samp1.mean.rich,
                              two.samp.mean.rich))
}

var.comp = as.data.frame(var.comp)
colnames(var.comp) = c("Shannon.CV","Bray.CV","Weighted.CV","Samp1.Mean.Rich","Two.Samp.Mean.Rich")

write.csv(var.comp,"CV_Richness_Compilation.csv",quote = F,row.names = F)

# plot coef of variance as function of richness

pdf("CV_vs_Richness.pdf",height=15)

  par(pty="s",cex.lab=2,cex.axis=1.5,mfrow=c(3,1))

  # shannon variance
  mod.to.plot = 100*var.comp$Shannon.CV ~ log10(var.comp$Samp1.Mean.Rich)
  plot(mod.to.plot,xlab="Log10(Mean Peak Richness)",ylab="Shannon CV (%)",cex=0.4)
  points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

  # bray variance
  mod.to.plot = 100*var.comp$Bray.CV ~ log10(var.comp$Two.Samp.Mean.Rich)
  plot(mod.to.plot,xlab="Log10(Mean Peak Richness)",ylab="Bray-Curtis CV (%)",cex=0.4)
  points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

  # weighted variance
  mod.to.plot = 100*var.comp$Weighted.CV ~ log10(var.comp$Samp1.Mean.Rich)
  plot(mod.to.plot,xlab="Log10(Mean Peak Richness)",ylab="Weighted Trait CV (%)",cex=0.4)
  points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

dev.off()
