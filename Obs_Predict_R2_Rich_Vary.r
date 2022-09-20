# this is for running more levels of peak richness to build a relationship between R2 of pred-vs-obs and peak richness. 
# expectation is the R2 will go up as peak richness goes up

rm(list=ls())

library('vegan')
fig.made = 0 # toggle setting whether to make the regression of observed and true abundance/differences
rsq.compiled = numeric()

for (peak.rich in c(10,20,40,60,80,100,150,200,250,300,350,400,500,600,700,800,900,1000,1500,2000)) {

bray.full.comp = numeric()

for (bray.it in 1:100) {

# generate two samples that differ by some random amount. goal is to have different levels of compositional differences across iterations
niche.breadth = sample(x = c(1,5,10,25,50,100,200),size = 1)
Samp1.niche = sample(x = seq(0,1000),size=1)
Samp2.niche = sample(x = seq(0,1000),size=1)

Samp1 = numeric()
Samp2 = numeric()
for (i in 1:peak.rich) {

  Samp1 = rbind(Samp1,c(max(0.1,rnorm(n = 1,mean = 100*min(i/Samp1.niche,Samp1.niche/i),sd = niche.breadth))))  
  Samp2 = rbind(Samp2,c(max(0.1,rnorm(n = 1,mean = 100*min(i/Samp2.niche,Samp2.niche/i),sd = niche.breadth))))  

}

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
  
  min.detect = sample(x = c(0),size = 1) # the minimum 'peaks.detect' value (think of it as ionization efficiency) below which we assume the peak can't be observed regardless of its concentration
  
  spxsite.mult = spxsite
  
  peaks.detect = runif(n = nrow(spxsite),min = 0,max = 100)
  peaks.detect[peaks.detect < min.detect] = 0
  
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

########################################################

######## change the error assumption
# based on assumption that a each peak has a different level of detectability across samples and they are uncorrelated across samples
bray.comp = numeric()

for (i in 1:100) {
  
  min.detect = sample(x = c(0),size = 1) # the minimum 'peaks.detect' value (think of it as ionization efficiency) below which we assume the peak can't be observed regardless of its concentration
  
  spxsite.mult = spxsite
  
  samp1.peaks.detect = runif(n = nrow(spxsite),min = 0,max = 100)
  samp2.peaks.detect = runif(n = nrow(spxsite),min = 0,max = 100)
  samp1.peaks.detect[samp1.peaks.detect < min.detect] = 0
  samp2.peaks.detect[samp2.peaks.detect < min.detect] = 0
  
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
  ########
  
bray.comp = as.data.frame(bray.comp)
colnames(bray.comp) = c("Bray.Peak.Samp.Var","Samp1_Shannon.Peak.Samp.Var","Samp2_Shannon.Peak.Samp.Var","Samp1_Weighted.Peak.Samp.Var","Samp2_Weighted.Peak.Samp.Var")

peak.and.samp.var = bray.comp

# make full matrix
bray.full.comp = rbind(bray.full.comp,cbind(rep(vegdist(x = t(spxsite),method = "bray"),nrow(peak.var)),
                                            rep(diversity(spxsite$Samp1, index = "shannon"),nrow(peak.var)),
                                            rep(diversity(spxsite$Samp2, index = "shannon"),nrow(peak.var)),
                                            weighted.mean(as.numeric(row.names(spxsite)),spxsite$Samp1),
                                            weighted.mean(as.numeric(row.names(spxsite)),spxsite$Samp2),
                                            peak.var,peak.and.samp.var))

#check niche breadth and order of magnitude variation
print(c(peak.rich,niche.breadth,round(order.of.mag,digits=2),bray.it,date()))

} # end bray.it loop

bray.full.comp = as.data.frame(bray.full.comp)
colnames(bray.full.comp)[1:5] = c("Bray.True","Samp1.Shannon.True","Samp2.Shannon.True","Samp1.Weighted.True","Samp2.Weighted.True")

# extract R2 of obs vs. true

# bray
mod.to.plot = bray.full.comp$Bray.Peak.Var ~ bray.full.comp$Bray.True
mod.out = summary(lm(mod.to.plot))
Bray.Peak.Var.rsq.obs.vs.true = mod.out$r.squared

mod.to.plot = bray.full.comp$Bray.Peak.Samp.Var ~ bray.full.comp$Bray.True
mod.out = summary(lm(mod.to.plot))
Bray.Peak.Samp.Var.rsq.obs.vs.true = mod.out$r.squared

# shannon
mod.to.plot = bray.full.comp$Samp1_Shannon.Peak.Var ~ bray.full.comp$Samp1.Shannon.True
mod.out = summary(lm(mod.to.plot))
Shannon.Peak.Var.rsq.obs.vs.true = mod.out$r.squared

mod.to.plot = bray.full.comp$Samp1_Shannon.Peak.Samp.Var ~ bray.full.comp$Samp1.Shannon.True
mod.out = summary(lm(mod.to.plot))
Shannon.Peak.Samp.Var.rsq.obs.vs.true = mod.out$r.squared

# weighted trait
mod.to.plot = bray.full.comp$Samp1_Weighted.Peak.Var ~ bray.full.comp$Samp1.Weighted.True
mod.out = summary(lm(mod.to.plot))
Weighted.Peak.Var.rsq.obs.vs.true = mod.out$r.squared

mod.to.plot = bray.full.comp$Samp1_Weighted.Peak.Samp.Var ~ bray.full.comp$Samp1.Weighted.True
mod.out = summary(lm(mod.to.plot))
Weighted.Peak.Samp.Var.rsq.obs.vs.true = mod.out$r.squared

rsq.compiled = rbind(rsq.compiled,c(
  peak.rich,
  Bray.Peak.Var.rsq.obs.vs.true,
  Bray.Peak.Samp.Var.rsq.obs.vs.true,
  Shannon.Peak.Var.rsq.obs.vs.true,
  Shannon.Peak.Samp.Var.rsq.obs.vs.true,
  Weighted.Peak.Var.rsq.obs.vs.true,
  Weighted.Peak.Samp.Var.rsq.obs.vs.true))

} # end richness loop

rsq.compiled = as.data.frame(rsq.compiled)
colnames(rsq.compiled) = c("Peak_Richness","Bray.Peak.Var.rsq.obs.vs.true",
                           "Bray.Peak.Samp.Var.rsq.obs.vs.true",
                           "Shannon.Peak.Var.rsq.obs.vs.true",
                           "Shannon.Peak.Samp.Var.rsq.obs.vs.true",
                           "Weighted.Peak.Var.rsq.obs.vs.true",
                           "Weighted.Peak.Samp.Var.rsq.obs.vs.true")

write.csv(rsq.compiled,"Pred.v.Obs_R2_vs_Rich.csv",quote=F,row.names = F)

#### make figures

pdf("Obs.v.True.Rsq_v_Rich.pdf",height = 15,width = 10)
par(pty="s",mfrow = c(3,2),cex.lab=2,cex.axis=1.5)

# Bray
mod.to.plot = rsq.compiled$Bray.Peak.Var.rsq.obs.vs.true ~ rsq.compiled$Peak_Richness
plot(mod.to.plot,typ="p",ylab="Observed vs. True R2 (Bray-Curtis)",xlab="Peak Richness",ylim=c(0.6,1),main="Among Peak Variation")
points(lowess(mod.to.plot,f = 0.1),typ="l",col=4)

mod.to.plot = rsq.compiled$Bray.Peak.Samp.Var.rsq.obs.vs.true ~ rsq.compiled$Peak_Richness
plot(mod.to.plot,typ="p",ylim=c(0.6,1),xlab="Peak Richness",ylab="Observed vs. True R2 (Bray-Curtis)",main="Among Peak & Sample Variation")
points(lowess(mod.to.plot,f = 0.1),typ="l",col=4)

# Shannon
mod.to.plot = rsq.compiled$Shannon.Peak.Var.rsq.obs.vs.true ~ rsq.compiled$Peak_Richness
plot(mod.to.plot,typ="p",ylim=c(0.6,1),xlab="Peak Richness",ylab="Observed vs. True R2 (Shannon)")
points(lowess(mod.to.plot,f = 0.1),typ="l",col=4)

mod.to.plot = rsq.compiled$Shannon.Peak.Samp.Var.rsq.obs.vs.true ~ rsq.compiled$Peak_Richness
plot(mod.to.plot,typ="p",ylim=c(0.6,1),xlab="Peak Richness",ylab="Observed vs. True R2 (Shannon)")
points(lowess(mod.to.plot,f = 0.1),typ="l",col=4)

# Weighted trait
mod.to.plot = rsq.compiled$Weighted.Peak.Var.rsq.obs.vs.true ~ rsq.compiled$Peak_Richness
plot(mod.to.plot,typ="p",ylim=c(0.6,1),xlab="Peak Richness",ylab="Observed vs. True R2 (Weighted Trait)")
points(lowess(mod.to.plot,f = 0.1),typ="l",col=4)

mod.to.plot = rsq.compiled$Weighted.Peak.Samp.Var.rsq.obs.vs.true ~ rsq.compiled$Peak_Richness
plot(mod.to.plot,typ="p",ylim=c(0.6,1),xlab="Peak Richness",ylab="Observed vs. True R2 (Weighted Trait)")
points(lowess(mod.to.plot,f = 0.1),typ="l",col=4)

dev.off()
