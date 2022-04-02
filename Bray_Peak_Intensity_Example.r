rm(list=ls())

library('vegan')
fig.made = 0 # toggle setting whether to make the regression of observed and true abundance/differences
rsq.compiled = numeric()

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

  Samp1 = rbind(Samp1,c(max(0.1,rnorm(n = 1,mean = 100*min(i/Samp1.niche,Samp1.niche/i),sd = niche.breadth))))  
  Samp2 = rbind(Samp2,c(max(0.1,rnorm(n = 1,mean = 100*min(i/Samp2.niche,Samp2.niche/i),sd = niche.breadth))))  

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
  
  # maybe build another setup that systematically varies the prob of zero)
  #min.detect = sample(x = c(1:10),size = 1) # the minimum 'peaks.detect' value (think of it as ionization efficiency) below which we assume the peak can't be observed regardless of its concentration
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
  
  if (i %in% seq(10,100,by = 10)) { # only doing the below for 10 iterations
  
  ########
  #find r.sq of within and between peak relationships
  
  # first remove peaks that had true zero abundance in either sample
  #spxsite.error.1 = spxsite.mult[-which(spxsite$Samp1 == 0 | spxsite$Samp2 == 0),]
  #spxsite.true = spxsite[-which(spxsite$Samp1 == 0 | spxsite$Samp2 == 0),]
  spxsite.error.1 = spxsite.mult
  spxsite.true = spxsite
  if (identical(rownames(spxsite.error.1),rownames(spxsite.true)) != T) {
    
    print("Error: subsets for spxsite do not match")
    break()
    
  }
  
  # find rsq for observed vs. true abundance
  mod.to.plot = I(spxsite.error.1$Samp1/max(spxsite.error.1$Samp1)) ~ I(spxsite.true$Samp1/max(spxsite.true$Samp1))
  mod = summary(lm(mod.to.plot))
  obs.vs.true.abund.rsq = mod$r.squared
  rm('mod','mod.to.plot')
  
  # find rsq for between-peak differences in observed abundance vs. between-peak differences in true abundance (using only 1 sample)
  temp.y = numeric()
  for (first.peak in 1:I(nrow(spxsite.error.1)-1)) {
    
    temp.y = c(temp.y, I( spxsite.error.1$Samp1[first.peak] - spxsite.error.1$Samp1[I(first.peak+1):nrow(spxsite.error.1)]))
    #print(first.peak)
  }
  
  temp.x = numeric()
  for (first.peak in 1:I(nrow(spxsite.true)-1)) {
    
    temp.x = c(temp.x, I( spxsite.true$Samp1[first.peak] - spxsite.true$Samp1[I(first.peak+1):nrow(spxsite.true)]))
    
  }
  
  if (length(temp.x)-I(nrow(spxsite.true)*(nrow(spxsite.true)-1)/2) != 0) {
    
    print("Error in between peak comparisons")
    break()
    
  }
  
  mod.to.plot = temp.y ~ temp.x 
  full.mod = summary(lm(mod.to.plot))
  btw.peak.obs.diff.vs.true.diff.rsq = full.mod$r.squared
  rm('mod.to.plot','full.mod','temp.y','temp.x')
  
  # find rsq for within-peak differences in observed abundance vs. within-peak differences in true abundance (between samples)
  true.within.diff = spxsite.true$Samp1 - spxsite.true$Samp2
  error.within.diff = spxsite.error.1$Samp1 - spxsite.error.1$Samp2
  temp.y = 2*(true.within.diff - min(true.within.diff))/(max(true.within.diff)-min(true.within.diff)) - 1
  temp.x = 2*(error.within.diff - min(error.within.diff))/(max(error.within.diff)-min(error.within.diff)) - 1
  #temp.y = true.within.diff/max(true.within.diff)
  #temp.x = error.within.diff/max(error.within.diff)
  mod.to.plot = temp.y ~ temp.x 
  mod = summary(lm(mod.to.plot))
  wth.peak.obs.diff.vs.true.diff.rsq = mod$r.squared
  rm('temp.y','temp.x','mod','true.within.diff','error.within.diff')
  
  rsq.compiled = rbind(rsq.compiled,c(
    "Detectability.Same.Between.Samples",
    peak.rich,
    nrow(spxsite),
    nrow(spxsite.true),
    niche.breadth,
    bray.it,
    min.detect,
    Samp1.niche,
    Samp2.niche,
    wth.peak.obs.diff.vs.true.diff.rsq,
    btw.peak.obs.diff.vs.true.diff.rsq,
    obs.vs.true.abund.rsq
  ))
  
  }
  ########
}

bray.comp = as.data.frame(bray.comp)
colnames(bray.comp) = c("Bray.Peak.Var","Samp1_Shannon.Peak.Var","Samp2_Shannon.Peak.Var","Samp1_Weighted.Peak.Var","Samp2_Weighted.Peak.Var")

peak.var = bray.comp

##################################################
if (obs.vs.true.abund.rsq > 0.5 & obs.vs.true.abund.rsq < 0.6 &
    btw.peak.obs.diff.vs.true.diff.rsq > 0.47 & btw.peak.obs.diff.vs.true.diff.rsq < 0.57 &
    wth.peak.obs.diff.vs.true.diff.rsq > 0.69 & wth.peak.obs.diff.vs.true.diff.rsq < 0.79 &
    fig.made == 0 & peak.rich == 1000
    ) {
  
  fig.made = 1 # changes the toggle to make the second panels
  
  # first remove peaks that had true zero abundance in either sample
  #spxsite.error.1 = spxsite.mult[-which(spxsite$Samp1 == 0 | spxsite$Samp2 == 0),]
  #spxsite.true = spxsite[-which(spxsite$Samp1 == 0 | spxsite$Samp2 == 0),]
  spxsite.error.1 = spxsite.mult
  spxsite.true = spxsite
  
  if (identical(rownames(spxsite.error.1),rownames(spxsite.true)) != T) {
    
    print("Error: subsets for spxsite do not match")
    break()
    
  }
  
  pdf("Obs.Abund_vs_True.Abund.pdf",height=5,width=5)
  # plot observed vs. true
  mod.to.plot = I(spxsite.error.1$Samp1/max(spxsite.error.1$Samp1)) ~ I(spxsite.true$Samp1/max(spxsite.true$Samp1))
  par(pty="s")
  plot(mod.to.plot,xlab="True Abundance",ylab="Observed Intensity",cex.lab=2,cex.axis=1.5,cex=0.3,main="Wth-peak; Wth-sample")
  #abline(0,1,col=2,lwd=2)
  mod = summary(lm(mod.to.plot))
  r.sq = round(mod$r.squared,digits = 2)
  abline(mod,lwd=2,col=4)
  mtext(text = substitute(paste(" ",R^2," = ", r.sq), list(r.sq=r.sq)),line = -2,adj = 0,side = 3,cex=1.5)
  
  #mod.to.plot = log10(I(spxsite.error.1$Samp1/max(spxsite.error.1$Samp1))) ~ log10(I(spxsite.true$Samp1/max(spxsite.true$Samp1)))
  #plot(mod.to.plot,xlab="Log10(True Abundance)",ylab="Log10(Observed Intensity)",cex.lab=2,cex.axis=1.5,cex=0.3)
  #abline(0,1,col=2,lwd=2)
  #mod = summary(lm(mod.to.plot))
  #abline(mod,lwd=2,col=4)
  rm('mod','r.sq')
  dev.off()
  
  pdf("Obs_v_True_Diff.pdf",height = 10,width = 10)
### make figures for
# between-peak differences in intensity vs. between-peak differences in true abundance (could use both samples)
# within-peak differences in intensity vs. within-peak differences in true abundance (between samples)

# between-peak differences in intensity vs. between-peak differences in true abundance (using only 1 sample)
temp.y = numeric()
for (first.peak in 1:I(nrow(spxsite.error.1)-1)) {
  
  temp.y = c(temp.y, I( spxsite.error.1$Samp1[first.peak] - spxsite.error.1$Samp1[I(first.peak+1):nrow(spxsite.error.1)]))
    
}

temp.x = numeric()
for (first.peak in 1:I(nrow(spxsite.true)-1)) {
  
  temp.x = c(temp.x, I( spxsite.true$Samp1[first.peak] - spxsite.true$Samp1[I(first.peak+1):nrow(spxsite.true)]))

}

print(c(length(temp.x)-I(nrow(spxsite.true)*(nrow(spxsite.true)-1)/2)," This should be zero"))

#mod.to.plot = temp.y ~ temp.x 
#full.mod = summary(lm(mod.to.plot))

rand.subset = sample(x = 1:I((nrow(spxsite.true)-1)*nrow(spxsite.true)/2),size = 1000,replace = F)
temp.y = 2*(temp.y - min(temp.y))/(max(temp.y)-min(temp.y)) - 1
temp.x = 2*(temp.x - min(temp.x))/(max(temp.x)-min(temp.x)) - 1
temp.y = temp.y[rand.subset]
temp.x = temp.x[rand.subset]
mod.to.plot = temp.y ~ temp.x 
par(pty="s",mfrow=c(2,2))
plot(mod.to.plot,xlab="True Difference",ylab="Observed Difference",cex.lab=2,cex.axis=1.5,cex=0.3,main="Btw-peak diff.; Same error btw samples")
mod = summary(lm(mod.to.plot))
r.sq = round(mod$r.squared,digits = 2)
abline(mod,col=2,lwd=2)
mtext(text = "A ",line = -2,adj = 1,side = 1,cex = 1.5)
mtext(text = substitute(paste(" ",R^2," = ", r.sq), list(r.sq=r.sq)),line = -2,adj = 0,side = 3,cex=1.5)
rm('temp.y','temp.x','mod','rand.subset','r.sq')

# within-peak differences in intensity vs. within-peak differences in true abundance (between samples)
true.within.diff = spxsite.true$Samp1 - spxsite.true$Samp2
error.within.diff = spxsite.error.1$Samp1 - spxsite.error.1$Samp2
temp.x = 2*(true.within.diff - min(true.within.diff))/(max(true.within.diff)-min(true.within.diff)) - 1
temp.y = 2*(error.within.diff - min(error.within.diff))/(max(error.within.diff)-min(error.within.diff)) - 1
#temp.y = true.within.diff/max(true.within.diff)
#temp.x = error.within.diff/max(error.within.diff)
mod.to.plot = temp.y ~ temp.x 
#par(pty="s",cex.lab=2,cex.axis=1.5)
plot(mod.to.plot,xlab="True Difference",ylab="Observed Difference",cex.lab=2,cex.axis=1.5,cex=0.3,main="Wth-peak diff.; Same error btw samples")
mod = summary(lm(mod.to.plot))
r.sq = round(mod$r.squared,digits = 2)
abline(mod,col=2,lwd=2)
mtext(text = "B ",line = -2,adj = 1,side = 1,cex = 1.5)
mtext(text = substitute(paste(" ",R^2," = ", r.sq), list(r.sq=r.sq)),line = -2,adj = 0,side = 3,cex=1.5)
rm('temp.y','temp.x','mod','true.within.diff','error.within.diff','r.sq')

}
########################################################

######## change the error assumption
# based on assumption that a each peak has a different level of detectability across samples and they are uncorrelated across samples
bray.comp = numeric()

for (i in 1:100) {
  
  #min.detect = sample(x = c(1:10),size = 1) # the minimum 'peaks.detect' value (think of it as ionization efficiency) below which we assume the peak can't be observed regardless of its concentration
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
  
  if (i %in% seq(10,100,by = 10)) { # only doing the below for 10 iterations
    
  ########
  #find r.sq of within and between peak relationships
  
  # first remove peaks that had true zero abundance in either sample
  #spxsite.error.2 = spxsite.mult[-which(spxsite$Samp1 == 0 | spxsite$Samp2 == 0),]
  #spxsite.true = spxsite[-which(spxsite$Samp1 == 0 | spxsite$Samp2 == 0),]
  spxsite.error.2 = spxsite.mult
  spxsite.true = spxsite
  
  if (identical(rownames(spxsite.error.2),rownames(spxsite.true)) != T) {
    
    print("Error: subsets for spxsite do not match")
    break()
    
  }
  
  # find rsq for observed vs. true abundance
  mod.to.plot = I(spxsite.error.2$Samp1/max(spxsite.error.2$Samp1)) ~ I(spxsite.true$Samp1/max(spxsite.true$Samp1))
  mod = summary(lm(mod.to.plot))
  obs.vs.true.abund.rsq = mod$r.squared
  rm('mod','mod.to.plot')
  
  # find rsq for between-peak differences in observed abundance vs. between-peak differences in true abundance (using only 1 sample)
  temp.y = numeric()
  for (first.peak in 1:I(nrow(spxsite.error.2)-1)) {
    
    temp.y = c(temp.y, I( spxsite.error.2$Samp1[first.peak] - spxsite.error.2$Samp1[I(first.peak+1):nrow(spxsite.error.2)]))
    
  }
  
  temp.x = numeric()
  for (first.peak in 1:I(nrow(spxsite.true)-1)) {
    
    temp.x = c(temp.x, I( spxsite.true$Samp1[first.peak] - spxsite.true$Samp1[I(first.peak+1):nrow(spxsite.true)]))
    
  }
  
  if (length(temp.x)-I(nrow(spxsite.true)*(nrow(spxsite.true)-1)/2) != 0) {
    
    print("Error in between peak comparisons")
    break()
    
  }
  
  mod.to.plot = temp.y ~ temp.x 
  full.mod = summary(lm(mod.to.plot))
  btw.peak.obs.diff.vs.true.diff.rsq = full.mod$r.squared
  rm('mod.to.plot','full.mod','temp.y','temp.x')
  
  # find rsq for within-peak differences in observed abundance vs. within-peak differences in true abundance (between samples)
  true.within.diff = spxsite.true$Samp1 - spxsite.true$Samp2
  error.within.diff = spxsite.error.2$Samp1 - spxsite.error.2$Samp2
  temp.y = 2*(true.within.diff - min(true.within.diff))/(max(true.within.diff)-min(true.within.diff)) - 1
  temp.x = 2*(error.within.diff - min(error.within.diff))/(max(error.within.diff)-min(error.within.diff)) - 1
  #temp.y = true.within.diff/max(true.within.diff)
  #temp.x = error.within.diff/max(error.within.diff)
  mod.to.plot = temp.y ~ temp.x 
  mod = summary(lm(mod.to.plot))
  wth.peak.obs.diff.vs.true.diff.rsq = mod$r.squared
  rm('temp.y','temp.x','mod','true.within.diff','error.within.diff')
  
  rsq.compiled = rbind(rsq.compiled,c(
    "Detectability.Differs.Between.Samples",
    peak.rich,
    nrow(spxsite),
    nrow(spxsite.true),
    niche.breadth,
    bray.it,
    min.detect,
    Samp1.niche,
    Samp2.niche,
    wth.peak.obs.diff.vs.true.diff.rsq,
    btw.peak.obs.diff.vs.true.diff.rsq,
    obs.vs.true.abund.rsq
  ))
  
  }
  ########
  
}

bray.comp = as.data.frame(bray.comp)
colnames(bray.comp) = c("Bray.Peak.Samp.Var","Samp1_Shannon.Peak.Samp.Var","Samp2_Shannon.Peak.Samp.Var","Samp1_Weighted.Peak.Samp.Var","Samp2_Weighted.Peak.Samp.Var")

peak.and.samp.var = bray.comp

##################################################
if (fig.made == 1 & 
    wth.peak.obs.diff.vs.true.diff.rsq > 0.48 & wth.peak.obs.diff.vs.true.diff.rsq < 0.58) { #toggle set automatically
### make figures for
# between-peak differences in intensity vs. between-peak differences in true abundance (could use both samples)
# within-peak differences in intensity vs. within-peak differences in true abundance (between samples)

fig.made = 2 # this toggles these figures off for the rest of the iterations
# first remove peaks that had true zero abundance in either sample
#spxsite.error.2 = spxsite.mult[-which(spxsite$Samp1 == 0 | spxsite$Samp2 == 0),]
#spxsite.true = spxsite[-which(spxsite$Samp1 == 0 | spxsite$Samp2 == 0),]
spxsite.error.2 = spxsite.mult
spxsite.true = spxsite

if (identical(rownames(spxsite.error.2),rownames(spxsite.true)) != T) {
  
  print("Error: subsets for spxsite do not match")
  break()
  
}

# between-peak differences in intensity vs. between-peak differences in true abundance (using only 1 sample)
temp.y = numeric()
for (first.peak in 1:I(nrow(spxsite.error.2)-1)) {
  
  temp.y = c(temp.y, I( spxsite.error.2$Samp1[first.peak] - spxsite.error.2$Samp1[I(first.peak+1):nrow(spxsite.error.2)]))
  
}

temp.x = numeric()
for (first.peak in 1:I(nrow(spxsite.true)-1)) {
  
  temp.x = c(temp.x, I( spxsite.true$Samp1[first.peak] - spxsite.true$Samp1[I(first.peak+1):nrow(spxsite.true)]))
  
}

if (length(temp.x)-I(nrow(spxsite.true)*(nrow(spxsite.true)-1)/2) != 0) {
  
  print("Error in between peak comparisons")
  break()
  
}

#mod.to.plot = temp.y ~ temp.x 
#full.mod = summary(lm(mod.to.plot))

print(c(length(temp.x)-I(nrow(spxsite.true)*(nrow(spxsite.true)-1)/2)," This should be zero"))

#mod.to.plot = temp.y ~ temp.x 
#full.mod = summary(lm(mod.to.plot))

rand.subset = sample(x = 1:I((nrow(spxsite.true)-1)*nrow(spxsite.true)/2),size = 1000,replace = F)
temp.y = 2*(temp.y - min(temp.y))/(max(temp.y)-min(temp.y)) - 1
temp.x = 2*(temp.x - min(temp.x))/(max(temp.x)-min(temp.x)) - 1
temp.y = temp.y[rand.subset]
temp.x = temp.x[rand.subset]
mod.to.plot = temp.y ~ temp.x 
#par(pty="s",cex.lab=2,cex.axis=1.5,mfrow=c(2,2),cex=0.3)
plot(mod.to.plot,xlab="True Difference",ylab="Observed Difference",cex.lab=2,cex.axis=1.5,cex=0.3,main="Btw-peak diff.; Diff. error btw samples")
mod = summary(lm(mod.to.plot))
r.sq = round(mod$r.squared,digits = 2)
abline(mod,col=2,lwd=2)
mtext(text = "C ",line = -2,adj = 1,side = 1,cex = 1.5)
mtext(text = substitute(paste(" ",R^2," = ", r.sq), list(r.sq=r.sq)),line = -2,adj = 0,side = 3,cex=1.5)
rm('temp.y','temp.x','mod','rand.subset','r.sq')

# within-peak differences in intensity vs. within-peak differences in true abundance (between samples)
true.within.diff = spxsite.true$Samp1 - spxsite.true$Samp2
error.within.diff = spxsite.error.2$Samp1 - spxsite.error.2$Samp2
temp.x = 2*(true.within.diff - min(true.within.diff))/(max(true.within.diff)-min(true.within.diff)) - 1
temp.y = 2*(error.within.diff - min(error.within.diff))/(max(error.within.diff)-min(error.within.diff)) - 1
mod.to.plot = temp.y ~ temp.x 
#par(pty="s",cex.lab=2,cex.axis=1.5)
plot(mod.to.plot,xlab="True Difference",ylab="Observed Difference",cex.lab=2,cex.axis=1.5,cex=0.3,main="Wth-peak diff.; Diff. error btw samples")
mod = summary(lm(mod.to.plot))
r.sq = round(mod$r.squared,digits = 2)
abline(mod,col=2,lwd=2)
mtext(text = "D ",line = -2,adj = 1,side = 1,cex = 1.5)
mtext(text = substitute(paste(" ",R^2," = ", r.sq), list(r.sq=r.sq)),line = -2,adj = 0,side = 3,cex=1.5)
rm('temp.y','temp.x','mod','r.sq')

dev.off()

} 

if (fig.made == 1) {
  
  fig.made = 0
  dev.off()
  
}
########################################################

##########
# based on assumption that a each peak has a different level of detectability across samples and they are uncorrelated across samples, but the two samples have different ranges in detectability
bray.comp = numeric()

for (i in 1:100) {
  
  #min.detect = sample(x = c(1:10),size = 1) # the minimum 'peaks.detect' value (think of it as ionization efficiency) below which we assume the peak can't be observed regardless of its concentration
  min.detect = sample(x = c(0),size = 1) # the minimum 'peaks.detect' value (think of it as ionization efficiency) below which we assume the peak can't be observed regardless of its concentration
  
  spxsite.mult = spxsite
  
  samp1.peaks.detect = runif(n = nrow(spxsite),min = 0,max = 10)
  samp2.peaks.detect = runif(n = nrow(spxsite),min = 0,max = 100)
  samp1.peaks.detect[samp1.peaks.detect < I(min.detect/10)] = 0 # this is hard coded to rescale the min.detect because the detect range is 10x lower for samp1 in this case
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

bray.comp = as.data.frame(bray.comp)
colnames(bray.comp) = c("Bray.Samp.Detect.Var","Samp1_Shannon.Samp.Detect.Var","Samp2_Shannon.Samp.Detect.Var","Samp1_Weighted.Samp.Detect.Var","Samp2_Weighted.Samp.Detect.Var")

samples.diff.detect.range = bray.comp

# based on assumption that the level of detectability varies with the trait value, but does not vary across samples
bray.comp = numeric()

for (i in 1:1) {
  
  # note that this does not include a min.detect step
  
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
print(c(niche.breadth,order.of.mag,bray.it,date(),fig.made))

}

bray.full.comp = as.data.frame(bray.full.comp)
colnames(bray.full.comp)[1:5] = c("Bray.True","Samp1.Shannon.True","Samp2.Shannon.True","Samp1.Weighted.True","Samp2.Weighted.True")

# make plots

pdf(paste0("Bray_Simulations_",peak.rich,".pdf"),height = 10)
par(pty="s",cex.lab=2,cex.axis=1.5,mfrow=c(2,1))

mod.to.plot = bray.full.comp$Bray.Peak.Var ~ bray.full.comp$Bray.True
plot(mod.to.plot,ylab="Simulated Bray-Curtis",xlab="True Bray-Curtis",main="Peaks Vary in Detectability",ylim=c(0,1),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Bray.Peak.Samp.Var ~ bray.full.comp$Bray.True
plot(mod.to.plot,ylab="Simulated Bray-Curtis",xlab="True Bray-Curtis",main="Peak Detectability Varies Across Peaks and Samples",ylim=c(0,1),cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

#mod.to.plot = bray.full.comp$Bray.Samp.Detect.Var ~ bray.full.comp$Bray.True
#plot(mod.to.plot,ylab="Simulated Bray-Curtis",xlab="True Bray-Curtis",main="Range in Peak Detectability Varies Across Samples",ylim=c(0,1),cex=0.4)
#abline(0,1,col=2,lwd=2)
#points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

#mod.to.plot = bray.full.comp$Bray.Trait.Corr ~ bray.full.comp$Bray.True
#plot(mod.to.plot,ylab="Simulated Bray-Curtis",xlab="True Bray-Curtis",main="Peak Detectability Correlates with Trait Value",ylim=c(0,1),cex=0.4)
#abline(0,1,col=2,lwd=2)
#points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

dev.off()

# make plots

pdf(paste0("Shannon_Simulations_",peak.rich,".pdf"),height = 10)
par(pty="s",cex.lab=2,cex.axis=1.5,mfrow=c(2,1))

mod.to.plot = bray.full.comp$Samp1_Shannon.Peak.Var ~ bray.full.comp$Samp1.Shannon.True
plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main="Peaks Vary in Detectability",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

#mod.to.plot = bray.full.comp$Samp2_Shannon.Peak.Var ~ bray.full.comp$Samp2.Shannon.True
#plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main="Peaks Vary in Detectability (Samp2)",cex=0.4)
#abline(0,1,col=2,lwd=2)
#points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Samp1_Shannon.Peak.Samp.Var ~ bray.full.comp$Samp1.Shannon.True
plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main="Peak Detectability Varies Across Peaks and Samples",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

#mod.to.plot = bray.full.comp$Samp2_Shannon.Peak.Samp.Var ~ bray.full.comp$Samp2.Shannon.True
#plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main="Peak Detectability Varies Across Peaks and Samples (Samp2)",cex=0.4)
#abline(0,1,col=2,lwd=2)
#points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

#mod.to.plot = bray.full.comp$Samp1_Shannon.Samp.Detect.Var ~ bray.full.comp$Samp1.Shannon.True
#plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main="Range in Peak Detectability Varies Across Samples (Samp1)",cex=0.4)
#abline(0,1,col=2,lwd=2)
#points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

#mod.to.plot = bray.full.comp$Samp2_Shannon.Samp.Detect.Var ~ bray.full.comp$Samp2.Shannon.True
#plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main="Range in Peak Detectability Varies Across Samples (Samp2)",cex=0.4)
#abline(0,1,col=2,lwd=2)
#points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

#mod.to.plot = bray.full.comp$Samp1_Shannon.Trait.Corr ~ bray.full.comp$Samp1.Shannon.True
#plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main="Peak Detectability Correlates with Trait Value (Samp1)",cex=0.4)
#abline(0,1,col=2,lwd=2)
#points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

#mod.to.plot = bray.full.comp$Samp2_Shannon.Trait.Corr ~ bray.full.comp$Samp2.Shannon.True
#plot(mod.to.plot,ylab="Simulated Shannon",xlab="True Shannon",main="Peak Detectability Correlates with Trait Value (Samp2)",cex=0.4)
#abline(0,1,col=2,lwd=2)
#points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

dev.off()

#

pdf(paste0("Weighted_Simulations_",peak.rich,".pdf"),height = 10)
par(pty="s",cex.lab=2,cex.axis=1.5,mfrow=c(2,1))

mod.to.plot = bray.full.comp$Samp1_Weighted.Peak.Var ~ bray.full.comp$Samp1.Weighted.True
plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main="Peaks Vary in Detectability",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

#mod.to.plot = bray.full.comp$Samp2_Weighted.Peak.Var ~ bray.full.comp$Samp2.Weighted.True
#plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main="Peaks Vary in Detectability (Samp2)",cex=0.4)
#abline(0,1,col=2,lwd=2)
#points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

mod.to.plot = bray.full.comp$Samp1_Weighted.Peak.Samp.Var ~ bray.full.comp$Samp1.Weighted.True
plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main="Peak Detectability Varies Across Peaks and Samples",cex=0.4)
abline(0,1,col=2,lwd=2)
points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

#mod.to.plot = bray.full.comp$Samp2_Weighted.Peak.Samp.Var ~ bray.full.comp$Samp2.Weighted.True
#plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main="Peak Detectability Varies Across Peaks and Samples (Samp2)",cex=0.4)
#abline(0,1,col=2,lwd=2)
#points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

#mod.to.plot = bray.full.comp$Samp1_Weighted.Samp.Detect.Var ~ bray.full.comp$Samp1.Weighted.True
#plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main="Range in Peak Detectability Varies Across Samples (Samp1)",cex=0.4)
#abline(0,1,col=2,lwd=2)
#points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

#mod.to.plot = bray.full.comp$Samp1_Weighted.Samp.Detect.Var ~ bray.full.comp$Samp1.Weighted.True
#plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main="Range in Peak Detectability Varies Across Samples (Samp2)",cex=0.4)
#abline(0,1,col=2,lwd=2)
#points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

#mod.to.plot = bray.full.comp$Samp1_Weighted.Trait.Corr ~ bray.full.comp$Samp1.Weighted.True
#plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main="Peak Detectability Correlates with Trait Value (Samp1)",cex=0.4)
#abline(0,1,col=2,lwd=2)
#points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

#mod.to.plot = bray.full.comp$Samp1_Weighted.Trait.Corr ~ bray.full.comp$Samp1.Weighted.True
#plot(mod.to.plot,ylab="Simulated Weighted Trait",xlab="True Weighted Trait",main="Peak Detectability Correlates with Trait Value (Samp2)",cex=0.4)
#abline(0,1,col=2,lwd=2)
#points(lowess(mod.to.plot,f = 0.3),typ="l",lwd=2,lty=2,col=5)

dev.off()

}

#######
# kernal density plot of r.sq values

# change to data frame
rsq.compiled = as.data.frame(rsq.compiled)
# add col names
colnames(rsq.compiled) = c("Error.Type",
                           "Starting.Peak.Richness",
                           "Actual.Peak.Richness",
                           "True.Zero.Removed.Peak.Richness",
                           "Niche.Breadth",
                           "Iteration",
                           "Min.Detect",
                           "Samp1.niche",
                           "Samp2.niche",
                           "Within.Peak.Diff.Rsq",
                           "Between.Peak.Diff.Rsq",
                           "Obs.vs.True.Rsq")
# changing variable types
cols.to.change = sapply(rsq.compiled,is.factor)
rsq.compiled[cols.to.change] = lapply(rsq.compiled[cols.to.change],as.character)
rsq.compiled[2:ncol(rsq.compiled)] = lapply(rsq.compiled[2:ncol(rsq.compiled)],as.numeric)

write.csv(x = rsq.compiled,file = "Rsq_Compiled.csv",quote = F,row.names = F)

# generating density functions
same.rsq.den.obs.v.true.abund = density(rsq.compiled$Obs.vs.True.Rsq[which(
  rsq.compiled$Error.Type == "Detectability.Same.Between.Samples"
)],from = 0,to = 1)

same.rsq.median.obs.v.true.abund = median(rsq.compiled$Obs.vs.True.Rsq[which(
  rsq.compiled$Error.Type == "Detectability.Same.Between.Samples"
)]) # 0.5525339

differ.rsq.den.obs.v.true.abund = density(rsq.compiled$Obs.vs.True.Rsq[which(
  rsq.compiled$Error.Type == "Detectability.Differs.Between.Samples"
)],from = 0,to = 1)

differ.rsq.median.obs.v.true.abund = median(rsq.compiled$Obs.vs.True.Rsq[which(
  rsq.compiled$Error.Type == "Detectability.Differs.Between.Samples"
)]) # 0.5574674

same.rsq.wth.obs.v.true.diff = density(rsq.compiled$Within.Peak.Diff.Rsq[which(
  rsq.compiled$Error.Type == "Detectability.Same.Between.Samples"
)],from = 0,to = 1)

same.rsq.median.wth.obs.v.true.diff = median(rsq.compiled$Within.Peak.Diff.Rsq[which(
  rsq.compiled$Error.Type == "Detectability.Same.Between.Samples"
)]) # 0.7359672

differ.rsq.wth.obs.v.true.diff = density(rsq.compiled$Within.Peak.Diff.Rsq[which(
  rsq.compiled$Error.Type == "Detectability.Differs.Between.Samples"
)],from = 0,to = 1)

differ.rsq.median.wth.obs.v.true.diff = median(rsq.compiled$Within.Peak.Diff.Rsq[which(
  rsq.compiled$Error.Type == "Detectability.Differs.Between.Samples"
)]) # 0.5329667

same.rsq.btw.obs.v.true.diff = density(rsq.compiled$Between.Peak.Diff.Rsq[which(
  rsq.compiled$Error.Type == "Detectability.Same.Between.Samples"
)],from = 0,to = 1)

same.rsq.median.btw.obs.v.true.diff = median(rsq.compiled$Between.Peak.Diff.Rsq[which(
  rsq.compiled$Error.Type == "Detectability.Same.Between.Samples"
)]) # 0.5177061

differ.rsq.btw.obs.v.true.diff = density(rsq.compiled$Between.Peak.Diff.Rsq[which(
  rsq.compiled$Error.Type == "Detectability.Differs.Between.Samples"
)],from = 0,to = 1)

differ.rsq.median.btw.obs.v.true.diff = median(rsq.compiled$Between.Peak.Diff.Rsq[which(
  rsq.compiled$Error.Type == "Detectability.Differs.Between.Samples"
)]) # 0.5232878

# make the plots
# need one panel for obs vs true abund
# one panel for within peak difference for each type of error
# one panel for between peak for one type of error (doesn't matter which)
# either 4 panels or one panel with 4 lines

pdf("Rsq_Density.pdf") # trying a single panel first

par(pty="s",cex.lab=2,cex.axis=1.5)
plot(same.rsq.den.obs.v.true.abund,typ="l",lwd=3,col=1,xlab=substitute(paste(R^2," Values")),ylab="Density",xlim=c(0,1),main="",ylim=c(0,11.5))
points(same.rsq.wth.obs.v.true.diff,typ="l",col=8,lwd=3)
points(differ.rsq.wth.obs.v.true.diff,typ="l",col=2,lwd=3)
points(same.rsq.btw.obs.v.true.diff,typ="l",col=4,lwd=3)
#points(differ.rsq.den.obs.v.true.abund,typ="l",col=3,lwd=3) # not needed
#points(differ.rsq.btw.obs.v.true.diff,typ="l",col=5,lwd=3) # not needed

legend(x = 0,y=11,legend = c("Obs. vs. True Abund.","Within-Peak Diff. (Same)","Within-Peak Diff. (Differ)","Between-Peak Diff."),col = c(1,8,2,4),lwd = 2,lty=1)

dev.off()


