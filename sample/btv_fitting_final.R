# Load libraries and source files
library(IDPmisc)
library(rjags)
library(readxl)
setwd("~/Dropbox/bluetongue/code/Final Code")
source("temp_functions.R")
source("mcmc_utils_all.R")
source("temp_deriv_functions.R")

# Set working directory and seed
bt <- read.csv("Bluetongue_Params_0723.csv")
colnames(bt)[1] <- "Trait"
set.seed(5)

# Set up MCMC params.
n.chains<-5
n.adapt<-50000
n.samps<-70000
Temps2<-seq(0,50, by=0.1)


# Subset all traits into own data frames
a <- subset(bt, bt$Trait == 'a')
b <- subset(bt, bt$Trait == 'b')
efd <- subset(bt, bt$Trait == 'EFD')
edt <- subset(bt, bt$Trait == 'EDT')
ldt <- subset(bt, bt$Trait == 'LDT')
pudt <- subset(bt, bt$Trait == 'PuDT') 
mu <- subset(bt, bt$Trait == 'mu') 
pdr <- subset(bt, bt$Trait == 'PDR')
pe <- subset(bt, bt$Trait == 'pE')
pl <- subset(bt, bt$Trait == 'pL')
pp <- subset(bt, bt$Trait == 'pP')


###################################### Mu Quad Model ####################################################
n.chains<-5
n.adapt<-5000
n.samps<-10000
mu.model.quad <- jags.model("mu-quad.bug", data = list('T' = mu$AmbientTemp, 'Y' = mu$TraitValueSI,
                                                       'N' = length(mu$AmbientTemp)), inits = list('inter' = 0.01,
                                                                                                   'n.slope' = 0.05,
                                                                                                   'qd' = 2),
                            n.chains = n.chains, n.adapt = n.adapt)
update(mu.model.quad)

coda.mu.quad <- coda.samples(mu.model.quad, c('inter', 'n.slope', 'qd', 'tau'), n.samps)

#par(mar=c(1,1,1,1))
#plot(coda.mu.quad, ask=F)

samps.mu.qd <- make.pos.quad.samps(coda.mu.quad, nchains = n.chains, samp.lims = c(1, n.samps))
samps.mu.qd$sigma <- 1/samps.mu.qd$tau
#saveRDS(samps.mu.qd, file = "mu_samps.rds")
#saveRDS(samps.mu.qd , file = "../../../Final Results/mu/mu_samps.rds")

priors.mu.qd <- list()
priors.mu.qd$names <- c("inter", "n.slope", "qd", "tau","sigma", n.samps)
priors.mu.qd$fun <- c("gamma", "gamma", "gamma", "normal", "normal")
priors.mu.qd$hyper<-matrix(NA, ncol=4, nrow=3)
priors.mu.qd$hyper[,1] <- c(1,4,NA)
priors.mu.qd$hyper[,2] <- c(4,4,NA)
priors.mu.qd$hyper[,3] <- c(2,2,NA)
priors.mu.qd$hyper[,4] <- c(400,1/200,NA)
priors.mu.qd$hyper[,5] <- c(1/400,200,NA)


par(mar = c(6.5, 6.5, 5, 5), mgp = c(5, 1, 0))
plot.hists(samps.mu.qd, my.par = c(2,2), n.hists = 5, priors = priors.mu.qd)


out.mu.qd <- make.sims.temp.resp(sim="quad.pos.trunc", samps.mu.qd, Temps2, thinned=seq(1,n.samps, length=5000))
q.mu.qd <- temp.sim.quants(out.mu.qd$fits, length(Temps2))


par(mar=c(5,5,4,1)+.1)
#dev.off()
plot(mu$AmbientTemp, mu$TraitValueSI, xlim=c(0,50), ylim=c(0,0.25),
     pch = 16,
     xlab="T (ºC)",
     ylab = "Mortality Rate (mu)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)

add.sim.lines(Temps2, sim.data=out.mu.qd$fits, q=q.mu.qd, mycol=1, lwd = 2)

# Plot mu trait with uncertainty
par(mar=c(5,5,4,1)+.1)
plot(mu$AmbientTemp, mu$TraitValueSI, xlim=c(0,50), ylim=c(0,0.25),
     pch = 16,
     xlab="T (ºC)",
     ylab="Mortality Rate (mu)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)
lines(Temps2, rowMeans(out.mu.qd$fits),lwd=4, col = 1, lty = "dashed")
for (i in seq(1,1000,100))
{lines(Temps2, out.mu.qd$fits[,i], col = "darkgrey", lwd =.5)}

#saveRDS(out.mu.qd, file = "mu_fit.rds")
#saveRDS(q.mu.qd, file = "mu_q.rds")



# Set up MCMC params.
n.chains<-5
n.adapt<-50000
n.samps<-70000
Temps2<-seq(0,50, by=0.1)


###################################### a Briere Model ##########################################################
a.fit <- jags.model('a_briere.bug', data = list('T' = a$AmbientTemp,
                                                 'Y' = a$TraitValueSI,
                                                 'N' = length(a$AmbientTemp)),inits = list('c' = 1, 'Tm' = 35, 'T0' = 3),
                      n.chains = n.chains,
                      n.adapt = n.adapt)

update(a.fit)
# Check to see if we have convergence in catepillar plots.
coda.a <- coda.samples(a.fit, c('c', 'Tm', 'T0', 'tau'), n.samps)
par(mar=c(1,1,1,1))
#plot(coda.a, ask=F)
samps.a <- make.briere.samps(coda.a, nchains=n.chains, samp.lims=c(1, n.samps))
priors.a <- list()
priors.a$names <- c("Tmin", "Tmax", "k", "tau")
priors.a$fun <- c( "uniform", "uniform","gamma", 'gamma')
priors.a$hyper <- matrix(NA, ncol=4, nrow=3)
priors.a$hyper[,1] <- c(1, 20, NA)
priors.a$hyper[,2] <- c(20, 40, NA)
priors.a$hyper[,3] <- c(1, 20, NA)
priors.a$hyper[,4] <- c(0.01, 0.01, NA)

## Plot histograms
par(mar = c(6.5, 6.5, 1.5, 1.5), mgp = c(5, 1, 0))
plot.hists(samps.a, my.par=c(2,2), n.hists=4, priors=priors.a)

out.a <- make.sims.temp.resp(sim="briere", samps.a, Temps2, thinned=seq(1,n.samps, length=5000))
q.a <- temp.sim.quants(out.a$fits, length(Temps2))

# Plot a trait
par(mar=c(5,5,4,1)+.1)
#dev.off()
plot(a$AmbientTemp, a$TraitValueSI, xlim=c(0,50), ylim=c(0,1),
     pch = 16,
     xlab="T (ºC)",
     ylab = "Biting Rate (a)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)

add.sim.lines(Temps2, sim.data=out.a$fits, q=q.a, mycol=1, lwd = 2)

# Plot a trait with uncertainty
par(mfrow=c(1,1),mar=c(5,5,4,1)+.1)
#dev.off()
plot(a$AmbientTemp, a$TraitValueSI, xlim=c(0,50), ylim=c(0,1),
     pch = 16,
     xlab="T (ºC)",
     ylab="Biting Rate (a)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis =2)
lines(Temps2, rowMeans(out.a$fits),lwd=4, col = 1, lty = "dashed")
for (i in seq(1,1000,100))
{lines(Temps2, out.a$fits[,i], col = "darkgray", lwd =.5)}

#To save
saveRDS(out.a, file = "a_fit.rds")
saveRDS(q.a, file = "a_q.rds")



###################################### EFD Briere Model ##########################################################
efd.fit <- jags.model('EFD_briere.bug',
                     data = list('Y' = efd$TraitValueSI,'T' = efd$AmbientTemp,
                                 'N' = length(efd$AmbientTemp)),
                     inits =list('c' = 100,'Tm' = 31, 'T0' = 1.2),
                     n.chains = n.chains,
                     n.adapt = n.adapt)

update(efd.fit)
efd.coda <- coda.samples(efd.fit, c('c','Tm','T0', 'tau'), n.samps)
#par(mar=c(1,1,1,1))
#plot(efd.coda)

efd.samps <- make.briere.samps(efd.coda, nchains=n.chains, samp.lims=c(1, n.samps))
#names(efd.samps) <- c("Tmin", "Tmax", "k", "tau")


priors.efd <- list()
priors.efd$names <- c("c", "Tm","T0", "tau")
priors.efd$fun <- c( "gamma", "uniform","uniform", "gamma")
priors.efd$hyper <- matrix(NA, ncol=4, nrow=3)
priors.efd$hyper[,1] <- c(1, 1, NA)
priors.efd$hyper[,2] <- c(29, 35, NA)
priors.efd$hyper[,3] <- c(1, 2, NA)
priors.efd$hyper[,4] <- c(9, 0.0005, NA)

# par(mar = c(6.5, 6.5, 1.5, 1.5), mgp = c(5, 1, 0))
# plot.hists(efd.samps[,c(1:4)], my.par=c(2,2), n.hists=4, priors=priors)

efd.temps <- seq(0,50, by=0.1)
efd.out <- make.sims.temp.resp(sim="briere", efd.samps, efd.temps, 
                              thinned=seq(1,n.samps,length=5000))
efd.q <- temp.sim.quants(efd.out$fits, length(efd.temps))

# Plot a trait
par(mar=c(5,5,4,1)+.1)
#dev.off()
plot(efd$AmbientTemp, efd$TraitValueSI, xlim=c(0,50), ylim=c(0,80),
     pch = 16,
     xlab="T (ºC)",
     ylab = "Eggs per Female per Day (EFD)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)

add.sim.lines(Temps2, sim.data=efd.out$fits, q=efd.q, mycol=1, lwd = 2)

# Plot efd trait with uncertainty
par(mfrow=c(1,1),mar=c(5,5,4,1)+.1)
plot(efd$AmbientTemp, efd$TraitValueSI, xlim=c(0,50), ylim=c(0,80),
     pch = 16,
     xlab="T (ºC)",
     ylab="Eggs per Female per Day (EFD)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis =2)
lines(efd.temps, rowMeans(efd.out$fits),lwd=4, col = 1, lty = "dashed")
for (i in seq(1,1000,100))
{lines(efd.temps, efd.out$fits[,i], col = "darkgray", lwd =.5)}

saveRDS(efd.out, file = "efd_fit.rds")
saveRDS(efd.q, file = "efd_q.rds")

###################################### EDT Quad Case ####################################################
edt.model.quad <- jags.model("edt-quad.bug", data = list('T' = edt$AmbientTemp, 'Y' = edt$TraitValueSI,
                                                         'N' = length(edt$AmbientTemp)),
                             inits = list('inter' = 35,
                                          'n.slope' = 1,
                                          'qd' = 0.02),
                             n.chains = n.chains, n.adapt = n.adapt)

update(edt.model.quad)
coda.edt.quad <- coda.samples(edt.model.quad, c('inter', 'n.slope', 'qd', 'tau'), thinned=seq(1,75000, length=10), 75000)

samps.edt.qd <- make.pos.quad.samps(coda.edt.quad, nchains = n.chains, samp.lims = c(1,n.samps))
names(samps.edt.qd)[4] = "tau"
#par(mar=c(1,1,1,1))
#plot(coda.edt.quad, ask=F)

priors.edt.qd <- list()
priors.edt.qd$names <- c("inter", "n.slope", "qd", "tau", n.samps)
priors.edt.qd$fun <- c("gamma", "gamma", "gamma", "normal")
priors.edt.qd$hyper<-matrix(NA, ncol=4, nrow=3)
priors.edt.qd$hyper[,1] <- c(1,0.01,NA)
priors.edt.qd$hyper[,2] <- c(1,0.5,NA)
priors.edt.qd$hyper[,3] <- c(4,28,NA)
priors.edt.qd$hyper[,4] <- c(3,1/800,NA)

#source("mcmc_utils_all.R") # Adjust limits for aesthetics
#par(mar = c(6.5, 6.5, 1.5, 1.5), mgp = c(5, 1, 0))
#plot.hists(samps.edt.qd, my.par = c(2,2), n.hists = 4, priors = priors.edt.qd)

source("mcmc_utils_all.R")
source("temp_deriv_functions.R")
out.edt.qd <- make.sims.temp.resp(sim="quad.pos.trunc", samps.edt.qd, Temps2, thinned=seq(1,n.samps, length=5000))
q.edt.qd <- temp.sim.quants(out.edt.qd$fits, length(Temps2))

par(mar=c(5,5,4,1)+.1)
#dev.off()
plot(edt$AmbientTemp, edt$TraitValueSI, xlim=c(0,50), ylim=c(0,150),
     pch = 16,
     xlab="T (ºC)",
     ylab = "Eggs Development Time (EDT)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)

add.sim.lines(Temps2, sim.data=out.edt.qd$fits, q=q.edt.qd, mycol=1, lwd = 2)


# Plot edt trait with uncertainty
par(mar=c(5,5,4,1)+.1)
plot(edt$AmbientTemp, edt$TraitValueSI, xlim=c(0,50), ylim=c(0,150),
     pch = 16,
     xlab="T (ºC)",
     ylab="Eggs Development Time (EDT)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)
lines(Temps2, rowMeans(out.edt.qd$fits),lwd=4, col = 1, lty = "dashed")
for (i in seq(1,1000,100))
{lines(Temps2, out.edt.qd$fits[,i], col = "darkgrey", lwd =.5)}

saveRDS(out.edt.qd, file = "edt_fit.rds")
saveRDS(q.edt.qd, file = "edt_q.rds")


###################################### LDT Quad Case ####################################################
ldt.model.quad <- jags.model("ldt-quad.bug", data = list('T' = ldt$AmbientTemp, 'Y' = ldt$TraitValueSI,
                                                         'N' = length(ldt$AmbientTemp)),
                             inits = list('inter' = 35,
                                          'n.slope' = 1,
                                          'qd' = 0.02),
                             n.chains = n.chains, n.adapt = n.adapt)

update(ldt.model.quad)
coda.ldt.quad <- coda.samples(ldt.model.quad, c('inter', 'n.slope', 'qd', 'tau'),thinned=seq(1,75000, length=150), 75000)

#par(mar=c(1,1,1,1))
#plot(coda.ldt.quad, ask=F)

samps.ldt.qd <- make.pos.quad.samps(coda.ldt.quad, nchains = n.chains, samp.lims = c(1,n.samps))
names(samps.ldt.qd)[4] = "tau"
saveRDS(samps.ldt.qd, file = "ldt_samps.rds")

priors.ldt.qd <- list()
priors.ldt.qd$names <- c("inter", "n.slope", "qd", "tau", n.samps)
priors.ldt.qd$fun <- c("gamma", "gamma", "gamma", "normal")
priors.ldt.qd$hyper<-matrix(NA, ncol=4, nrow=3)
priors.ldt.qd$hyper[,1] <- c(2,2,NA)
priors.ldt.qd$hyper[,2] <- c(3,3,NA)
priors.ldt.qd$hyper[,3] <- c(2,2,NA)
priors.ldt.qd$hyper[,4] <- c(1000,1/500,NA)

#par(mar = c(6.5, 6.5, 1.5, 1.5), mgp = c(5, 1, 0))
#plot.hists(samps.ldt.qd, my.par = c(2,2), n.hists = 4, priors = priors.ldt.qd)

out.ldt.qd <- make.sims.temp.resp(sim="quad.pos.trunc", samps.ldt.qd, Temps2, thinned=seq(1,n.samps, length=5000))
q.ldt.qd <- temp.sim.quants(out.ldt.qd$fits, length(Temps2))

par(mar=c(5,5,4,1)+.1)
#dev.off()
plot(ldt$AmbientTemp, ldt$TraitValueSI, xlim=c(0,50), ylim=c(0,60),
     pch = 16,
     xlab="T (ºC)",
     ylab = "Larva Development Time (LDT)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)

add.sim.lines(Temps2, sim.data=out.ldt.qd$fits, q=q.ldt.qd, mycol=1, lwd = 2)


# Plot a trait with uncertainty
par(mar=c(5,5,4,1)+.1)
plot(ldt$AmbientTemp, ldt$TraitValueSI, xlim=c(0,50), ylim=c(0,60),
     pch = 16,
     xlab="T (ºC)",
     ylab="Larva Development Time (LDT)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)
lines(Temps2, rowMeans(out.ldt.qd$fits),lwd=4, col = 1, lty = "dashed")
for (i in seq(1,1000,100))
{lines(Temps2, out.ldt.qd$fits[,i], col = "darkgrey", lwd =.5)}

saveRDS(out.ldt.qd, file = "ldt_fit.rds")
saveRDS(q.ldt.qd, file = "ldt_q.rds")


###################################### PuDT Quad Case ####################################################
pudt.model.quad <- jags.model("pudt-quad.bug", data = list('T' = pudt$AmbientTemp, 'Y' = pudt$TraitValueSI,
                                                           'N' = length(pudt$AmbientTemp)),
                              inits = list('inter' = 20,
                                           'n.slope' = 5,
                                           'qd' = 0.2),
                              n.chains = n.chains, n.adapt = n.adapt)

update(pudt.model.quad)
coda.pudt.quad <- coda.samples(pudt.model.quad, c('inter', 'n.slope', 'qd', 'tau'), thinned=seq(1,100000, length=500), 100000)

#par(mar=c(1,1,1,1))
#plot(coda.pudt.quad, ask=F)

samps.pudt.qd <- make.pos.quad.samps(coda.pudt.quad, nchains = n.chains, samp.lims = c(1,n.samps))
names(samps.pudt.qd)[4] = "tau"
saveRDS(samps.pudt.qd, file = "pudt_samps.rds")

priors.pudt.qd <- list()
priors.pudt.qd$names <- c("inter", "n.slope", "qd", "tau", n.samps)
priors.pudt.qd$fun <- c("gamma", "gamma", "gamma", "normal")
priors.pudt.qd$hyper<-matrix(NA, ncol=4, nrow=3)
priors.pudt.qd$hyper[,1] <- c(1,0.01,NA)
priors.pudt.qd$hyper[,2] <- c(1,0.5,NA)
priors.pudt.qd$hyper[,3] <- c(4,28,NA)
priors.pudt.qd$hyper[,4] <- c(3,1/200,NA)

#par(mar = c(6.5, 6.5, 1.5, 1.5), mgp = c(5, 1, 0))
#plot.hists(samps.pudt.qd, my.par = c(2,2), n.hists = 4, priors = priors.pudt.qd)

out.pudt.qd <- make.sims.temp.resp(sim="quad.pos.trunc", samps.pudt.qd, Temps2, thinned=seq(1,n.samps, length=5000))
q.pudt.qd <- temp.sim.quants(out.pudt.qd$fits, length(Temps2))


par(mar=c(5,5,4,1)+.1)
#dev.off()
plot(pudt$AmbientTemp, pudt$TraitValueSI, xlim=c(0,50), ylim=c(0,150),
     pch = 16,
     xlab="T (ºC)",
     ylab = "Pupa Development Time (PuDT)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)

add.sim.lines(Temps2, sim.data=out.pudt.qd$fits, q=q.pudt.qd, mycol=1, lwd = 2)


# Plot a trait with uncertainty
par(mar=c(5,5,4,1)+.1)
plot(pudt$AmbientTemp, pudt$TraitValueSI, xlim=c(0,50), ylim=c(0,150),
     pch = 16,
     xlab="T (ºC)",
     ylab="Pupa Development Time (PuDT)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)
lines(Temps2, rowMeans(out.pudt.qd$fits),lwd=4, col = 1, lty = "dashed")
for (i in seq(1,1000,100))
{lines(Temps2, out.pudt.qd$fits[,i], col = "darkgrey", lwd =.5)}

saveRDS(out.pudt.qd, file = "pudt_fit.rds")
saveRDS(q.pudt.qd, file = "pudt_q.rds")



# New temp functinos for these traits
source("temp_functions.R")


###################################### b Briere (Binomial) Model ##########################################################
b.fit <- jags.model('b_brierebinom.bug',
                    data = list('Y' = b$TraitValueSI*25,'T' = b$AmbientTemp,
                                'N'=length(b$AmbientTemp),
                                'n'=rep(25, length(b[,1]))),
                    inits =list('c' = 0.0008,'Tm' = 35, 'T0' = 12),
                    n.chains = n.chains,
                    n.adapt = n.adapt)

update(b.fit)
b.coda <- coda.samples(b.fit, c('c','Tm','T0', 'tau'), n.samps)
#par(mar=c(1,1,1,1))
#plot(b.coda)

b.samps <- make.briere.samps(b.coda, nchains=n.chains, samp.lims=c(1, n.samps))
#names(b.samps) = c("Tmin", "Tmax", "k", "tau")



priors.b <- list()
priors.b$names <- c("c", "Tm", "T0", "sigma")
priors.b$fun <- c("gamma", "uniform", "uniform", "gamma")
priors.b$hyper <- matrix(NA, ncol=4, nrow=3)
priors.b$hyper[,1] <- c(1, 10, NA)
priors.b$hyper[,2] <- c(25, 35, NA)
priors.b$hyper[,3] <- c(10, 24, NA)
priors.b$hyper[,4] <- c(5.5, 0.025, NA)

#par(mar = c(6.5, 6.5, 1.5, 1.5), mgp = c(5, 1, 0))
#plot.hists(b.samps[,], my.par=c(2,2), n.hists=3, priors=priors.b)

b.temps <- seq(0,50, by=0.1)
b.out <- make.sims.temp.resp(sim="briere", b.samps, b.temps, 
                             thinned=seq(1,n.samps,length=5000))
b.q <- temp.sim.quants(b.out$fits, length(b.temps))


par(mar=c(5,5,4,1)+.1)
#dev.off()
plot(b$AmbientTemp, b$TraitValueSI, xlim=c(0,50), ylim=c(0,1.2),
     pch = 16,
     xlab="T (ºC)",
     ylab = "Vector Competence (b)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)

add.sim.lines(b.temps, sim.data=b.out$fits, q=b.q, mycol=1, lwd = 2)


# Plot b trait with uncertainty
par(mar=c(5,5,4,1)+.1)
plot(b$AmbientTemp, b$TraitValueSI, xlim=c(0,50), ylim=c(0,1.2),
     pch = 16,
     xlab="T (ºC)",
     ylab="Vector Competence (b)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)
lines(b.temps, rowMeans(b.out$fits),lwd=4, col = 1, lty = "dashed")
for (i in seq(1,1000,100))
{lines(b.temps, b.out$fits[,i], col = "darkgray", lwd =.5)}

saveRDS(b.out, file = "b_fit.rds")
saveRDS(b.q, file = "b_q.rds")


###################################### pE Briere Model ####################################################
pe.fit <- jags.model('pE_briere.bug',
                    data = list('Y' = pe$TraitValueSI,'T' = pe$AmbientTemp,
                                'N'=length(pe$AmbientTemp)),
                    inits =list('c' = 0.1,'Tm' = 35,'T0' = 5),
                    n.chains = n.chains,
                    n.adapt = n.adapt)

pe.coda <- coda.samples(pe.fit, c('c','Tm','T0', 'tau'), n.samps)
#par(mar=c(1,1,1,1))
#plot(pe.coda)

pe.samps <- make.briere.samps(pe.coda, nchains=n.chains, samp.lims=c(1, n.samps))
## Rename histogram titles
# names(pe.samps)[1] = "Tmin"
# names(pe.samps)[2] = "Tmax"
# names(pe.samps)[3] = "k"
# names(pe.samps)[4] = "tau"
#saveRDS(pe.samps, file = "pE_samps.rds")

priors.pe <- list()
priors.pe$names <- c("c", "Tm","T0", 'tau')
priors.pe$fun <- c( "gamma", "uniform","uniform", "normal")
priors.pe$hyper <- matrix(NA, ncol=4, nrow=3)
priors.pe$hyper[,1] <- c(1,1/5, NA)
priors.pe$hyper[,2] <- c(34, 40, NA)
priors.pe$hyper[,3] <- c(0, 10, NA)
priors.pe$hyper[,4] <- c(7, 0.000000005, NA)

# source("mcmc_utils_all.R") # Here, we toggle with ylim for a few histograms for aesthetics
# par(mar = c(6.5, 6.5, 1.5, 1.5), mgp = c(5, 1, 0))
# plot.hists(pe.samps[,c(1:4)], my.par=c(2,2), n.hists=4, priors=priors)

pe.temps <- seq(0,50, by=0.1)
pe.out <- make.sims.temp.resp(sim="briere", pe.samps, pe.temps, 
                             thinned=seq(1,n.samps,length=5000))
pe.q <- temp.sim.quants(pe.out$fits, length(pe.temps))


par(mar=c(5,5,4,1)+.1)
#dev.off()
plot(pe$AmbientTemp, pe$TraitValueSI, xlim=c(0,50), ylim=c(0,1),
     pch = 16,
     xlab="T (ºC)",
     ylab = "Eggs Survival Probability (pE)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)

add.sim.lines(pe.temps, sim.data=pe.out$fits, q=pe.q, mycol=1, lwd = 2)

# Plot pe trait with uncertainty
par(mar=c(5,5,4,1)+.1)
plot(pe$AmbientTemp, pe$TraitValueSI, xlim=c(0,50), ylim=c(0,1),
     pch = 16,
     xlab="T (ºC)",
     ylab="Eggs Survival Probability (pE)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)
lines(pe.temps, rowMeans(pe.out$fits),lwd=4, col = 1, lty = "dashed")
for (i in seq(1,1000,100))
{lines(pe.temps, pe.out$fits[,i], col = "darkgrey", lwd =.5)}

saveRDS(pe.out, file = "pE_fit.rds")
saveRDS(pe.q, file = "pE_q.rds")

###################################### pL Briere Model ####################################################
pl.fit <- jags.model('pL_briere.bug',
                    data = list('Y' = pl$TraitValueSI,'T' = pl$AmbientTemp,
                                'N'=length(pl$AmbientTemp)),
                    inits =list('c' = 0.01,'Tm' = 35,'T0' = 5),
                    n.chains = n.chains,
                    n.adapt = n.adapt)

pl.coda <- coda.samples(pl.fit, c('c','Tm','T0', 'tau'), n.samps)
#par(mar=c(1,1,1,1))
#plot(pl.coda)

pl.samps <- make.briere.samps(pl.coda, nchains=n.chains, samp.lims=c(1, n.samps))
## Rename histogram titles
# names(pl.samps)[1] = "Tmin"
# names(pl.samps)[2] = "Tmax"
# names(pl.samps)[3] = "k"
# names(pl.samps)[4] = "tau"
saveRDS(pl.samps, file = "pL_samps.rds")


priors.pl <- list()
priors.pl$names <- c("c", "Tm","T0", "tau")
priors.pl$fun <- c( "gamma", "uniform","uniform", "gamma")
priors.pl$hyper = matrix(NA, ncol=4, nrow=3)
priors.pl$hyper[,1] <- c(1,1/2, NA)
priors.pl$hyper[,2] <- c(30, 40, NA)
priors.pl$hyper[,3] <- c(0, 8, NA)
priors.pl$hyper[,4] <- c(1.5, 0.001, NA)

# source("mcmc_utils_all.R")
# par(mar = c(6.5, 6.5, 1.5, 1.5), mgp = c(5, 1, 0))
# plot.hists(pl.samps[,c(1:4)], my.par=c(2,2), n.hists=4, priors=priors)

pl.temps <- seq(0,50, by=0.1)
pl.out <- make.sims.temp.resp(sim="briere", pl.samps, pl.temps, 
                             thinned=seq(1,n.samps,length=5000))
pl.q <- temp.sim.quants(pl.out$fits, length(pl.temps))

par(mar=c(5,5,4,1)+.1)
#dev.off()
plot(pl$AmbientTemp, pl$TraitValueSI, xlim=c(0,50), ylim=c(0,1),
     pch = 16,
     xlab="T (ºC)",
     ylab = "Larva Survival Probability (pL)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)

add.sim.lines(pl.temps, sim.data=pl.out$fits, q=pl.q, mycol=1, lwd = 2)

# Plot b trait with uncertainty
par(mar=c(5,5,4,1)+.1)
plot(pl$AmbientTemp, pl$TraitValueSI, xlim=c(0,50), ylim=c(0,0.3),
     pch = 16,
     xlab="T (ºC)",
     ylab="Larva Survival Probability (pL)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)
lines(pl.temps, rowMeans(pl.out$fits),lwd=4, col = 1, lty = "dashed")
for (i in seq(1,1000,100))
{lines(pl.temps, pl.out$fits[,i], col = "darkgrey", lwd =.5)}

saveRDS(pl.out, file = "pL_fit.rds")
saveRDS(pl.q, file = "pL_q.rds")

###################################### pP Briere Model ####################################################
pp.fit <- jags.model('pP_briere.bug',
                    data = list('Y' = pp$TraitValueSI,'T' = pp$AmbientTemp,
                                'N'=length(pp$AmbientTemp)),
                    inits =list('c' = 0.01,'Tm' =37.5,'T0' = 1.5),
                    n.chains = n.chains,
                    n.adapt = n.adapt)

pp.coda <- coda.samples(pp.fit, c('c','Tm','T0', 'tau'), n.samps)
#par(mar=c(1,1,1,1))
#plot(pp.coda)

pp.samps <- make.briere.samps(pp.coda, nchains=n.chains, samp.lims=c(1, n.samps))
## Rename histogram titles
# names(pp.samps)[1] = "Tmin"
# names(pp.samps)[2] = "Tmax"
# names(pp.samps)[3] = "k"
# names(pp.samps)[4] = "tau"
saveRDS(pp.samps, file = "pP_samps.rds")

priors.pp <- list()
priors.pp$names <- c("c", "Tm","T0", "tau")
priors.pp$fun <- c( "gamma", "uniform","uniform", "gamma")
priors.pp$hyper <- matrix(NA, ncol=4, nrow=3)
priors.pp$hyper[,1] <- c(1,5, NA)
priors.pp$hyper[,2] <- c(35, 40, NA)
priors.pp$hyper[,3] <- c(1, 5, NA)
priors.pp$hyper[,4] <- c(10, 0.002, NA)

# source("mcmc_utils_all.R")
# par(mar = c(6.5, 6.5, 1.5, 1.5), mgp = c(5, 1, 0))
# plot.hists(pp.samps[,c(1:4)], my.par=c(2,2), n.hists=4, priors=priors)

pp.temps <- seq(0,50, by=0.1)
pp.out <- make.sims.temp.resp(sim="briere.trunc", pp.samps, pp.temps, 
                             thinned=seq(1,n.samps,length=5000))
pp.q <- temp.sim.quants(pp.out$fits, length(pp.temps)) 

par(mar=c(5,5,4,1)+.1)
#dev.off()
plot(pp$AmbientTemp, pp$TraitValueSI, xlim=c(0,50), ylim=c(0,1),
     pch = 16,
     xlab="T (ºC)",
     ylab = "Pupa Survival Probability (pP)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)

add.sim.lines(pp.temps, sim.data=pp.out$fits, q=pp.q, mycol=1, lwd = 2)


# Plot pP trait with uncertainty
par(mar=c(5,5,4,1)+.1)
plot(pp$AmbientTemp, pp$TraitValueSI, xlim=c(0,50), ylim=c(0,1.2),
     pch = 16,
     xlab="T (ºC)",
     ylab="Pupa Survival Probability (pP)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)
lines(pp.temps, rowMeans(pp.out$fits),lwd=4, col = 1, lty = "dashed")
for (i in seq(1,1000,100))
{lines(pp.temps, pp.out$fits[,i], col = "darkgrey", lwd =.5)}

saveRDS(pp.out, file = "pP_fit.rds")
saveRDS(pp.q, file = "pP_q.rds")


###################################### PDR Briere Model ####################################################
pdr.fit <- jags.model('PDR_briere.bug',
                     data = list('Y'=pdr$TraitValueSI,'T'=pdr$AmbientTemp,
                                 'N'=length(pdr$AmbientTemp)),
                     inits =list('c' = 3,'Tm' = 35, 'T0' = 15),
                     n.chains = n.chains,
                     n.adapt = n.adapt)

pdr.coda <- coda.samples(pdr.fit, c('c','Tm','T0', 'tau'), n.samps)
#par(mar=c(1,1,1,1))
#plot(pdr.coda)

pdr.samps <- make.briere.samps(pdr.coda, nchains=n.chains, samp.lims=c(1, n.samps))
## Rename histrogram titles
# names(pdr.samps)[1] = "Tmin"
# names(pdr.samps)[2] = "Tmax"
# names(pdr.samps)[3] = "k"
# names(pdr.samps)[4] = "tau"
#saveRDS(pdr.samps, file = "pdr_samps.rds")


priors.pdr <- list()
priors.pdr$names <- c("c", "Tm","T0", "tau")
priors.pdr$fun <- c( "gamma", "uniform","uniform", "gamma")
priors.pdr$hyper <- matrix(NA, ncol=4, nrow=3)
priors.pdr$hyper[,1] <- c(1, 10, NA)
priors.pdr$hyper[,2] <- c(18, 45, NA)
priors.pdr$hyper[,3] <- c(1, 17, NA)
priors.pdr$hyper[,4] <- c(9, 0.05, NA)

# par(mar = c(6.5, 6.5, 1.5, 1.5), mgp = c(5, 1, 0))
# source("mcmc_utils_all.R")
# plot.hists(pdr.samps[,c(1:4)], my.par=c(2,2), n.hists=4, priors=priors)

pdr.temps <- seq(0,50, by=0.1)
pdr.out <- make.sims.temp.resp(sim="briere", pdr.samps, pdr.temps, 
                              thinned=seq(1,n.samps,length=5000))
pdr.q <- temp.sim.quants(pdr.out$fits, length(pdr.temps))


par(mar=c(5,5,4,1)+.1)
#dev.off()
plot(pdr$AmbientTemp, pdr$TraitValueSI, xlim=c(0,50), ylim=c(0,0.35),
     pch = 16,
     xlab="T (ºC)",
     ylab = "Parasite Development Rate (P)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)

add.sim.lines(pdr.temps, sim.data=pdr.out$fits, q=pdr.q, mycol=1, lwd = 2)

# Plot pdr trait with uncertainty
par(mar=c(5,5,4,1)+.1)
plot(pdr$AmbientTemp, pdr$TraitValueSI, xlim=c(0,50), ylim=c(0,0.35),
     pch = 16,
     xlab="T (ºC)",
     ylab="Parasite Development Rate (P)",
     col="black", cex=1.5, lwd=2, cex.lab=2, cex.axis = 2)
lines(pdr.temps, rowMeans(pdr.out$fits),lwd=4, col = 1, lty = "dashed")
for (i in seq(1,1000,100))
{lines(pdr.temps, pdr.out$fits[,i], col = "darkgrey", lwd =.5)}

saveRDS(pdr.out, file = "pdr_fit.rds")
saveRDS(pdr.q, file = "pdr_q.rds")

