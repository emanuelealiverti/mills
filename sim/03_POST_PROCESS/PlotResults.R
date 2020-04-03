rm(list=ls())
emp_kl = function(x,y) sum(x*(log(x) - log(y)),na.rm = T)

cat("+++++++++++++++++++",'\n')
cat("POSTERIOR MEAN",'\n')
cat("+++++++++++++++++++",'\n')

Y = read.table("../02_COMPETITORS/data1.txt")
for(j in 1:NCOL(Y)) Y[,j] = factor(Y[,j], levels = 1:4, labels = 1:4)
p = NCOL(Y)
d=4
H = 5
all_pairs = combn(p,2)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Now compare goodness of fit (kl divergence, pearson residuals, wasserstein)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
load("../01_MILLS/SIM1_biv.RData")
load("../02_COMPETITORS/SIMPLEX_FACTOR/SIM1_SF_biv.RData")
load("../02_COMPETITORS/LATENT_CLASS//SIM1_LC_biv.RData")
## Empirical
emp_biv = apply(all_pairs, 2, function(x) 1/NROW(Y)*table(Y[,x[1]],Y[,x[2]])) 

#+++++++++++++++++
# Person residuals
#+++++++++++++++++
pears_mills = colSums(((apply(biv_pmean, 1, c) - emp_biv)^2) / apply(biv_pmean, 1, c))
pears_lc = colSums( ( (apply(biv_pmean_LC, 1, c) - emp_biv)^2 ) / apply(biv_pmean_LC, 1, c) )
pears_sf = colSums( ( (apply(biv_pmean_SF, 1, c) - emp_biv)^2 ) / apply(biv_pmean_SF, 1, c) )


#++++++++++++
# kl
#++++++++++++
kl_sf = kl_lc = numeric(NCOL(all_pairs))
kl_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
	kl_mills[k] = emp_kl( emp_biv[,k], c(biv_pmean[k,,]))
	kl_lc[k] = emp_kl( emp_biv[,k], c(biv_pmean_LC[k,,]))
	kl_sf[k] = emp_kl( emp_biv[,k], c(biv_pmean_SF[k,,])) 
}
#++++++++++++
# Wasserstein
#++++++++++++
require(transport)

wass_sf = wass_lc = numeric(NCOL(all_pairs))
wass_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
p1 = pgrid(mass = matrix(emp_biv[,k],4,4))
p2 = pgrid(mass =biv_pmean[k,,])
wass_mills[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_pmean_LC[k,,])
wass_lc[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_pmean_SF[k,,])
wass_sf[k] = wasserstein(p1,p2)

}

dfS1 = list("PEARSON" = list("EMPIRICAL" = data.frame("MILLS" = sqrt(pears_mills/3),"LATENT CLASS"= sqrt(pears_lc/3), "TUCKER"= sqrt(pears_sf/3))),
	  "KULLBACK–LEIBLER" = list("EMPIRICAL" = data.frame("MILLS" = kl_mills,"LATENT CLASS"= kl_lc,"TUCKER"= kl_sf)),
	  "WASSERSTEIN" = list("EMPIRICAL" = data.frame("MILLS" = wass_mills,"LATENT CLASS"= wass_lc,"TUCKER"= wass_sf)))

#+++++++
# SECOND
#+++++++
cat("SECOND SCENARIO", "\n")
Y = read.table("../02_COMPETITORS/data2.txt")
for(j in 1:NCOL(Y)) Y[,j] = factor(Y[,j], levels = 1:4, labels = 1:4)
p = NCOL(Y)
d=4
H = 5
load("../01_MILLS/SIM2_biv.RData")
load("../02_COMPETITORS/SIMPLEX_FACTOR/SIM2_SF_biv.RData")
load("../02_COMPETITORS/LATENT_CLASS//SIM2_LC_biv.RData")
## Empirical
emp_biv = apply(all_pairs, 2, function(x) 1/NROW(Y)*table(Y[,x[1]],Y[,x[2]])) 

#+++++++++++++++++
# Person residuals
#+++++++++++++++++
pears_mills = colSums(((apply(biv_pmean, 1, c) - emp_biv)^2) / apply(biv_pmean, 1, c))
pears_lc = colSums( ( (apply(biv_pmean_LC, 1, c) - emp_biv)^2 ) / apply(biv_pmean_LC, 1, c) )
pears_sf = colSums( ( (apply(biv_pmean_SF, 1, c) - emp_biv)^2 ) / apply(biv_pmean_SF, 1, c) )

#++++++++++++
# kl
#++++++++++++
kl_sf = kl_lc = numeric(NCOL(all_pairs))
kl_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
	kl_mills[k] = emp_kl( emp_biv[,k], c(biv_pmean[k,,]))
	kl_lc[k] = emp_kl( emp_biv[,k], c(biv_pmean_LC[k,,]))
	kl_sf[k] = emp_kl( emp_biv[,k], c(biv_pmean_SF[k,,]) )
}
#++++++++++++
# Wasserstein
#++++++++++++
require(transport)

wass_sf = wass_lc = numeric(NCOL(all_pairs))
wass_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
p1 = pgrid(mass = matrix(emp_biv[,k],4,4))
p2 = pgrid(mass =biv_pmean[k,,])
wass_mills[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_pmean_LC[k,,])
wass_lc[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_pmean_SF[k,,])
wass_sf[k] = wasserstein(p1,p2)

}

dfS2 = list("PEARSON" = list("EMPIRICAL" = data.frame("MILLS" = sqrt(pears_mills/3),"LATENT CLASS"= sqrt(pears_lc/3), "TUCKER"= sqrt(pears_sf/3))),
	  "KULLBACK–LEIBLER" = list("EMPIRICAL" = data.frame("MILLS" = kl_mills,"LATENT CLASS"= kl_lc,"TUCKER"= kl_sf)),
	  "WASSERSTEIN" = list("EMPIRICAL" = data.frame("MILLS" = wass_mills,"LATENT CLASS"= wass_lc,"TUCKER"= wass_sf)))

#+++++++
# THIRD
#+++++++
cat("THIRD SCENARIO", "\n")
Y = read.table("../02_COMPETITORS/data3.txt")
for(j in 1:NCOL(Y)) Y[,j] = factor(Y[,j], levels = 1:4, labels = 1:4)
all_pairs = combn(p,2)

load("../01_MILLS/SIM3_biv.RData")
load("../02_COMPETITORS/SIMPLEX_FACTOR/SIM3_SF_biv.RData")
load("../02_COMPETITORS/LATENT_CLASS//SIM3_LC_biv.RData")
## Empirical
emp_biv = apply(all_pairs, 2, function(x) 1/NROW(Y)*table(Y[,x[1]],Y[,x[2]])) 

#+++++++++++++++++
# Person residuals
#+++++++++++++++++
pears_mills = colSums(((apply(biv_pmean, 1, c) - emp_biv)^2) / apply(biv_pmean, 1, c))
pears_lc = colSums( ( (apply(biv_pmean_LC, 1, c) - emp_biv)^2 ) / apply(biv_pmean_LC, 1, c) )
pears_sf = colSums( ( (apply(biv_pmean_SF, 1, c) - emp_biv)^2 ) / apply(biv_pmean_SF, 1, c) )

#++++++++++++
# kl
#++++++++++++
kl_sf = kl_lc = numeric(NCOL(all_pairs))
kl_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
	kl_mills[k] = emp_kl( emp_biv[,k], c(biv_pmean[k,,]))
	kl_lc[k] = emp_kl( emp_biv[,k], c(biv_pmean_LC[k,,]))
	kl_sf[k] = emp_kl( emp_biv[,k], c(biv_pmean_SF[k,,]) )
}
#++++++++++++
# Wasserstein
#++++++++++++
require(transport)

wass_sf = wass_lc = numeric(NCOL(all_pairs))
wass_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
p1 = pgrid(mass = matrix(emp_biv[,k],4,4))
p2 = pgrid(mass =biv_pmean[k,,])
wass_mills[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_pmean_LC[k,,])
wass_lc[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_pmean_SF[k,,])
wass_sf[k] = wasserstein(p1,p2)

}

dfS3 = list("PEARSON" = list("EMPIRICAL" = data.frame("MILLS" = sqrt(pears_mills/3),"LATENT CLASS"= sqrt(pears_lc/3), "TUCKER"= sqrt(pears_sf/3))),
	  "KULLBACK–LEIBLER" = list("EMPIRICAL" = data.frame("MILLS" = kl_mills,"LATENT CLASS"= kl_lc,"TUCKER"= kl_sf)),
	  "WASSERSTEIN" = list("EMPIRICAL" = data.frame("MILLS" = wass_mills,"LATENT CLASS"= wass_lc,"TUCKER"= wass_sf)))


#+++++++
# FOURTH
#+++++++
cat("FOURTH SCENARIO", "\n")
Y = read.table("../02_COMPETITORS/data4.txt")
for(j in 1:NCOL(Y)) Y[,j] = factor(Y[,j], levels = 1:4, labels = 1:4)
p = NCOL(Y)
d=4
H = 5
load("../01_MILLS/SIM4_biv.RData")
load("../02_COMPETITORS/SIMPLEX_FACTOR/SIM4_SF_biv.RData")
load("../02_COMPETITORS/LATENT_CLASS//SIM4_LC_biv.RData")
emp_biv = apply(all_pairs, 2, function(x) 1/NROW(Y)*table(Y[,x[1]],Y[,x[2]])) 

#+++++++++++++++++
# Person residuals
#+++++++++++++++++
pears_mills = colSums(((apply(biv_pmean, 1, c) - emp_biv)^2) / apply(biv_pmean, 1, c))
pears_lc = colSums( ( (apply(biv_pmean_LC, 1, c) - emp_biv)^2 ) / apply(biv_pmean_LC, 1, c) )
pears_sf = colSums( ( (apply(biv_pmean_SF, 1, c) - emp_biv)^2 ) / apply(biv_pmean_SF, 1, c) )


#++++++++++++
# kl
#++++++++++++
kl_sf = kl_lc = numeric(NCOL(all_pairs))
kl_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
	kl_mills[k] = emp_kl( emp_biv[,k], c(biv_pmean[k,,]))
	kl_lc[k] = emp_kl( emp_biv[,k], c(biv_pmean_LC[k,,]))
	kl_sf[k] = emp_kl( emp_biv[,k], c(biv_pmean_SF[k,,])) 
}
#++++++++++++
# Wasserstein
#++++++++++++
require(transport)

wass_sf = wass_lc = numeric(NCOL(all_pairs))
wass_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
p1 = pgrid(mass = matrix(emp_biv[,k],4,4))
p2 = pgrid(mass =biv_pmean[k,,])
wass_mills[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_pmean_LC[k,,])
wass_lc[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_pmean_SF[k,,])
wass_sf[k] = wasserstein(p1,p2)

}



dfS4 = list("PEARSON" = list("EMPIRICAL" = data.frame("MILLS" = sqrt(pears_mills/3),"LATENT CLASS"= sqrt(pears_lc/3), "TUCKER"= sqrt(pears_sf/3))),
	  "KULLBACK–LEIBLER" = list("EMPIRICAL" = data.frame("MILLS" = kl_mills,"LATENT CLASS"= kl_lc,"TUCKER"= kl_sf)),
	  "WASSERSTEIN" = list("EMPIRICAL" = data.frame("MILLS" = wass_mills,"LATENT CLASS"= wass_lc,"TUCKER"= wass_sf)))



df_comb = reshape2::melt(list("Scenario 1" = dfS1, "Scenario 2" = dfS2, "Scenario 3" = dfS3, "Scenario 4" = dfS4))

df_comb$L2 = factor(df_comb$L2, labels = c("Kullback-Liebler","Pearson","Wasserstein" ))
df_comb$variable = factor(df_comb$variable, labels = c("MILLS", "Latent Class model", "Simplex Factor model"))

df_comb$Q = "Posterior\nmean"
df_mean = df_comb


rm(list = ls()[-which(ls()=="df_mean")])
emp_kl = function(x,y) sum(x*(log(x) - log(y)),na.rm = T)
cat("+++++++++++++++++++",'\n')
cat("POSTERIOR QUANTILE",'\n')
cat("+++++++++++++++++++",'\n')

Y = read.table("../02_COMPETITORS/data1.txt")
for(j in 1:NCOL(Y)) Y[,j] = factor(Y[,j], levels = 1:4, labels = 1:4)
p = NCOL(Y)
d=4
H = 5
all_pairs = combn(p,2)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Now compare goodness of fit (kl divergence, pearson residuals, wasserstein)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
load("../01_MILLS/SIM1_biv.RData")
load("../02_COMPETITORS/SIMPLEX_FACTOR/SIM1_SF_biv.RData")
load("../02_COMPETITORS/LATENT_CLASS//SIM1_LC_biv.RData")
## Empirical
emp_biv = apply(all_pairs, 2, function(x) 1/NROW(Y)*table(Y[,x[1]],Y[,x[2]])) 

for(bb in 1:NCOL(all_pairs)){
	biv_qLO[bb,,] = biv_qLO[bb,,] / sum(biv_qLO[bb,,])
	biv_qLO_LC[bb,,] = biv_qLO_LC[bb,,] / sum(biv_qLO_LC[bb,,])
	biv_qLO_SF[bb,,] = biv_qLO_SF[bb,,] / sum(biv_qLO_SF[bb,,])
}

#+++++++++++++++++
# Person residuals
#+++++++++++++++++
pears_mills = colSums(((apply(biv_qLO, 1, c) - emp_biv)^2) / apply(biv_qLO, 1, c))
pears_lc = colSums( ( (apply(biv_qLO_LC, 1, c) - emp_biv)^2 ) / apply(biv_qLO_LC, 1, c) )
pears_sf = colSums( ( (apply(biv_qLO_SF, 1, c) - emp_biv)^2 ) / apply(biv_qLO_SF, 1, c) )


#++++++++++++
# kl
#++++++++++++
kl_sf = kl_lc = numeric(NCOL(all_pairs))
kl_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
	kl_mills[k] = emp_kl( emp_biv[,k], c(biv_qLO[k,,]))
	kl_lc[k] = emp_kl( emp_biv[,k], c(biv_qLO_LC[k,,]))
	kl_sf[k] = emp_kl( emp_biv[,k], c(biv_qLO_SF[k,,])) 
}
#++++++++++++
# Wasserstein
#++++++++++++
require(transport)

wass_sf = wass_lc = numeric(NCOL(all_pairs))
wass_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
p1 = pgrid(mass = matrix(emp_biv[,k],4,4))
p2 = pgrid(mass =biv_qLO[k,,])
wass_mills[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_qLO_LC[k,,])
wass_lc[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_qLO_SF[k,,])
wass_sf[k] = wasserstein(p1,p2)

}

dfS1 = list("PEARSON" = list("EMPIRICAL" = data.frame("MILLS" = sqrt(pears_mills/3),"LATENT CLASS"= sqrt(pears_lc/3), "TUCKER"= sqrt(pears_sf/3))),
	  "KULLBACK–LEIBLER" = list("EMPIRICAL" = data.frame("MILLS" = kl_mills,"LATENT CLASS"= kl_lc,"TUCKER"= kl_sf)),
	  "WASSERSTEIN" = list("EMPIRICAL" = data.frame("MILLS" = wass_mills,"LATENT CLASS"= wass_lc,"TUCKER"= wass_sf)))

#+++++++
# SECOND
#+++++++
cat("SECOND SCENARIO", "\n")
Y = read.table("../02_COMPETITORS/data2.txt")
for(j in 1:NCOL(Y)) Y[,j] = factor(Y[,j], levels = 1:4, labels = 1:4)
p = NCOL(Y)
d=4
H = 5
load("../01_MILLS/SIM2_biv.RData")
load("../02_COMPETITORS/SIMPLEX_FACTOR/SIM2_SF_biv.RData")
load("../02_COMPETITORS/LATENT_CLASS//SIM2_LC_biv.RData")
## Empirical
emp_biv = apply(all_pairs, 2, function(x) 1/NROW(Y)*table(Y[,x[1]],Y[,x[2]])) 
for(bb in 1:NCOL(all_pairs)){
	biv_qLO[bb,,] = biv_qLO[bb,,] / sum(biv_qLO[bb,,])
	biv_qLO_LC[bb,,] = biv_qLO_LC[bb,,] / sum(biv_qLO_LC[bb,,])
	biv_qLO_SF[bb,,] = biv_qLO_SF[bb,,] / sum(biv_qLO_SF[bb,,])
}

#+++++++++++++++++
# Person residuals
#+++++++++++++++++
pears_mills = colSums(((apply(biv_qLO, 1, c) - emp_biv)^2) / apply(biv_qLO, 1, c))
pears_lc = colSums( ( (apply(biv_qLO_LC, 1, c) - emp_biv)^2 ) / apply(biv_qLO_LC, 1, c) )
pears_sf = colSums( ( (apply(biv_qLO_SF, 1, c) - emp_biv)^2 ) / apply(biv_qLO_SF, 1, c) )

#++++++++++++
# kl
#++++++++++++
kl_sf = kl_lc = numeric(NCOL(all_pairs))
kl_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
	kl_mills[k] = emp_kl( emp_biv[,k], c(biv_qLO[k,,]))
	kl_lc[k] = emp_kl( emp_biv[,k], c(biv_qLO_LC[k,,]))
	kl_sf[k] = emp_kl( emp_biv[,k], c(biv_qLO_SF[k,,]) )
}
#++++++++++++
# Wasserstein
#++++++++++++
require(transport)

wass_sf = wass_lc = numeric(NCOL(all_pairs))
wass_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
p1 = pgrid(mass = matrix(emp_biv[,k],4,4))
p2 = pgrid(mass =biv_qLO[k,,])
wass_mills[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_qLO_LC[k,,])
wass_lc[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_qLO_SF[k,,])
wass_sf[k] = wasserstein(p1,p2)

}

dfS2 = list("PEARSON" = list("EMPIRICAL" = data.frame("MILLS" = sqrt(pears_mills/3),"LATENT CLASS"= sqrt(pears_lc/3), "TUCKER"= sqrt(pears_sf/3))),
	  "KULLBACK–LEIBLER" = list("EMPIRICAL" = data.frame("MILLS" = kl_mills,"LATENT CLASS"= kl_lc,"TUCKER"= kl_sf)),
	  "WASSERSTEIN" = list("EMPIRICAL" = data.frame("MILLS" = wass_mills,"LATENT CLASS"= wass_lc,"TUCKER"= wass_sf)))

#+++++++
# THIRD
#+++++++
cat("THIRD SCENARIO", "\n")
Y = read.table("../02_COMPETITORS/data3.txt")
for(j in 1:NCOL(Y)) Y[,j] = factor(Y[,j], levels = 1:4, labels = 1:4)
all_pairs = combn(p,2)

load("../01_MILLS/SIM3_biv.RData")
load("../02_COMPETITORS/SIMPLEX_FACTOR/SIM3_SF_biv.RData")
load("../02_COMPETITORS/LATENT_CLASS//SIM3_LC_biv.RData")
## Empirical
emp_biv = apply(all_pairs, 2, function(x) 1/NROW(Y)*table(Y[,x[1]],Y[,x[2]])) 
for(bb in 1:NCOL(all_pairs)){
	biv_qLO[bb,,] = biv_qLO[bb,,] / sum(biv_qLO[bb,,])
	biv_qLO_LC[bb,,] = biv_qLO_LC[bb,,] / sum(biv_qLO_LC[bb,,])
	biv_qLO_SF[bb,,] = biv_qLO_SF[bb,,] / sum(biv_qLO_SF[bb,,])
}

#+++++++++++++++++
# Person residuals
#+++++++++++++++++
pears_mills = colSums(((apply(biv_qLO, 1, c) - emp_biv)^2) / apply(biv_qLO, 1, c))
pears_lc = colSums( ( (apply(biv_qLO_LC, 1, c) - emp_biv)^2 ) / apply(biv_qLO_LC, 1, c) )
pears_sf = colSums( ( (apply(biv_qLO_SF, 1, c) - emp_biv)^2 ) / apply(biv_qLO_SF, 1, c) )

#++++++++++++
# kl
#++++++++++++
kl_sf = kl_lc = numeric(NCOL(all_pairs))
kl_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
	kl_mills[k] = emp_kl( emp_biv[,k], c(biv_qLO[k,,]))
	kl_lc[k] = emp_kl( emp_biv[,k], c(biv_qLO_LC[k,,]))
	kl_sf[k] = emp_kl( emp_biv[,k], c(biv_qLO_SF[k,,]) )
}
#++++++++++++
# Wasserstein
#++++++++++++
require(transport)

wass_sf = wass_lc = numeric(NCOL(all_pairs))
wass_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
p1 = pgrid(mass = matrix(emp_biv[,k],4,4))
p2 = pgrid(mass =biv_qLO[k,,])
wass_mills[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_qLO_LC[k,,])
wass_lc[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_qLO_SF[k,,])
wass_sf[k] = wasserstein(p1,p2)

}

dfS3 = list("PEARSON" = list("EMPIRICAL" = data.frame("MILLS" = sqrt(pears_mills/3),"LATENT CLASS"= sqrt(pears_lc/3), "TUCKER"= sqrt(pears_sf/3))),
	  "KULLBACK–LEIBLER" = list("EMPIRICAL" = data.frame("MILLS" = kl_mills,"LATENT CLASS"= kl_lc,"TUCKER"= kl_sf)),
	  "WASSERSTEIN" = list("EMPIRICAL" = data.frame("MILLS" = wass_mills,"LATENT CLASS"= wass_lc,"TUCKER"= wass_sf)))


#+++++++
# FOURTH
#+++++++
cat("FOURTH SCENARIO", "\n")
Y = read.table("../02_COMPETITORS/data4.txt")
for(j in 1:NCOL(Y)) Y[,j] = factor(Y[,j], levels = 1:4, labels = 1:4)
p = NCOL(Y)
d=4
H = 5
load("../01_MILLS/SIM4_biv.RData")
load("../02_COMPETITORS/SIMPLEX_FACTOR/SIM4_SF_biv.RData")
load("../02_COMPETITORS/LATENT_CLASS//SIM4_LC_biv.RData")
emp_biv = apply(all_pairs, 2, function(x) 1/NROW(Y)*table(Y[,x[1]],Y[,x[2]])) 
for(bb in 1:NCOL(all_pairs)){
	biv_qLO[bb,,] = biv_qLO[bb,,] / sum(biv_qLO[bb,,])
	biv_qLO_LC[bb,,] = biv_qLO_LC[bb,,] / sum(biv_qLO_LC[bb,,])
	biv_qLO_SF[bb,,] = biv_qLO_SF[bb,,] / sum(biv_qLO_SF[bb,,])
}

#+++++++++++++++++
# Person residuals
#+++++++++++++++++
pears_mills = colSums(((apply(biv_qLO, 1, c) - emp_biv)^2) / apply(biv_qLO, 1, c))
pears_lc = colSums( ( (apply(biv_qLO_LC, 1, c) - emp_biv)^2 ) / apply(biv_qLO_LC, 1, c) )
pears_sf = colSums( ( (apply(biv_qLO_SF, 1, c) - emp_biv)^2 ) / apply(biv_qLO_SF, 1, c) )


#++++++++++++
# kl
#++++++++++++
kl_sf = kl_lc = numeric(NCOL(all_pairs))
kl_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
	kl_mills[k] = emp_kl( emp_biv[,k], c(biv_qLO[k,,]))
	kl_lc[k] = emp_kl( emp_biv[,k], c(biv_qLO_LC[k,,]))
	kl_sf[k] = emp_kl( emp_biv[,k], c(biv_qLO_SF[k,,])) 
}
#++++++++++++
# Wasserstein
#++++++++++++
require(transport)

wass_sf = wass_lc = numeric(NCOL(all_pairs))
wass_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
p1 = pgrid(mass = matrix(emp_biv[,k],4,4))
p2 = pgrid(mass =biv_qLO[k,,])
wass_mills[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_qLO_LC[k,,])
wass_lc[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_qLO_SF[k,,])
wass_sf[k] = wasserstein(p1,p2)

}



dfS4 = list("PEARSON" = list("EMPIRICAL" = data.frame("MILLS" = sqrt(pears_mills/3),"LATENT CLASS"= sqrt(pears_lc/3), "TUCKER"= sqrt(pears_sf/3))),
	  "KULLBACK–LEIBLER" = list("EMPIRICAL" = data.frame("MILLS" = kl_mills,"LATENT CLASS"= kl_lc,"TUCKER"= kl_sf)),
	  "WASSERSTEIN" = list("EMPIRICAL" = data.frame("MILLS" = wass_mills,"LATENT CLASS"= wass_lc,"TUCKER"= wass_sf)))



df_comb = reshape2::melt(list("Scenario 1" = dfS1, "Scenario 2" = dfS2, "Scenario 3" = dfS3, "Scenario 4" = dfS4))

df_comb$L2 = factor(df_comb$L2, labels = c("Kullback-Liebler","Pearson","Wasserstein" ))
df_comb$variable = factor(df_comb$variable, labels = c("MILLS", "Latent Class model", "Simplex Factor model"))

df_comb$Q = "Posterior\n0.025 quantile"
df_LO = df_comb


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list = ls()[-which(ls() %in% c("df_mean","df_LO"))])
emp_kl = function(x,y) sum(x*(log(x) - log(y)),na.rm = T)
cat("+++++++++++++++++++",'\n')
cat("POSTERIOR QUANTILE",'\n')
cat("+++++++++++++++++++",'\n')

Y = read.table("../02_COMPETITORS/data1.txt")
for(j in 1:NCOL(Y)) Y[,j] = factor(Y[,j], levels = 1:4, labels = 1:4)
p = NCOL(Y)
d=4
H = 5
all_pairs = combn(p,2)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Now compare goodness of fit (kl divergence, pearson residuals, wasserstein)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
load("../01_MILLS/SIM1_biv.RData")
load("../02_COMPETITORS/SIMPLEX_FACTOR/SIM1_SF_biv.RData")
load("../02_COMPETITORS/LATENT_CLASS//SIM1_LC_biv.RData")
## Empirical
emp_biv = apply(all_pairs, 2, function(x) 1/NROW(Y)*table(Y[,x[1]],Y[,x[2]])) 

for(bb in 1:NCOL(all_pairs)){
	biv_qUP[bb,,] = biv_qUP[bb,,] / sum(biv_qUP[bb,,])
	biv_qUP_LC[bb,,] = biv_qUP_LC[bb,,] / sum(biv_qUP_LC[bb,,])
	biv_qUP_SF[bb,,] = biv_qUP_SF[bb,,] / sum(biv_qUP_SF[bb,,])
}

#+++++++++++++++++
# Person residuals
#+++++++++++++++++
pears_mills = colSums(((apply(biv_qUP, 1, c) - emp_biv)^2) / apply(biv_qUP, 1, c))
pears_lc = colSums( ( (apply(biv_qUP_LC, 1, c) - emp_biv)^2 ) / apply(biv_qUP_LC, 1, c) )
pears_sf = colSums( ( (apply(biv_qUP_SF, 1, c) - emp_biv)^2 ) / apply(biv_qUP_SF, 1, c) )


#++++++++++++
# kl
#++++++++++++
kl_sf = kl_lc = numeric(NCOL(all_pairs))
kl_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
	kl_mills[k] = emp_kl( emp_biv[,k], c(biv_qUP[k,,]))
	kl_lc[k] = emp_kl( emp_biv[,k], c(biv_qUP_LC[k,,]))
	kl_sf[k] = emp_kl( emp_biv[,k], c(biv_qUP_SF[k,,])) 
}
#++++++++++++
# Wasserstein
#++++++++++++
require(transport)

wass_sf = wass_lc = numeric(NCOL(all_pairs))
wass_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
p1 = pgrid(mass = matrix(emp_biv[,k],4,4))
p2 = pgrid(mass =biv_qUP[k,,])
wass_mills[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_qUP_LC[k,,])
wass_lc[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_qUP_SF[k,,])
wass_sf[k] = wasserstein(p1,p2)

}

dfS1 = list("PEARSON" = list("EMPIRICAL" = data.frame("MILLS" = sqrt(pears_mills/3),"LATENT CLASS"= sqrt(pears_lc/3), "TUCKER"= sqrt(pears_sf/3))),
	  "KULLBACK–LEIBLER" = list("EMPIRICAL" = data.frame("MILLS" = kl_mills,"LATENT CLASS"= kl_lc,"TUCKER"= kl_sf)),
	  "WASSERSTEIN" = list("EMPIRICAL" = data.frame("MILLS" = wass_mills,"LATENT CLASS"= wass_lc,"TUCKER"= wass_sf)))

#+++++++
# SECOND
#+++++++
cat("SECOND SCENARIO", "\n")
Y = read.table("../02_COMPETITORS/data2.txt")
for(j in 1:NCOL(Y)) Y[,j] = factor(Y[,j], levels = 1:4, labels = 1:4)
p = NCOL(Y)
d=4
H = 5
load("../01_MILLS/SIM2_biv.RData")
load("../02_COMPETITORS/SIMPLEX_FACTOR/SIM2_SF_biv.RData")
load("../02_COMPETITORS/LATENT_CLASS//SIM2_LC_biv.RData")
## Empirical
emp_biv = apply(all_pairs, 2, function(x) 1/NROW(Y)*table(Y[,x[1]],Y[,x[2]])) 
for(bb in 1:NCOL(all_pairs)){
	biv_qUP[bb,,] = biv_qUP[bb,,] / sum(biv_qUP[bb,,])
	biv_qUP_LC[bb,,] = biv_qUP_LC[bb,,] / sum(biv_qUP_LC[bb,,])
	biv_qUP_SF[bb,,] = biv_qUP_SF[bb,,] / sum(biv_qUP_SF[bb,,])
}

#+++++++++++++++++
# Person residuals
#+++++++++++++++++
pears_mills = colSums(((apply(biv_qUP, 1, c) - emp_biv)^2) / apply(biv_qUP, 1, c))
pears_lc = colSums( ( (apply(biv_qUP_LC, 1, c) - emp_biv)^2 ) / apply(biv_qUP_LC, 1, c) )
pears_sf = colSums( ( (apply(biv_qUP_SF, 1, c) - emp_biv)^2 ) / apply(biv_qUP_SF, 1, c) )

#++++++++++++
# kl
#++++++++++++
kl_sf = kl_lc = numeric(NCOL(all_pairs))
kl_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
	kl_mills[k] = emp_kl( emp_biv[,k], c(biv_qUP[k,,]))
	kl_lc[k] = emp_kl( emp_biv[,k], c(biv_qUP_LC[k,,]))
	kl_sf[k] = emp_kl( emp_biv[,k], c(biv_qUP_SF[k,,]) )
}
#++++++++++++
# Wasserstein
#++++++++++++
require(transport)

wass_sf = wass_lc = numeric(NCOL(all_pairs))
wass_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
p1 = pgrid(mass = matrix(emp_biv[,k],4,4))
p2 = pgrid(mass =biv_qUP[k,,])
wass_mills[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_qUP_LC[k,,])
wass_lc[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_qUP_SF[k,,])
wass_sf[k] = wasserstein(p1,p2)

}

dfS2 = list("PEARSON" = list("EMPIRICAL" = data.frame("MILLS" = sqrt(pears_mills/3),"LATENT CLASS"= sqrt(pears_lc/3), "TUCKER"= sqrt(pears_sf/3))),
	  "KULLBACK–LEIBLER" = list("EMPIRICAL" = data.frame("MILLS" = kl_mills,"LATENT CLASS"= kl_lc,"TUCKER"= kl_sf)),
	  "WASSERSTEIN" = list("EMPIRICAL" = data.frame("MILLS" = wass_mills,"LATENT CLASS"= wass_lc,"TUCKER"= wass_sf)))

#+++++++
# THIRD
#+++++++
cat("THIRD SCENARIO", "\n")
Y = read.table("../02_COMPETITORS/data3.txt")
for(j in 1:NCOL(Y)) Y[,j] = factor(Y[,j], levels = 1:4, labels = 1:4)
all_pairs = combn(p,2)

load("../01_MILLS/SIM3_biv.RData")
load("../02_COMPETITORS/SIMPLEX_FACTOR/SIM3_SF_biv.RData")
load("../02_COMPETITORS/LATENT_CLASS//SIM3_LC_biv.RData")
## Empirical
emp_biv = apply(all_pairs, 2, function(x) 1/NROW(Y)*table(Y[,x[1]],Y[,x[2]])) 
for(bb in 1:NCOL(all_pairs)){
	biv_qUP[bb,,] = biv_qUP[bb,,] / sum(biv_qUP[bb,,])
	biv_qUP_LC[bb,,] = biv_qUP_LC[bb,,] / sum(biv_qUP_LC[bb,,])
	biv_qUP_SF[bb,,] = biv_qUP_SF[bb,,] / sum(biv_qUP_SF[bb,,])
}

#+++++++++++++++++
# Person residuals
#+++++++++++++++++
pears_mills = colSums(((apply(biv_qUP, 1, c) - emp_biv)^2) / apply(biv_qUP, 1, c))
pears_lc = colSums( ( (apply(biv_qUP_LC, 1, c) - emp_biv)^2 ) / apply(biv_qUP_LC, 1, c) )
pears_sf = colSums( ( (apply(biv_qUP_SF, 1, c) - emp_biv)^2 ) / apply(biv_qUP_SF, 1, c) )

#++++++++++++
# kl
#++++++++++++
kl_sf = kl_lc = numeric(NCOL(all_pairs))
kl_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
	kl_mills[k] = emp_kl( emp_biv[,k], c(biv_qUP[k,,]))
	kl_lc[k] = emp_kl( emp_biv[,k], c(biv_qUP_LC[k,,]))
	kl_sf[k] = emp_kl( emp_biv[,k], c(biv_qUP_SF[k,,]) )
}
#++++++++++++
# Wasserstein
#++++++++++++
require(transport)

wass_sf = wass_lc = numeric(NCOL(all_pairs))
wass_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
p1 = pgrid(mass = matrix(emp_biv[,k],4,4))
p2 = pgrid(mass =biv_qUP[k,,])
wass_mills[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_qUP_LC[k,,])
wass_lc[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_qUP_SF[k,,])
wass_sf[k] = wasserstein(p1,p2)

}

dfS3 = list("PEARSON" = list("EMPIRICAL" = data.frame("MILLS" = sqrt(pears_mills/3),"LATENT CLASS"= sqrt(pears_lc/3), "TUCKER"= sqrt(pears_sf/3))),
	  "KULLBACK–LEIBLER" = list("EMPIRICAL" = data.frame("MILLS" = kl_mills,"LATENT CLASS"= kl_lc,"TUCKER"= kl_sf)),
	  "WASSERSTEIN" = list("EMPIRICAL" = data.frame("MILLS" = wass_mills,"LATENT CLASS"= wass_lc,"TUCKER"= wass_sf)))


#+++++++
# FOURTH
#+++++++
cat("FOURTH SCENARIO", "\n")
Y = read.table("../02_COMPETITORS/data4.txt")
for(j in 1:NCOL(Y)) Y[,j] = factor(Y[,j], levels = 1:4, labels = 1:4)
p = NCOL(Y)
d=4
H = 5
load("../01_MILLS/SIM4_biv.RData")
load("../02_COMPETITORS/SIMPLEX_FACTOR/SIM4_SF_biv.RData")
load("../02_COMPETITORS/LATENT_CLASS//SIM4_LC_biv.RData")
emp_biv = apply(all_pairs, 2, function(x) 1/NROW(Y)*table(Y[,x[1]],Y[,x[2]])) 
for(bb in 1:NCOL(all_pairs)){
	biv_qUP[bb,,] = biv_qUP[bb,,] / sum(biv_qUP[bb,,])
	biv_qUP_LC[bb,,] = biv_qUP_LC[bb,,] / sum(biv_qUP_LC[bb,,])
	biv_qUP_SF[bb,,] = biv_qUP_SF[bb,,] / sum(biv_qUP_SF[bb,,])
}

#+++++++++++++++++
# Person residuals
#+++++++++++++++++
pears_mills = colSums(((apply(biv_qUP, 1, c) - emp_biv)^2) / apply(biv_qUP, 1, c))
pears_lc = colSums( ( (apply(biv_qUP_LC, 1, c) - emp_biv)^2 ) / apply(biv_qUP_LC, 1, c) )
pears_sf = colSums( ( (apply(biv_qUP_SF, 1, c) - emp_biv)^2 ) / apply(biv_qUP_SF, 1, c) )


#++++++++++++
# kl
#++++++++++++
kl_sf = kl_lc = numeric(NCOL(all_pairs))
kl_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
	kl_mills[k] = emp_kl( emp_biv[,k], c(biv_qUP[k,,]))
	kl_lc[k] = emp_kl( emp_biv[,k], c(biv_qUP_LC[k,,]))
	kl_sf[k] = emp_kl( emp_biv[,k], c(biv_qUP_SF[k,,])) 
}
#++++++++++++
# Wasserstein
#++++++++++++
require(transport)

wass_sf = wass_lc = numeric(NCOL(all_pairs))
wass_mills = numeric(NCOL(all_pairs))
for (k in 1:length(kl_lc)) {
p1 = pgrid(mass = matrix(emp_biv[,k],4,4))
p2 = pgrid(mass =biv_qUP[k,,])
wass_mills[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_qUP_LC[k,,])
wass_lc[k] = wasserstein(p1,p2)

p2 = pgrid(mass = biv_qUP_SF[k,,])
wass_sf[k] = wasserstein(p1,p2)

}



dfS4 = list("PEARSON" = list("EMPIRICAL" = data.frame("MILLS" = sqrt(pears_mills/3),"LATENT CLASS"= sqrt(pears_lc/3), "TUCKER"= sqrt(pears_sf/3))),
	  "KULLBACK–LEIBLER" = list("EMPIRICAL" = data.frame("MILLS" = kl_mills,"LATENT CLASS"= kl_lc,"TUCKER"= kl_sf)),
	  "WASSERSTEIN" = list("EMPIRICAL" = data.frame("MILLS" = wass_mills,"LATENT CLASS"= wass_lc,"TUCKER"= wass_sf)))



df_comb = reshape2::melt(list("Scenario 1" = dfS1, "Scenario 2" = dfS2, "Scenario 3" = dfS3, "Scenario 4" = dfS4))

df_comb$L2 = factor(df_comb$L2, labels = c("Kullback-Liebler","Pearson","Wasserstein" ))
df_comb$variable = factor(df_comb$variable, labels = c("MILLS", "Latent Class model", "Simplex Factor model"))

df_comb$Q = "Posterior\n0.975 quantile"
df_UP = df_comb


require(tidyverse)
df_comb = bind_rows(df_mean,df_LO,df_UP)
ll = levels(factor(df_comb$Q))
df_comb$Q = ordered(df_comb$Q,levels=c(ll[3],ll[1],ll[2]))
df_comb$N = "Posterior"
pl_tot = ggplot(df_comb) + 
	geom_boxplot(aes(y=value,x=L2,fill=variable),alpha = .75,show.legend=T) + 
	scale_fill_manual(values = c("#FABD2F","#FB4934", "#83a598"),name="Method")+
	#geom_jitter(aes(y=value,x=variable,fill=L2),alpha = .1,width = .2) + 
	xlab("")+ylab("")+
	ylim(0,0.2)+
	#geom_point(aes(L2,value,fill=variable,col=variable), size=-1)+
	scale_x_discrete(labels=ggplot2:::parse_safe)+
	#scale_y_continuous(breaks = seq(0,0.35,by=0.05))+
	facet_grid(rows=vars(Q), col=vars(L1),scales = "fixed") + 
	theme_light(base_size = 25)+
	coord_flip()+
	theme(legend.title = element_text(face = "bold"),
	      strip.text.x = element_text(size = 25, face = "bold"),
	      strip.text.y = element_text(size = 25, face = "bold"),
	      axis.text.x = element_text(angle=90),
	      legend.position = "bottom",
	      legend.text = element_text(size=25),
	      plot.margin=unit(c(0.3,0.1,0.3,-0.2), "cm"),
	      legend.box = "horizontal",
	      legend.background = element_rect(fill="lightgray", size=0.5, linetype="solid")
	) 
pl_tot
	

#ggsave(pl_tot, file = "simulation_complete.pdf",width = 20,height = 13)




