#++++++++++++++++++++++++++++++++++++++++
#  Generate data from a latent class model
#++++++++++++++++++++++++++++++++++++++++
rm(list = ls())
library(mills, lib.loc = "../../lib/")

set.seed(123)
n = 300
d = 4
p = 15
H = 4
nu = rdirichlet(1,rep(1/H,H))
psi = array(NA, dim = c(H,p,d))
for(h in 1:H) {for(j in 1:p) psi[h,j,] = rdirichlet(1,rep(3,d))}
z = sample(1:H,size=n,prob = nu,rep=T)
Y = matrix(NA,n,p)
for(i in 1:n) for(j in 1:p) Y[i,j] = sample(1:d,1,prob = psi[z[i],j,])
write.table(Y,file = "../02_COMPETITORS/data1.txt",row.names=F,col.names=F)
Y = data.frame(Y)
for(j in 1:NCOL(Y)) Y[,j] = factor(Y[,j], levels = 1:4, labels = 1:4)

#++++++
# GIBBS
#++++++
set.seed(123)

prior = list(s0=1, a0=100,a1=10)
res = gibbsPPssl(Y = Y, d = 4, H = 5,prior = prior)
save(res, file = 'Scen1.RData')
