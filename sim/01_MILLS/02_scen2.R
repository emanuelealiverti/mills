#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Generate data from a simple LL model over a small subset of variables
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(mills, lib.loc = "../../lib/")
rm(list = ls())

#+++++++++++++++++++++++++
# More challenging setting
#+++++++++++++++++++++++++
set.seed(31)
n = 400
p = 5
d = 4
H = 1 
betas = rnorm(1+(d-1)*p + (d-1)*(d-1) *choose(p,2),sd=.1)
tab = re_par_LogLinear_FULL(betas,d = d,p = p)

prob_vec  =  data.frame(as.table(tab))
id_seq    =  c(1:NROW(prob_vec))
id        =  sample(x = id_seq, size = n, replace  = T, prob = prob_vec$Freq)
Y         =  prob_vec[id,1:p]
p2 =  10 
Y2 = matrix(NA,n,p2)
pr_ind = matrix(NA,d,p2)
for(j in 1:p2) {
	pr_ind[,j] = rdirichlet(alpha = rep(3,d))
	Y2[,j] = sample(LETTERS[1:d],n, rep=T, prob = pr_ind[,j])
}
Y = data.frame(cbind(Y,Y2))
row.names(Y) = NULL
for(j in 1:NCOL(Y)) Y[,j] = factor(Y[,j], levels = LETTERS[1:4], labels = 1:4)
write.table(Y,file = "../02_COMPETITORS/data2.txt",row.names = F,col.names = F)

#+++++++++++++++++++++
# Posterior simulation
#+++++++++++++++++++++
H = 5
prior = list(s0=1, a0=100,a1=10)
res = gibbsPPssl(Y = Y, d = 4, H = 5,prior = prior)

save(res,file = 'Scen2.RData')
