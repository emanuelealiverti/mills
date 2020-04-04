#+++++++++++++++++++++++++++++++++++++++++++++
# Generate from two blocks of loglinear models
#+++++++++++++++++++++++++++++++++++++++++++++
library(mills, lib.loc = "../../lib/")
rm(list = ls())

#+++++++++++++++++++++++++
# More challenging setting
#+++++++++++++++++++++++++
set.seed(3)
n = 400
p = 5
d = 4
H = 1 
betas = rnorm(1+(d-1)*p + (d-1)*(d-1) *choose(p,2),sd = .1)
tab1 = re_par_LogLinear_FULL(betas,d = d,p = p)

prob_vec  =  data.frame(as.table(tab1))
id_seq    =  c(1:NROW(prob_vec))
id        =  sample(x = id_seq, size = n, replace  = T, prob = prob_vec$Freq)
Y         =  prob_vec[id,1:p]
p2 =  10 

betas = rnorm(1+(d-1)*p + (d-1)*(d-1) *choose(p,2) + (d-1)^3*choose(p,3),sd=.1)
tab2 = re_par_LogLinear_FULL(betas,d = d,p = p)

prob_vec  =  data.frame(as.table(tab2))
id_seq    =  c(1:NROW(prob_vec))
id        =  sample(x = id_seq, size = n, replace  = T, prob = prob_vec$Freq)
Y2         =  prob_vec[id,1:p]

p2 =  5 

pr_ind = matrix(NA, p2,d)
Y3 = matrix(NA,n,p2)
for(j in 1:p2) {
pr_ind[j,] = rdirichlet(alpha = rep(3,d))
	Y3[,j] = sample(LETTERS[1:d],n, rep=T, prob = pr_ind[j,])
}
Y = data.frame(cbind(Y,Y2,Y3))
row.names(Y) = NULL
for(j in 1:NCOL(Y)) Y[,j] = factor(Y[,j], levels = LETTERS[1:4], labels = 1:4)
write.table(Y,file = "../02_COMPETITORS/data4.txt",col.names = F,row.names = F)
#+++++++++++++++++++++
# Posterior simulation
#+++++++++++++++++++++
H = 5
prior = list(s0=1, a0=100,a1=10)
res = gibbsPPssl(Y = Y, d = 4, H = H,prior = prior)
save(res,file = 'Scen4.RData')
