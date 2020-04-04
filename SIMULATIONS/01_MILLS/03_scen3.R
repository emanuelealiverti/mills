#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Generate 5 variable from a joint pdf and the others independently
# (see Russo, Durante Scarpa (2018, CSDA)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list = ls())
library(mills, lib.loc = "../../lib/")
set.seed(123)
n = 300
p = 15
group_id = c(6:15)
group_dep = c(1,2,3,4,5)
kern = matrix(0.25,15,4)
for (j in 1:length(group_id)){
	kern[group_id[j],] = rdirichlet(1,c(10,10,10,10))
}
pmf_joint = array(0.6/(4^5-4),c(rep(4,5)))
pmf_joint[1,1,1,1,1] = 0.1
pmf_joint[2,2,2,2,2] = 0.1
pmf_joint[3,3,3,3,3] = 0.1
pmf_joint[4,4,4,4,4] = 0.1

#Vectorized probability table
(vec_pi_Y_0_joint = as.matrix(reshape2::melt(pmf_joint)))
tensor_data = matrix(0,n,p)

for (i in 1:n){
	for (j in 1:(length(group_id))){
		tensor_data[i,group_id[j]] =
			sample(c(1:4),1,replace=TRUE,prob=kern[group_id[j],])
	}

	tensor_data[i,group_dep] = c(vec_pi_Y_0_joint[sample(c(1:dim(vec_pi_Y_0_joint)[1]),1,replace=TRUE,prob=vec_pi_Y_0_joint[,6]),1:5])
}

Y = data.frame(tensor_data)
write.table(Y,file = "../02_COMPETITORS/data3.txt",row.names=F,col.names=F)
for(j in 1:NCOL(Y)) Y[,j] = factor(Y[,j], levels = 1:4, labels = 1:4)
#++++++
# GIBBS
#++++++
set.seed(123)
prior = list(s0=1, a0=100,a1=10)
res = gibbsPPssl(Y = Y, d = 4, H = 5,prior = prior)
save(res, file = 'Scen3.RData')
