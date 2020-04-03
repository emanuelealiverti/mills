rm(list=ls())
library(mills, lib.loc = "../../lib/")


for(ss in as.character(1:4)){
	load(file=sprintf("Scen%s.RData",ss))

	# invert binomial coefficient
	p    =  which(choose(1:100,2) == 105)
	d    =  sqrt(dim(res$theta)[4])
	H    =  dim(res$nu)[2]
	N_it =  max(dim(res$nu))

	biv = array(NA, dim = c(N_it,dim(res$theta)[3], rep(d,2)))
	dim(biv)
	all_pairs = combn(p,2)
	tmp = matrix(NA,H,d^2)
	for(it in 1:N_it) {
		for(pp in 1:NCOL(all_pairs)){
			for(h in 1:H){
				tmp[h,] = 		c(lo2pr( ( res$theta[it,h,pp,] - log(sum(exp(res$theta[it,h,pp,]))))))

			}
			biv[it,pp,,] = matrix(colSums(tmp*res$nu[it,]),d)
		}
		if(it %% 25 == 0) cat(it,"\n")
	}

	biv_pmean =  apply(biv,c(2,3,4),mean)
	biv_qLO   =  apply(biv,c(2,3,4),quantile,.025)
	biv_qUP   =  apply(biv,c(2,3,4),quantile,.975)

	save(list = c(	"biv_pmean" , "biv_qLO"   , "biv_qUP"    ), file=sprintf("SIM%s_biv.RData",ss))
}


