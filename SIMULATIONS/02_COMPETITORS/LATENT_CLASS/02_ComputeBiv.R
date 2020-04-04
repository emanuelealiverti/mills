#+++++++++++++++++++++++++++++++++++++++++++++++++++
# Compare now with LC bivariate 
#+++++++++++++++++++++++++++++++++++++++++++++++++++
for(ss in as.character(1:4)){
	load(file=sprintf("SIM%s_LC.RData",ss))

	N_it = max(dim(tmp$nu))
	H = min(dim(tmp$nu))
	p = dim(tmp$prob)[2]
	d = dim(tmp$prob)[4]

	bivLC = array(NA, dim = c(N_it,choose(p,2), rep(d,2)))
	all_pairs = combn(p,2)
	ker = matrix(NA,H,d^2)
	for(it in 1:N_it) {
		for(pp in 1:NCOL(all_pairs)){
			for(h in 1:H){
				ker[h,] = c(tcrossprod(tmp$prob[it,all_pairs[1,pp],h,], tmp$prob[it,all_pairs[2,pp],h,]))
			}
			bivLC[it,pp,,] = matrix(colSums(ker*tmp$nu[it,]),d)
		}
		if(it %% 25 == 0) cat(it,"\n")
	}


	biv_pmean_LC =  apply(bivLC,c(2,3,4),mean)
	biv_qLO_LC   =  apply(bivLC,c(2,3,4),quantile, 0.025)
	biv_qUP_LC    =  apply(bivLC,c(2,3,4),quantile, .975)

	save(list = c(	"biv_pmean_LC" , "biv_qLO_LC"   , "biv_qUP_LC"    ), file=sprintf("SIM%s_LC_biv.RData",ss))
}

