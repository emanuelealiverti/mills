## Generic function to compute bivariates under the simplex factor model.
## See Bhattacharya & Dunson (2012, JASA)
## we implement a finite approximation with a Dirichlet(1/H, ..., 1/H) distribution instead of a DP. Therefore, \alpha = 1, \nu_h = 1/H for h = 1, \dots, H
biv_dist = function(p1, p2, mw){
	# p1: matrix H x D (variable 1)
	# p2: matrix H x D (varaible 2)
	D=ncol(p1)
	out = matrix(0,D,D)
	for(d1 in 1:D){
		for(d2 in 1:D){
			out[d1,d2] = .5 * (mw %*% p1[,d1] )*(mw %*% p2[,d2])+
				.5 * mw %*% (p1[,d1] * p2[,d2] ) 

		}
	}
	return(out)
}	

## For each scenario, compute posterior mean and quantiles.
for(ss in as.character(1:4)){
	load(file=sprintf("SIM%s_SF.RData",ss))

	nit = dim(tmp$eta)[1]
	D = dim(tmp$lambda)[4]
	P = dim(tmp$lambda)[2]
	np = choose(P,2)
	H = dim(tmp$lambda)[3]
	out = array(NA, dim = c(nit,choose(P,2),D,D))

	for (it in 1:nit) {
		bb = 1
		for (p1 in 1:(P-1)) {
			m1 = tmp$lambda[it, p1, , ]
			for (p2 in (p1+1):P)
			{
				m2 = tmp$lambda[it, p2, , ]
				out[it,bb,,] = biv_dist(m1,m2,rep(1/H,H))
				bb = bb+1
			}
		}
		if(it %% 100 == 0) cat(it,'\n')
	}

	biv_pmean_SF =  apply(out, c(2,3,4),mean)
	biv_qLO_SF   =  apply(out, c(2,3,4),quantile,.025)
	biv_qUP_SF   =  apply(out, c(2,3,4),quantile,.975)

	save(list = c(	"biv_pmean_SF" , "biv_qLO_SF"   , "biv_qUP_SF"   ), file=sprintf("SIM%s_SF_biv.RData",ss))
}

