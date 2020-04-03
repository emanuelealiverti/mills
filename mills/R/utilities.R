#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Utilities - CORNER PARAMETR
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# transform probability tensor into log-linear model

#transoform log-linear model into prob-tensor. 
# Create a tensor as argument since it is SUPER useful to create the formula from it
re_par_LogLinear_to_Tensor = function(LL_coef, d, p, sparse = F){
	###one parameter is redundant, check if is there or not
	x = array(0, dim = rep(d, p))
	ll = lapply(1:p, function(x) as.character(1:d))
	# I like to have xyz and than others
	dimnames(x) = ll
	p_tab = as.table(x)
	df = as.data.frame(round(p_tab))
	pairwise = as.formula('Freq~(.)^2')
	matr = model.matrix(pairwise, data = df)

	out = array((exp(matr %*% LL_coef)/sum(exp(matr %*% LL_coef))), dim = rep(d,p))
	return(out)

}

re_par_LogLinear_FULL = function(LL_coef, d, p, sparse = F){
	###one parameter is redundant, check if is there or not
	x = array(0, dim = rep(d, p))
	ll = lapply(1:p, function(x) as.character(1:d))
	# I like to have xyz and than others
	dimnames(x) = ll
	p_tab = as.table(x)
	df = as.data.frame(round(p_tab))
	pairwise = as.formula(paste0('Freq~(.)^',p))
	matr = model.matrix(pairwise, data = df)
	# Sometimes we want to specify a lot of zero coefficents
	if(length(LL_coef) < NCOL(matr)) LL_coef = c(LL_coef, rep(0, NCOL(matr) - length(LL_coef)))

	out = array((exp(matr %*% LL_coef)/sum(exp(matr %*% LL_coef))), dim = rep(d,p))
	return(out)

}
#transoform prob-tensor into LL (inverse of above). 
re_par_Tensor_to_LogLinear = function(x, int = F){
	if(!is.array(x)) stop('X should be in array form')
	d = dim(x)[1]
	p = length(dim(x))
	ll = lapply(1:p, function(x) as.character(1:d))
	# I like to have xyz and than others
	lab = c(c('x','y','z'), rev(letters)[-c(1:3)])
	names(ll) = lab[1:p]
	dimnames(x) = ll
	p_tab = as.table(x)
	df = as.data.frame(p_tab)
	saturated = as.formula(paste0('Freq~(.)^',p))
	matr = model.matrix(saturated, data = df)
	lab = colnames(matr)

	out = (solve(matr) %*% log(c(x)/c(x)[1]))[-1,1]
	if(int) out = (solve(matr) %*% log(c(x)/c(x)[1]))[,1]
	return(out)

}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Number of parameter 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
how_many_pars = function(P,d,H) H * (1 + P * (d-1) + P * (P-1) / 2 * (d - 1)^2)

#++++++++++++++++++++++++++++++++++
# Conditional gaussian mean and var
#++++++++++++++++++++++++++++++++++
cond_gauss = function(theta, mu, sigma_i, j){
	m = mu[j] - 1/sigma_i[j,j] * sigma_i[j, -j] %*% (theta[-j] - mu[-j] )
	return(list("m" = m, "p" = sigma_i[j,j]))
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Plots
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
trace_pl = function(dat, red = T)
{
	pl=NULL
	require(tidyverse)
	df = reshape2::melt(dat)
	if(red) rID = round(seq(1,max(dim(dat)),l=1000))
	if(NCOL(df) == 3){
		df = reshape2::melt(dat[rID,])
		pl = ggplot(df) + geom_line(aes(y = value, x = Var1 )) +
			facet_wrap(~Var2) + theme_bw()
	} 
	if(NCOL(df) == 4){
		df = reshape2::melt(dat[rID,,])
		P = seq(1,20)[choose(1:20,2) == max(df$Var2)]
		lab= apply(combn(P,2), 2, function(x) paste0('(', paste(x, collapse = ','),')' ))
		df$Var2 = factor(df$Var2, labels = lab)
		df = df %>% group_by(Var2,Var3) %>% mutate(cmean = cumsum(value)/(1:max(Var1)))
		df = df %>% mutate(cvar = cumsum(value^2)/(1:max(Var1)) - cmean^2)
		df = df %>% mutate(up = cmean + 1.96*sqrt(cvar), low = cmean - 1.96*(cvar)) 
		#df = df %>% group_mutate(up = quantile(value, .975), low = quantile(value, .025)) 
		pl=ggplot(df) + geom_line(aes(y = value, x = Var1 )) +
			geom_line(aes(y = cmean, x = Var1), col = 'red', alpha = .8)+
			geom_ribbon(aes(ymin = low, ymax = up, x = Var1), alpha = .1)+
			facet_wrap(~Var2+Var3) + theme_bw()
	}
	pl
}
matr_pl = function(dat){
	pl=NULL
	require(ggplot2)
	df = reshape2::melt(dat)
	P = seq(1,20)[choose(1:20,2) == max(df$Var1)]
	lab= apply(combn(P,2), 2, function(x) paste0('(', paste(x, collapse = ','),')' ))
	df$Var1 = factor(df$Var1, levels = 1:max(df$Var1), labels = lab)
	df$Var2 = factor(df$Var2, levels = 1:max(df$Var2))
	pl=ggplot(df) + geom_tile(aes(fill = value, x = Var1, y = Var2 )) +
		facet_wrap(~L1, scales = 'free') +
		scale_x_discrete(expand = c(0,0))+
		scale_y_discrete(expand = c(0,0))+
		theme_bw()
	pl
}

pr2lo = function(probs){
	log(probs/probs[1])
}
lo2pr = function(log_odds){
	tmp = exp(log_odds)
	tmp/sum(tmp)
}

triv_psi = function(theta, theta_m, d, p, ppid = NULL){

	triv    =  combn(p, 3)
	biv = combn(p,2)
	Ntri    =  NCOL(triv)
	tri_psi =  array(0, dim = c(dim(theta)[1], Ntri, d^3))

	if(is.null(ppid)) {
		ppid = array(NA, dim = c(d^3,3,3))
		x=1
		for(i in 1:d) {
			for(j in 1:d) {
				for(k in 1:d) {
					ppid[x,,] = t(rbind(1:3, combn(c(i,j,k), 2)))
					x = x+1

		}}}
	}

	#+++++++++++++++++++++++++++++
	# Loop over kernels 
	#+++++++++++++++++++++++++++++

	for(h in 1:dim(theta)[1]) {
		theta_c = theta[h,,]
		theta_mc = theta_m[h,,]

		for(tr in 1:Ntri) {

			#++++++++++++++++++++++++++++++++++++++++++++++++
			# Index of bivaraites involved in the 3 variate
			#++++++++++++++++++++++++++++++++++++++++++++++++
			id_biv = sapply(1:3, function(x) 
					which(colSums(combn(triv[, tr], 2)[,x] == biv) == 2))

			#++++++++++++++++++++++++++++++++++++++++
			# Extract the 3 bivariates of interest
			#++++++++++++++++++++++++++++++++++++++++

			temp = array(theta_c[id_biv, 1:(d^2)], dim = c(3,d,d))
			tri_arr = array(0, dim = rep(d,3))
			x=1
			for(i in 1:d) {
				for(j in 1:d) {
					for(k in 1:d) {
						cell_pairs = ppid[x,,]
						tri_arr[i,j,k] = sum(temp[cell_pairs]) - p/2*sum(theta_mc[cbind(triv[,tr],c(i,j,k))])
						x=x+1
					}
				}
			}

			tri_psi[h, tr, ] = c(tri_arr)

		}

	}
	return(tri_psi)
}

cramer_triv = function(tab)
{
	#+++++++++++++++++++++++++++
	# Compure theoretichal freq
	#+++++++++++++++++++++++++++
	marginals = sapply(1:3,function(l)  margin.table(tab,l))
	tmp = outer(outer(marginals[,1],marginals[,2], "*"), marginals[,3], "*")
	chis = sum ( (tab - tmp)^2/tmp )
	sqrt(chis / min(dim(tab)-1))
}

#+++++++++++++++++++++++++
# Generic functions to do stuff 
#+++++++++++++++++++++++++
log_exp_norm = function(x, is.log = T) {
	log_x = x
	if(!is.log) log_x = log(x)
	idx = !is.infinite(log_x)
	tmp = rep(-Inf, length(x))
	tmp[idx] = log_x[idx] - mean(log_x[idx])
	out = numeric(length(x))
	for(i in 1:length(x)) {
		out[i] = exp( -log(1 + sum(exp (tmp[-i] -tmp[i])) ))
	}
	out[is.na(out)] = 0
	out
}

emp_eta = function(start_at, Nit, burnin,it,int = 100){
	time_100 = as.double(difftime(Sys.time(), start_at, units='sec'))
	str_l = paste(sprintf("Timing (based on %g it):",int), round(time_100,4), "secs")
	print_plus_row(nchar(str_l))
	cat(str_l, '\n')
	tu = ( (Nit+burnin - it) / int) * time_100
	cat("Remaining:",tu , "secs", '\n')
	if( tu > 100) cat("Remaining:", tu/60, "min", '\n')
	print_plus_row(nchar(str_l))
}

print_plus_row = function(l){
	cat(paste(rep('+',l), collapse=''),'\n')
}

print_gibbs = function(every = 25) list(every = every)
#biv counts
suff_stat = function( Y_all, d){
	P    =  ncol(Y_all)
	pair  =  t(combn(P,2))
	Npair =  NROW(pair)
	biv   =  matrix(NA, nrow = Npair, ncol = d^2)

	for(i in 1:Npair){
		biv[i,] = c(table(Y_all[,pair[i,1] ], Y_all[,pair[i,2]]))
	}
	return(biv)
}



id_vec = function(R, C, d) (C-1)*d + R

rdirichlet = function(n=1, alpha){
    l = length(alpha)
    x = matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm = x %*% rep(1, l)
    return(x/as.vector(sm))
}

free_prob_to_constrLL = function(coef, X_inv, free_par_test, d=d){
	#Transform fixed probababilities into full
	nFix   =  length(free_par_test) + 1
	constr =  which(colSums(abs(X_inv[1:nFix,])) == 0)
	tmp    =  -X_inv[( nFix + 1):nrow(X_inv), free_par_test]
	tmpC   =  tmp %*% coef

	vec_constr        =  numeric(d^2)
	vec_constr[free_par_test] =  coef
	vec_constr[constr]  =  tmpC
	vec_constr
}

nrap_pois = function(y,X,w,theta_k,lambda = 1 / (2*s0), maxit=10,tol=1e-5) {
	it=1
	dd = 99
	while(it < maxit & dd > tol){
		eta = exp(X %*% theta_k)
		score = crossprod( w*X, (y-eta)) -c(lambda*theta_k)/2
		hess = t(X) %*% diag(c(w * eta), NROW(X)) %*% X + diag(lambda/2, NCOL(X))
		theta = theta_k + solve(hess) %*% score
		it = it+1
		dd = abs(sum(theta-theta_k))
		theta_k = theta
		#cat(it,"\n")
		#cat(theta,"\n")
	}
	theta
}


# Truncated exponential
rtexp = function(n = 1, lambda = 1) {
	tmp = runif(n)
	out = tmp * (exp(lambda)-1) + 1 
	return (log(out)/lambda)
}
	
dtexp = function(x, lambda = 1) {
	ll = lambda
	if(length(ll) != length(x)) { ll = rep(lambda, length(x)) }
	out = numeric(length(x))
	ks = (x >= 0 & x <= 1)
	out = ( lambda * exp(lambda*x) )  / (-1 + exp(lambda))
	out[!ks] = 0
	return(out)
}
	


