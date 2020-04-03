#+++++++++++++++++++++++++
# Gibbs 
#+++++++++++++++++++++++++
gibbsPPssl = function(Y, d, 
		   H = 2, 
		   Nit = 1e3, burnin = 500, 
		   prior = list(s0 = 1), 
		   verbose = print_gibbs(), ...) {


	#++++++++++++++++++++++++++++++++
	# Extract quantities from data 
	#++++++++++++++++++++++++++++++++

	N       =  NROW(Y)
	P       =  NCOL(Y)
	P_p     =  t(combn(P, 2))
	n_pairs =  NROW(P_p)
	Y_n = apply(Y,2,as.numeric) 
	X             =  model.matrix(Freq~(.)^2, as.data.frame(table(Y[,1:2])))
	X_inv         =  solve(X)

	#+++++++++++++++++++++
	# PRIOR specification
	#+++++++++++++++++++++

	# the model (likelihood) is parametrized as a function of the probabilities. 
	# s_0 referes to the canonical parameters.
	# Therefore, we are imposing a slightly different prior on the probs (simply more inflated)
	s0         =  prior$s0 *  diag(solve(crossprod(X_inv)))
	nu0        =  rep(1/H,H)
	prior_biv  =  rep(1/d^2,d^2)
	prior_marg =  rep(1/d,d)
	#l0         =  prior$l0
	a0         =  prior$a0
	a1         =  prior$a1
	#gamma0     =  prior$gamma0

	#++++++++++++++++++++++++++++++++++++
	# Quantities used inside the sampler
	#Trasnsform into N*Ptot id. Every row which interactions contributes
	#Y_df = data.frame(Y)
	#for(j in 1:NCOL(Y_df)) Y_df[,j] = factor(Y_df[,j], levels = 1:d)
	#Y_mat_id = model.matrix(~(.)^2, data = Y_df[,1:2])

	#Every observation contributes to all the pairs. Which cell?5
	Y_wp = matrix(0, N, n_pairs)
	for(i in 1:N){
		for(j in 1:n_pairs){
			id = P_p[j, ]
			obs = as.numeric(Y[i, c(id)])
			Y_wp[i,j] = id_vec(R = obs[1], C = obs[2], d = d)
		}
	}


	Npar = d^2

	P_group = matrix(NA, N, H)
	Group_stat = matrix(0, nrow = Npar, ncol = H)


	#+++++++++++++++++++++++++
	# Starting points 
	#+++++++++++++++++++++++++
	theta      =  array(rep(1/d^2), dim = c(H, n_pairs, Npar))
	theta[,,1] =  0
	weight     =  matrix(runif(n_pairs*H), n_pairs,H)
	Z          =  numeric(N)
	nu         =  nu0
	pr_slab   =  rep(.1,H)
	prs        =  matrix(.5,n_pairs,H)
	id_w       =  matrix(rbinom(n_pairs*H, 1, prob = .5 ),n_pairs,H)

	#+++++++++++++++++++++++++
	# Storage 
	#+++++++++++++++++++++++++
	theta_out   =  array(0, dim = c(Nit, H, n_pairs, Npar))
	nu_out      =  array(0, dim = c(Nit, H))
	N_clust_out =  array(0, dim = c(Nit, H))
	weight_out  =  array(0, dim = c(Nit, n_pairs, H))
	prs_out     =  array(0, dim = c(Nit, n_pairs, H))

	#+++++++++++++++++++++++++
	# Strings 
	#+++++++++++++++++++++++++
	time_str = paste("Sampling starts now.", format(Sys.time(), "%H:%M:%S"))

	print_plus_row(nchar(time_str))
	cat(time_str,'\n')
	print_plus_row(nchar(time_str))
	start_at = Sys.time()

	for(it in 1:(Nit + burnin)) {

		## Compute likelihood contribution in each cluster
		for(h in 1:H){
			loglike = dlmnomPPw(Y_wp = Y_wp, theta_c = tcrossprod(theta[h,,],X_inv),w = weight[,h])
			#loglike = dlmnomPP(Y_wp = Y_wp, theta_c = tcrossprod(theta[h,,],X_inv))
			P_group[,h] = log(nu[h]) + loglike
		}

		#Cluster assignment
		# 1] transform into prob. Avoid numerical instability
		pp = t( apply(P_group, 1, log_exp_norm) )

		#2] initial iteratrions leads to very low probabilities. We sample from the prior
		pp[rowSums(pp) == 0, ] = rep(1/H,H)

		for(i in 1:N){
			Z[i] = sample(1:H, 1, prob = pp[i,])
		}

		N_clust = tabulate(Z, H)

		## Update cluster specific sufficient statistics
		kappa_s    =  array(0, dim = c(H, n_pairs, d^2))
		Group_stat =  array(0, dim = c(H, n_pairs, d^2))
		for(h in 1:H){
			if(N_clust[h] > 1)
			{
				Group_stat[h,,] =  suff_stat( Y[Z == h, ], d)
				kappa_s [h,,]   =  Group_stat[h,,] - N_clust[h]/2
			}
		}

		#Compute PG sufficient stat
		#kappa = Group_stat - N/2

		# Small note: the prior on the log-linear log-odds induces a different prior on  
		for(h in 1:H){
			for(k in 1:n_pairs){
				#The first theta is fixed to 0 for identif
				for(j in 2:Npar){
					c_j =  log( sum( exp( theta[h, k, -j] ) ))
					psi =  theta[h, k, j] - c_j
					#w   =  BayesLogit::rpg.devroye(num = 1, n = N_clust[h], z = psi)
					w   =  BayesLogit::rpg.devroye(h = N_clust[h], z = psi)
					s_j =  1/(w + 1/s0[j] )
					m_j =  kappa_s[h, k, j ] + w*c_j

					theta[h, k , j] = rnorm(1, mean = s_j*m_j, sd = sqrt(s_j))
				}
			}
		}


		#sample cluster weights
		nu = rdirichlet(n=1,nu0 + N_clust)

		#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		# Update of composite likelihood weights from spike and slab exponential
		#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		for(h in 1:H) {
			#Compute pairs contribution
			#log_like_PP = rowSums(theta[h,,] * suff_stat(Y,d)/N) - log(rowSums(exp(tcrossprod(theta[h,,],X_inv))))
			log_like_PP = rowSums(theta[h,,] * suff_stat(Y[Z == h,],d)) - N*log(rowSums(exp(tcrossprod(theta[h,,],X_inv))))
			#log_like_PP = log_like_PP/sum(log_like_PP)

			c_spike = dexp(x = weight[,h], a1 - log_like_PP)
			c_slab = dgamma(x = weight[,h],shape = (1 + a0), rate = a1 - log_like_PP)

			prs[,h] = ( c_slab * pr_slab[h] ) / ( c_slab * pr_slab[h] + c_spike * (1-pr_slab[h]))
		# Avoids numerical underflows
			prs[!is.finite(prs[,h]),h] = 0
			id_w[,h] = rbinom(n_pairs, 1, prob = (prs[,h]) )

			ss_pars = a1 - log_like_PP
			shapes = 1 + id_w[,h] * a0

			# sample both from a gamma relyin on the fact that a gamma with shape = 1 is an exponential
			weight[,h] = rgamma(n = length(ss_pars),shape = shapes,rate = ss_pars)

			# Update (global) spike proabilities
			pr_slab[h] = rbeta(n = 1, shape1 = 2 + sum(id_w[,h]), 
					   shape2 = 20 + n_pairs - sum(id_w[,h]))
		}



		if(it < burnin & it %% verbose$every == 0) {
			cat(sprintf("BURNIN:  It %g over %g ", it, burnin), '\n')
		}
		#Estimate ETA
		if(it == 99){
			emp_eta(start_at,Nit,burnin,it)
		}


		if(it > burnin){
			if(it == burnin) cat("SAMPLING STARTED.")
			if(it %% verbose$every == 0) cat('\t', sprintf("It %g over %g ", it-burnin, Nit), '\n')

			theta_out[it - burnin,,,] =  theta
			nu_out[it-burnin,]        =  nu
			N_clust_out[it-burnin,]   =  N_clust
			weight_out[it-burnin,,]   =  weight
			prs_out[it-burnin,,]       =  prs
		}
	}

	out         =  list()
	out$'theta' =  theta_out
	out$'nu'    =  nu_out
	out$'N_c'   =  N_clust_out
	out$'ww'    =  weight_out
	out$'prs'   =  prs_out

	cat("Sampling done", format(Sys.time(), "%H:%M:%S"), '\n')
	return(out)
}

