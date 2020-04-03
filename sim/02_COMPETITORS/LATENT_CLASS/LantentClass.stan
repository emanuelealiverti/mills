data {
	int<lower=1> N; // num observations
	int<lower=1> P; // Variables
	int<lower=1> D; // Categories
	int<lower=2> H; // Components
	vector<lower=0>[H] alpha;
	vector<lower=0>[D] alpha_kern;
	int<lower=1> X[N,P];
}

parameters{
	simplex[H] nu;
	simplex[D] prob[P,H];
}

model {
	// vector[P] temp;
	real temp;
	real ps[H];
	real contr;
	matrix[D,H] prob_resh;
	matrix[D,H] prob_z;
	vector[P-1] mi_vec;
	// vector[D] prob_z;

	nu ~ dirichlet(alpha);

	// Marginal specific probability

	for(h in 1:H){
		for(p in 1:P){
			prob[p,h] ~ dirichlet(alpha_kern);
		}
	}

	// loop over observartion
	for(i in 1:N){
		for(h in 1:H)
		{
			temp  = 0.0;
			for(p in 1:P)
			{temp  +=  log(prob[p,h][X[i,p]]);
			}

			ps[h] = log(nu[h]) + temp;

		}
		target += log_sum_exp(ps);
	}
}
