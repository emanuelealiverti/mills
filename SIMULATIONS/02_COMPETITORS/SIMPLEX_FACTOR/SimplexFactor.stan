
data {
	int<lower=1> N; // num observations
	int<lower=1> P; // Variables
	int<lower=1> D; // Categories

	int<lower=2> H; // Components

	//vector<lower=0> [H] eta_prior;
	//vector<lower=0> [D] a_lambda;
	//real pen;

	int<lower=1> X[N,P];
}

parameters{
	simplex[D] lambda[P,H];
	simplex[H] eta[N];
}

model {
	vector[H] temp;

	for(i in 1:N) {
		eta[i] ~ dirichlet(rep_vector(1,H));
	}

	for(h in 1:H){
		for(p in 1:P){
			lambda[p,h] ~ dirichlet(rep_vector(1,D));
		}
	}

	// loop over observartion
	for(i in 1:N){
		for(p in 1:P){
			for(h in 1:H)
			{
				temp[h]  =  log(lambda[p,h][X[i,p]]) + log(eta[i,h]);
			}
			target += log_sum_exp(temp);
		}
	}

}
