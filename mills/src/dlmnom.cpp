#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//+++++++++++++++++++
//  Matrix subsetting
//+++++++++++++++++++

arma::colvec matrix_locs(arma::mat M, arma::umat locs) {
  
    arma::uvec eids = sub2ind( size(M), locs ); // Obtain Element IDs
    arma::vec v  = M.elem( eids );              // Values of the Elements
    
    return v;                               // Transpose to mimic R
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// General comment: there is no n\times k(\theta) since it computes the 
// likelihood for each observation!
// Comment 2: there is no 1+ in the normalising constant since the first
// parameter is already constrained at 0
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// [[Rcpp::export]]
arma::vec dlmnomPP(arma::mat Y_wp,  arma::mat theta_c){

	int n = Y_wp.n_rows;
	int p2 = Y_wp.n_cols;

	arma::mat like_2(n,p2,fill::zeros);
	arma::umat id_mat_2(p2,2);
	id_mat_2.col(0) = conv_to<arma::uvec>::from(linspace(0,p2-1,p2));

	arma::vec logL(n,fill::zeros);
	arma::vec theta_N(p2);

	for(int i = 0; i<n;i++){
		id_mat_2.col(1) = conv_to<arma::uvec>::from(Y_wp.row(i)-1);
		like_2(i) = accu(matrix_locs(theta_c, id_mat_2.t() ));
	}
	theta_N = sum(exp(theta_c),1);

	
	logL = sum(like_2,1) - accu(log(theta_N));


	return logL;
}
// [[Rcpp::export]]
arma::vec dlmnomPPw(arma::mat Y_wp,  arma::mat theta_c, arma::vec w){

	int n = Y_wp.n_rows;
	int p2 = Y_wp.n_cols;
	//arma::mat w_mat(theta_c.n_rows, theta_c.n_cols, fill::zeros);
	//theta_c.each_col() += w;

	arma::mat like_2(n,p2,fill::zeros);
	arma::umat id_mat_2(p2,2);
	id_mat_2.col(0) = conv_to<arma::uvec>::from(linspace(0,p2-1,p2));

	arma::vec logL(n,fill::zeros);
	arma::vec theta_N(p2);

	for(int i = 0; i<n;i++){
		id_mat_2.col(1) = conv_to<arma::uvec>::from(Y_wp.row(i)-1);
		like_2(i) = accu(w % matrix_locs(theta_c, id_mat_2.t() ));
	}
	theta_N = sum(exp(theta_c),1);
	
	logL = sum(like_2,1) - accu(w % log(theta_N));


	return logL;
}

// [[Rcpp::export]]
arma::vec dlmnom_sep(arma::mat Y_wp,  arma::mat theta_c){

	int n = Y_wp.n_rows;
	int p2 = Y_wp.n_cols;
	//arma::mat w_mat(theta_c.n_rows, theta_c.n_cols, fill::zeros);
	//theta_c.each_col() += w;

	arma::mat like_2(n,p2,fill::zeros);
	arma::umat id_mat_2(p2,2);
	id_mat_2.col(0) = conv_to<arma::uvec>::from(linspace(0,p2-1,p2));

	arma::vec logL(n,fill::zeros);
	arma::vec theta_N(p2);

	for(int i = 0; i<n;i++){
		id_mat_2.col(1) = conv_to<arma::uvec>::from(Y_wp.row(i)-1);
		like_2(i) = accu(matrix_locs(theta_c, id_mat_2.t() ));
	}
	theta_N = sum(exp(theta_c),1);
	
	logL = sum(like_2,1) - accu(log(theta_N));


	return logL;
}



