# Composite mixture of log-linear models for categorical data

This repository is associated with the paper [Aliverti, E. and Dunson, D. (2020) Composite mixture of log-linear models for categorical data.](https://arxiv.org/abs/2004.01462)

#### Abstract
Multivariate categorical data are routinely collected in many application areas. As the number of cells in the table grows exponentially with the number of variables, many or even most cells will contain zero observations. This severe sparsity motivates appropriate statistical methodologies that effectively reduce the number of free parameters, with penalized log-linear models and latent structure analysis being popular options. This article proposes a fundamentally new class of methods, which we refer to as Mixture of Log Linear models (mills). Combining latent class analysis and log-linear models, mills defines a novel Bayesian methodology to model complex multivariate categorical with flexibility and interpretability.



## Main contents
`mills/` R package implementing the methods, also includes different utilities for re-parametrising probability tensor into log-linear coefficients with corner parametrisation. Requires `RcppArmadillo` and `BayesLogit`.

`SIMULATIONS/` R code to reproduce the simulation studies included in Section 3 of the paper (specifically, Figure 2).
Simulations for the proposed approach require the included `mills` package. 
Simulations for the competitors require the [`rstan`](https://github.com/stan-dev/rstan/).


