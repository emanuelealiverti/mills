# Composite mixture of log-linear models for categorical data

Multivariate categorical data are routinely collected in many application areas. As the number of cells in the table grows exponentially with the number of variables, many or even most cells will contain zero observations. This severe sparsity motivates appropriate statistical methodologies that effectively reduce the number of free parameters, with penalized log-linear models and latent structure analysis being popular options. This article proposes a fundamentally new class of methods, which we refer to as Mixture of Log Linear models (mills). Combining latent class analysis and log-linear models, mills defines a novel Bayesian methodology to model complex multivariate categorical with flexibility and interpretability. mills is shown to have key advantages over alternative methods for contingency tables in simulations and an application investigating the relation among suicide attempts and empathy.

This repo is associarted with the paper [Aliverti, E. and Dunson, D. Composite mixture of log-linear models for categorical data](arxiv.org/)

## Main content
`mills/` R package implementing the methods
`sim/` R code to reproduce the simulation studies. Require `rstan` to impement the competior approaches.


