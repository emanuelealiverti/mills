# Simulation studies
This file illustrates how to reproduce results from [Aliverti, E. and Dunson, D. Composite mixture of log-linear models for categorical data](arxiv.org/), section 3.

## Set-up
We recommend to install the package locally, for example in `lib`.
Starting from the location of the repository (e.g. `~/GIT/mills/`
```
# if does
mkdir lib/

# compile attributes
R -q -e "setwd('mills'); Rcpp::compileAttributes(verbose=T)"

# fast build and install (no help or manuel for now)
R CMD build --no-build-vignettes --no-manual mills
R CMD INSTALL --no-docs --no-html -l lib mills
```

## Run simulations

Each folder contains different ordered scripts, supposed to be executed in order.

### Mills

```

01_MILLS/
├── 01_scen1.R
├── 02_scen2.R
├── 03_scen3.R
├── 04_scen4.R
└── 05_compute_biv.R
```

Each script `0%s_scen%s.R` simulates the data under scenario `%s`, save them as a `.txt` file and performs posterior sampling.
It can be executed, for example, using `Rscript 01_scen1.R`.
Posterior samples are saved in a separate `Scen%s.RData` files. 
Lastly `05_ComputeBiv.R` computes the posterior distribution of the estimated bivariate distributions, and summarise them via posterior mean and posterior $0.025$ and $0.975$ quantiles.
A one-line-command to run all the scripts  and compute the functional of interest is

```
for f in $(ls *.R); do Rscript $f;done
```

### Competitor

As a competitor approaches we have considered two popular latent structure models.
The methods are implemented in `stan` and, therefore, a working installation of the package [`rstan`](https://github.com/stan-dev/rstan/) is required.
Each file reads the respective `.txt` file (created by the associated `01_MILLS`
