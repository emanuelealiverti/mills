## Run posterior inference for the simplex factor model (Bhattacharya & Dunson, 2012, JASA)
require(rstan)
Y1 = read.table("../data1.txt")
H = 10
dataL = list(N = NROW(Y1),
	     P = NCOL(Y1),
	     D = length(unique(Y1[,1])),
	     H = H,
	     X = as.matrix(Y1))
set.seed(123)
m1 = stan(file = './SimplexFactor.stan',data = dataL,chains = 1)
tmp1 = extract(m1)
tmp1$lp__ = NULL
save(tmp1,file="SIM1_SF.RData")

set.seed(123)
Y2 = read.table("../data2.txt")
dataL$X = as.matrix(Y2)
dataL$N = NROW(Y2)
m2 = stan(file = './SimplexFactor.stan',data = dataL,chains = 1)
tmp2 = extract(m2)
tmp2$lp__ = NULL
save(tmp2,file="SIM2_SF.RData")

set.seed(123)
Y3 = read.table("../data3.txt")
dataL$X = as.matrix(Y3)
dataL$N = NROW(Y3)
m3 = stan(file = './SimplexFactor.stan',data = dataL,chains = 1)
tmp3 = extract(m3)
tmp3$lp__ = NULL
save(tmp3,file="SIM3_SF.RData")


set.seed(123)
Y4 = read.table("../data4.txt")
dataL$X = as.matrix(Y4)
dataL$N = NROW(Y4)
m4 = stan(file = './SimplexFactor.stan',data = dataL,chains = 1)
tmp4 = extract(m4)
tmp4$lp__ = NULL
save(tmp4,file="SIM4_SF.RData")

