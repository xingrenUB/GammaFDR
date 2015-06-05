# GammaFDR
Estimating empirical null distribution of Gamma statistics

install.packages("GammaFDR_1.0.tar.gz", repos = NULL, type = "source")

require(PoissonSeq)
require(fitdistrplus)
require(KernSmooth)

data(gender)
group = c(rep(1,17), rep(2,24))

geneID = gender[,1]
count = gender[,-1]
rownames(count) = geneID
data = list(n = count, y = group, type = "twoclass")
res = PS.Main(data)

s = res$tt
est = GammaFDR(s)
list(est$k0, est$theta0, est$p0)

