\name{GammaFDR}
\alias{GammaFDR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The main function of GammaFDR package
}
\description{
Estimate empirical null Gamma distribution parameters and false discovery rate.
}
\usage{
GammaFDR(s, method = c("MLE", "MM", "CF", "CF.smooth"), q = 0.75, d = 0.1, r = 0.05, J = 3, b = 0.1, plot = F, fdr.cutoff = 0.01, breaks = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{
 positive vector of Gamma mixture statistics
}
  \item{method}{
MLE: maximum likelihhod estimation. 
MME: mode matching.
CFE: characteristic function method.
CFE.smooth: smoothed characteristic function method.
}
  \item{q}{
 quantile of s for the zero interval region.
}
  \item{d}{
 bin width for mode matching method.
}
  \item{r}{
 tuning parameter for characteristic function method.
}
  \item{J}{
 degree of polynomials for smoothed characteristic function method.
}
  \item{b}{
 bin width for smoothed characteristic function method
}
  \item{plot}{
 display result in plot
}
  \item{fdr.cutoff}{
 false discovery rate
}
  \item{breaks}{
 number of breaks in histogram.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{ k0 }{ shape parameter of empirical null Gamma distribution }
  \item{ theta0 }{ scale parameter of empirical null Gamma distribution }
  \item{ p0 }{ proportion of null distribution in mixture data }
  \item{ fdr }{ false discovery rate}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
 Xing Ren
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
require(PoissonSeq)
require(fitdistrplus)
require(KernSmooth)

data(gender)
group = c(rep(1,17),rep(2,24))
geneID = gender[,1]
count = gender[,-1]
rownames(count) = geneID
data = list(n = count, y = group, type = "twoclass")
res = PS.Main(data)

s = res$tt
est = GammaFDR(s)
list(est$k0, est$theta0, est$p0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ empirical null distribution }
\keyword{ Gamma mixture }% __ONLY ONE__ keyword per line
