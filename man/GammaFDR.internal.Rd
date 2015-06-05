\name{GammaFDR-internal}
\title{Internal GammaFDR functions}
\alias{CF}
\alias{dCF}
\alias{eCF}
\alias{deCF}
\alias{dtgamma}
\alias{ptgamma}
\alias{MLE}
\alias{MME}
\alias{CFE}
\alias{t.opt}
\alias{CFE.smooth}
\alias{smoother}
\alias{p.est}
\alias{KLd}

\description{Internal functions}

\usage{
CF(t,k,theta)
dCF(t,k,theta)
eCF(t,s)
deCF(t,s)
dtgamma(s,k,theta,a,b)
ptgamma(s,k,theta,a,b)
MLE(s,q=0.75)
MME(s,q=0.75,d=0.1)
CFE(s,t,q=0.75)
t.opt(s,r=0.05)
CFE.smooth(s,sm,t,q=0.75)
smoother(s,J=3,b=0.1)
p.est(s,r,k,theta)
KLd(k0,theta0,k.hat,theta.hat)
}

\author{Xing Ren}

\details{
  These are not to be called by the user.
}

\keyword{internal}
