# pdf of truncated Gamma with shape k and scale theta between a and b
dtgamma = function(x,k,theta,a,b) {
  if(all(x>=a & x<=b)) dgamma(x,k,1/theta)/(pgamma(b,k,1/theta) - pgamma(a,k,1/theta))
  else 0
}

# cdf of truncated Gamma with shape k and scale theta between a and b
ptgamma = function(x,k,theta,a,b) {
  if(all(x>=a & x<=b)) (pgamma(x,k,1/theta) - pgamma(a,k,1/theta))/(pgamma(b,k,1/theta) - pgamma(a,k,1/theta))
  else if(all(x<a)) 0
  else if(all(x>b)) 1
}

# Maximum likelihood estimator for k, theta and p0.
# Fit the first q portion of the data to truncated Gamma between 0 and q-th quantile
# q: specified quantile of s to fit the dtgamma
MLE = function(s,q=0.75){
  qs = quantile(s,q)[[1]]
  x = s[s<=qs]
  f0= fitdist(x, "tgamma", method="mle", start=list(k=0.5, theta=2), fix.arg=list(a=0,b=qs))
  k = f0$estimate["k"][[1]]
  theta = f0$estimate["theta"][[1]]
  p = q/pgamma(qs,k,scale=theta)
  list(k=k, theta=theta, p=p)
}

# Moment matching estimator for k, theta and p0, Schwartzman (2008).
# q: specified quantile of s for moment matching
# d: bin width for moment matching
MME = function(s,q=0.75,d=0.1)
{ 
  N = length(s)
  qs = quantile(s,q)[[1]]
  s1 = s[s<=qs]
  bre = seq(0,max(s1),by=d)
  bre = c(bre,bre[length(bre)]+d)
  h = hist(s1,breaks=bre,plot=F)
  c = h$mids
  y = h$counts
  fit = glm(y~c+log(c),family=poisson)
  coef = unname(fit$coefficients)
  k = coef[3]+1
  theta = -1/coef[2]
  p = exp(coef[1])*gamma(k)*theta^k/(N*d)
  list(k=k,theta=theta,p=p)
}

# Estimator for k and theta based on characteristic functions.
CFE = function(s,t,q=0.75){
  qs = quantile(s,q)[[1]]
  a11 = sum(cos(t*s))
  a12 = sum(s*cos(t*s))
  a21 = sum(sin(t*s))
  a22 = sum(s*sin(t*s))
  k = -t/(a11^2+a21^2)*((a11*a12+a21*a22)^2/(a21*a12-a11*a22)+a21*a12-a11*a22)
  theta = (a11*a22-a21*a12)/(t*(a11*a12+a21*a22))
  p = q/pgamma(qs,k,scale=theta)
  list(k=k,theta=theta,p=p)
}

# Characteristic function for Gamma(k,theta)
CF = function(t,k,theta){
  (1-theta*t*1i)^-k
}

# Derivative of characteristic function for Gamma(k,theta)
dCF = function(t,k,theta){
  1i*k*theta*(1-theta*t*1i)^(-k-1)
}

# Empirical characteristic function for Gamma statistics s
eCF=function(t,s){
  sapply(t,function(t) mean(exp(t*s*1i)))
}

# Derivative of empirical characteristic function for Gamma statistics s
deCF=function(t,s){
  sapply(t,function(t) mean(exp(t*s*1i)*s*1i))
}

# Smoothing the empirical characteristic function using local polynomial regression
# J: degree of polynomial
# b: bandwidth of kernel smoothing
smoother = function(s,J=3,b=0.1){
  x = seq(0,log(length(s))/2,by=0.01)
  y = sapply(x, eCF, s)
  dfit.re = locpoly(x,Re(y),drv=1,degree=J,bandwidth=b,gridsize=length(x))
  dfit.im = locpoly(x,Im(y),drv=1,degree=J,bandwidth=b,gridsize=length(x))
  return(list(x=x,dre=dfit.re$y,dim=dfit.im$y))
}

# Estimator using smoothed charcteristic function
# sm: object return by smoother() function
# t: optimal t value returned by t.opt() function
CFE.smooth = function(s,sm,t,q=0.75){
  qs = quantile(s,q)[[1]]
  ind = which(abs(sm$x-t)==min(abs(sm$x-t)))[1]
  psi = eCF(t,s)
  re = Re(psi)
  dre = sm$dre[ind]
  im = Im(psi)
  dim = sm$dim[ind]
  m = Mod(psi)
  dm = (re*dre+im*dim)/m
  ri = re*dim-dre*im
  theta = -m*dm/(t*ri)
  k = -t/m*(ri^2/(m^2*dm)+dm)
  p = q/pgamma(qs,k,scale=theta)
  list(k=k,theta=theta,p=p)
}

# Find optimal t for the characteristic function method, Jin & Cai (2007).
# Look for minimum t that satisfies mod(eCF(t,s)) = N^-r. 
# Usually there is only one solution for small r, but no solution for large r.
t.opt = function(s,r){
  N = length(s)
  a = 0
  b = log(N)
  m = N^-r
  while(abs(a-b)>0.01){
    x = (a+b)/2
    y = Mod(eCF(x,s))
    if(y<=m) b=x
    else a=x
  }
  return(x)
}

# Estimator of p0 by characteristic function, Jin & Cai (2007).
p.est = function(s,r,k,theta){
  N=length(s)
  Omega = function(t){
    fz = function(z) {
      psi = eCF(t*z,s)
      psi0 = CF(t*z,k,theta)
      Re(psi)*Re(psi0)*(1-abs(z))
    }
    integrate(fz,-1,1)$value
  }
  t = seq(0,sqrt(2*r*log(N)),by=0.05)
  ome = sapply(t,Omega)
  min(ome)
  
  #optimize(Omega,c(0,sqrt(2*r*log(N))))$objective
  #curve(Omega,0,sqrt(2*r*log(N)))
}

# Kullback-Leibler Divergence for the true Gamma null and the estimated null.
KLd = function(k0,theta0,k.hat,theta.hat){
  (k0-1)*digamma(k0) - log(theta0) - k0 - log(gamma(k0)) + log(gamma(k.hat)) + k.hat*log(theta.hat) - (k.hat-1)*(digamma(k0)+log(theta0)) + k0*theta0/theta.hat
}

