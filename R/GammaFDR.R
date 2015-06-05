GammaFDR = function(s, method = c('MLE', 'MM', 'CF', 'CF.smooth'), 
                       q = 0.75, d = 0.1, r = 0.05, J = 3, b = 0.1, plot = F, fdr.cutoff = 0.01, breaks = 50){
  this.call = match.call()
  method = match.arg(method)
  if(any(s<0)) {
    print('Only non-negative values of s are allowed.')
    return
  }
  
  if(method=='MLE'){
    est = MLE(s,q)
  }
  
  if(method=='MM'){
    est = MME(s,q,d)
  }
  
  if(method=='CF'){
    t = t.opt(s,r)
    est = CFE(s,t)
  }
  
  if(method=='CF.smooth'){
    sm = smoother(s,J,b)
    t = t.opt(s,r)
    est = CFE.smooth(s,sm,t)
  }
  
  k0 = est$k
  theta0 = est$theta
  p0 = est$p
  fdr = sapply(s, function(x) p0*(1-pgamma(x,k0,scale=theta0))/(1-sum(s<x)/length(s)))
  s.cutoff = min(s[fdr<fdr.cutoff])  
  if(plot){
    hist(s,breaks)
    curve(p0*dgamma(x,k0,1/theta0),add==T,lty=2)
  }
  
  list(k0 = k0, theta0 = theta0, p0 = p0, fdr = fdr, method = method, call = this.call)
}