
stickbreak = function(alpha,N){
##### draw w_1,w_2, ..., w_N from stick break dist'n
  V = rbeta(N,1,alpha)
  w = V
  d = 1
  for(i in 2:N) {
    d = d*(1-V[i-1])
    w[i] = V[i]*d
  }
  return(w/sum(w))
}

PriorSampleDir = function(alpha,mu,sigma){
### draw F ~ DP(alpha,F0) where F0 = N(mu,sigma)
  N = 1000
  s = rnorm(N,mu,sigma)
  w = stickbreak(alpha,N)

  par(mfrow=c(1,2))
  # plot(s,w,type="h",xlab="",ylab="weights",yaxt="n")
  plot(s,w,type="h",xlab="",ylab="weights")
  o = order(s)
  s = s[o]
  w = w[o]
  F = cumsum(w)
  plot(s,F,type="s",xlab="",ylab="F",lwd=2)
  lines(s,pnorm(s,mean=mu,sd=sigma),lwd=2,col="green")
}


PosteriorSampleDir = function(x,alpha,mu,sigma,M){
### draw M samples F ~ posterior
  ## x is the observed data
  ## first draw F ~ DP(alpha + n, Fbar)


  n = length(x)
  N = 100
  # F = matrix(0,N,M)
  # print(dim(F))
  grid = seq(-5,10,length=100)
  grid = c(x,grid)
  grid = sort(grid)
  ng   = length(grid)
  F = rep(0,ng)
  for(i in 1:ng){
    F[i] = (n/(alpha+n))*mean(x <= grid[i]) + (alpha/(alpha+n))*pnorm(grid[i],mu,sigma)
    ## the posterior F is a combination of empirical cdf and prior
  }
  plot(grid,F,type="l",lwd=3,col="blue",xlab="",ylab="F")
  Fbar = F

  ## below I add a plot of the prior
  N = 1000
  s = rnorm(N,mu,sigma)
  w = stickbreak(alpha,N)
  o = order(s)
  s = s[o]
  w = w[o]
  F = cumsum(w)
  lines(s,pnorm(s,mean=mu,sd=sigma),lwd=2,col="green")


  ## below we generate new data using Chinese Restaurant process
  ## this bypasses the generation of F from posterior (whose result is Fbar)
  for(i in 1:M){
    ## draw new N data points and visualize their empirical cdf
    z = rnorm(N,mu,sigma)
    u = rbinom(N,1,n/(alpha+n))
    y = sample(x,size=N,replace=TRUE)
    z[u==1] = y[u==1]
    s = z
    w = stickbreak(alpha+n,N)
    o = order(s)
    s = s[o]
    w = w[o]
    F = cumsum(w)
    lines(s,F,type="s")
  }
  return(list(Fbar=Fbar,grid=grid))
}

posterior.demo = function (alpha, n, M=10) {
  #pdf(file="DirProcess-plot-tmp.pdf",width=10,height=5)
  par(mfrow=c(2,2),mai=rep(.5,4))

  mu = 0
  sigma = 1
  PriorSampleDir(alpha,mu,sigma)

  x = rnorm(n,5,1)
  out = PosteriorSampleDir(x,alpha,0,1,M)
  grid=seq(-10,10,length=50)
  lines(grid,pnorm(grid,5,1),lwd=3,col="red",type="l")

  grid=seq(-5,10,length=50)
  alpha = .05
  eps = sqrt((1/(2*n))*log(2/alpha))


  F = rep(0,50)
  for(i in 1:50) {
    F[i] = mean(x <= grid[i])
  }
  U = F + eps
  L = F - eps
  U[U>1] = 1
  L[L<0] = 0
  plot(grid,pnorm(grid,5,1),lwd=3,col="red",type="l",xlab="",ylab="")
  lines(grid,F,type="s")
  lines(grid,U,type="s")
  lines(grid,L,type="s")
  polygon(c(grid,rev(grid)),c(L,rev(U)),col="mistyrose2")
  lines(grid,pnorm(grid,5,1),lwd=3,col="red")
  lines(grid,F,type="s",lwd=3)
  lines(grid,U,type="s")
  lines(grid,L,type="s")
  lines(out$grid,out$Fbar,lwd=3,col="blue")
  #dev.off()
}


#while(1) {PriorSampleDir(alpha=30, mu=0, sigma=1); Sys.sleep(1.5)}
while(1) {posterior.demo(alpha=30, n=100, M=20); Sys.sleep(1.5)}

