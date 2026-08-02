# validate sampling distributions for E[z] & p(L<z<U)
# given sample size of size ns from known p(z|a) <- p(x)

# utils ------------------------------------------------------------------------

source('utils.r')
set.seed(666)

pbb = function(x,bb){ mean(bb[1] < x & x < bb[2]) }
labb = function(bb){ str('p',bb[1],'-',bb[2]) }
wtd.sample = function(x,n,w){ # faster
  # https://stackoverflow.com/a/15205104/5228288
  x[order(runif(len(x))^(1/w),decreasing=TRUE)][1:n] }

# config -----------------------------------------------------------------------

G = list( # source pop total duration
  n    = 1e5,           # source pop size
  mx   = 5,             # E[x]
  cvx  = c(.5,1,1.5,3), # CV[x]
  xmax = 50,            # max(x)
  fam  = c('gamma','lognormal')) # x ~ family

bbs = list(c(0,1),c(1,2),c(0,5),c(1,5),c(1,9),c(5,9)) # prop bounds
nsam = c(30,100,300,1000) # sample sizes
nrep = 1e4                # num replicates
uvec = seq(0,1,.01)       # dummy prop variable
zvec = uvec * G$mx * 1.2  # dummy mean variable

# data gen ---------------------------------------------------------------------

data.gen = function(n,mx,cvx,xmax,du,fam){
  x = pmin(xmax,het.funs[[fam]]$r(n,mx,cvx))
  z = runif(n,0,wtd.sample(x,n,w=x))
}

f.mean = function(m,cv,ns){ het.funs$gamma$d(zvec,m,cv/sqrt(ns)) }
f.prop = function(m,   ns){ dbeta(uvec,ns*m,ns*(1-m)) }
meas.se = function(z,ns,nr){
  mz = mean(z); cvz = sd(z)/mz
  zr = lapply(1:nr,function(r){ sample(z,ns) }) # replicate samples
  Xs = rbind(
    df(meas='mean',data='math',u=zvec,value=f.mean(mz,cvz,ns)),
    df(meas='mean',data='sample',u=NA,value=vapply(zr,mean,0)),
    rbind.lapply(bbs,function(bb){ rbind(
      df(meas=labb(bb),data='math',u=uvec,value=f.prop(pbb(z,bb),ns)),
      df(meas=labb(bb),data='sample',u=NA,value=vapply(zr,pbb,0,bb=bb)))
  }))
}

X = grid.apply(G,function(...){
  z = data.gen(...) # true z data
  Xi = grid.apply(list(ns=nsam),function(ns){
    Xs = meas.se(z,ns=ns,nr=nrep) # replicates & math distribution
},.par=0) })

# plots ------------------------------------------------------------------------

for (m in unique(X$meas)){
  Xm = subset(X, meas==m & data=='math')
  Xs = subset(X, meas==m & data=='sample')
  g = ggplot(map=aes(color=as.factor(ns),lty=data)) +
    facet_grid('fam ~ cvx') +
    geom_line(data=Xm,aes(x=u,y=value)) +
    geom_line(data=Xs,aes(x=value),stat='density',alpha=.5) +
    scale_clr_viridis(discrete=1,option='plasma') +
    scale_linetype_manual(values=c('22','solid')) +
    coord_cartesian(xlim=range(Xs$value)) +
    scale_y_continuous(trans='sqrt') +
    labs(title=m,x='Value',y='Cumulative Density',color='Ns',lty='Data')
  plot.save(g,'sam.distr',m)
}
