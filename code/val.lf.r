# validate likelihood functions for E[z] & P[L<z<U]
# given sample size of size ns from known p(z|a) <- p(x)
# we have 2 plot options: PDF or CDF

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
  cvx  = c(.5,1,1.5),   # CV[x]
  fam  = names(distrs)) # x ~ family

xmax = 50
bbs = list(c(0,1),c(1,5),c(5,9)) # prop bounds
nsam = c(30,100,300,1000) # sample sizes
nrep = 1e3                # num replicates
uvec = seq(0,1,.01)       # dummy prop variable
zvec = uvec * G$mx * 1.2  # dummy mean variable

# data gen ---------------------------------------------------------------------

data.gen = function(n,mx,cvx,du,fam){
  args = fit.distr(fam,m=mx,cv=cvx,u=xmax)
  x = distr.call(distrs[[fam]]$r,args,n=n)
  z = runif(n,0,wtd.sample(x,n,w=x))
}

# f.mean = function(m,cv,ns){ dgamma(zvec,n/cv^2,n/cv^2/m) } # PDF
# f.prop = function(m,   ns){ dbeta(uvec,ns*m,ns*(1-m)) } # PDF
f.mean = function(m,cv,ns){ pgamma(zvec,ns/cv^2,ns/cv^2/m) } # CDF
f.prop = function(m,   ns){ pbeta(uvec,ns*m,ns*(1-m)) } # CDF
stat.distr = function(z,ns,nr){
  mz = mean(z); cvz = sd(z)/mz
  zr = lapply(1:nr,function(r){ sample(z,ns) }) # replicate samples
  Xs = rbind(
    df(stat='mean',data='analytic',  u=zvec,value=f.mean(mz,cvz,ns)),
    df(stat='mean',data='stochastic',u=NA,value=vapply(zr,mean,0)),
    rbind.lapply(bbs,function(bb){ rbind(
      df(stat=labb(bb),data='analytic',  u=uvec,value=f.prop(pbb(z,bb),ns)),
      df(stat=labb(bb),data='stochastic',u=NA,value=vapply(zr,pbb,0,bb=bb)))
  }))
}

X = grid.apply(G,function(...){
  z = data.gen(...) # true z data
  Xi = grid.apply(list(ns=nsam),function(ns){
    Xs = stat.distr(z,ns=ns,nr=nrep) # replicates & analytic distribution
}) })
X$fam = factor(X$fam,names(distrs),lapply(distrs,`[[`,'l'))
X$cvx = str('CV[x]: ',X$cvx)

# plots ------------------------------------------------------------------------

plot.1o = list(w1=1.5,h1=1,wo=1.5,ho=1) # plot sizing
for (s in unique(X$stat)){
  Xf = subset(X, stat==s & data=='analytic')
  Xr = subset(X, stat==s & data=='stochastic')
  g = ggplot(map=aes(color=as.factor(ns),lty=data)) +
    facet_grid('fam ~ cvx') +
    geom_line(data=Xf,aes(x=u,y=value)) +
    # geom_line(data=Xr,aes(x=value),stat='density',alpha=.5) + # PDF
    geom_step(data=Xr,aes(x=value),stat='ecdf',alpha=.5) + # CDF
    scale_colorfill('plasma') +
    scale_linetype_manual(values=c('22','solid')) +
    coord_cartesian(xlim=range(Xr$value)) +
    # scale_y_continuous(trans='sqrt') + # PDF
    labs(x='Value',color='Ni',lty='Source',
      # y='Density') # PDF
      y='Cumulative Density') # CDF
  plot.save(g,'val',str('val.lf.',s))
}
