# ==============================================================================
# config

source('utils.r')

l = list() # labels
l$data = 'Duration: Pop'
l$dur  = 'Duration (years)'
fl = list() # factor labels
fl$data = list(
  x  = 'total (source) [x]',
  xa = 'total (active) [x|a]',
  za = 'trunc (active) [z|a]',
  zs = 'trunc (sampled) [z|s]')
fl$fam = list(
  gamma = 'Gamma',
  weibull = 'Weibull',
  lnorm = 'Log-Normal')

map = list() # aes maps
map$est  = c(mean='solid',median='31',turnover='11')
map$data = c(x='#999',xa='#3cc',za='#066',zs='#f90')
names(map$data) = fl$data[names(map$data)]

p.sam = function(z,sp0,st95){ 1+(sp0-1)*exp(z*log(.05/(1-sp0))/st95) }

# ==============================================================================
# toy simulation model

G0 = list(fam='gamma',m=5,cv=1,dmax=50,sp0=.5,st95=1)

run.toy = function(fam,m,cv,b,sp0,st95,n=1e5,zm=1,pdf=0){
  efun = function(...){ err <<- 1 }; err <<- 0 # classic R :)
  xfun = tryCatch(mcv.fun$r(fam,m=m,cv=cv,b=dmax),warning=efun,error=efun)
  if (err){ return(df(est=NA,data=NA,u=NA,value=NA)) }
  x  = xfun(n)                                 # total durs among source
  xa = sample(x,n,rep=1,p=x)                   # total durs among active
  za = runif(n,0,xa)                           # truncated durs among active
  zs = sample(za,n,rep=1,p=p.sam(za,sp0,st95)) # truncated durs among sampled
  S = rbind(
    df(est='mean',     data='x', d=NA,value=mean(x)),
    df(est='mean',     data='xa',d=NA,value=mean(xa)),
    df(est='mean',     data='za',d=NA,value=mean(za)*zm),
    df(est='mean',     data='zs',d=NA,value=mean(zs)*zm),
    df(est='median',   data='x', d=NA,value=median(x)),
    df(est='median',   data='xa',d=NA,value=median(xa)),
    df(est='median',   data='za',d=NA,value=median(za)*zm),
    df(est='median',   data='zs',d=NA,value=median(zs)*zm),
    df(est='turnover', data='za',d=NA,value=1/mean(za<1)),
    df(est='turnover', data='zs',d=NA,value=1/mean(zs<1)))
  if (pdf){
    pfun = function(x){ h = hist(x,c(d,dmax)); list(d=h$mids,value=h$dens) }
    d = seq(0,25,.2)
    S = rbind(S,
      df(est='pdf',data='x',  pfun(x)),
      df(est='pdf',data='xa', pfun(xa)),
      df(est='pdf',data='za', pfun(za)),
      df(est='pdf',data='zs', pfun(zs)) )}
  return(S)
}

clean.toy = function(S){
  S = subset(S,!is.na(est))
  S$data = factor(S$data,names(fl$data),fl$data)
  S$fam = factor(S$fam,names(fl$fam),fl$fam)
  S$CV = str('CV: ',S$cv)
  return(S)
}

plot.toy.distr = function(){
  G = ulist(G0,cv=c(.5,1,1.5),fam=names(fl$fam),n=1e6)
  S = clean.toy(grid.apply(G,run.toy,pdf=1))
  g = ggplot(S,aes(x=d,y=value,color=data)) +
    facet_grid('fam ~ CV') +
    geom_line() +
    geom_point(data=subset(S,est=='mean'),aes(y=0,x=value),shape=1) +
    coord_cartesian(xlim=c(0,25),ylim=c(0,.3)) +
    scale_color_manual(values=map$data) +
    labs(x=l$dur,y='Density',color=l$data)
  plot.save(g,'toy','toy.distrs')
}

plot.toy.est = function(){
  G = ulist(G0,cv=seq(0,2,.1),fam=names(fl$fam),n=1e6,zm=2)
  S = clean.toy(grid.apply(G,run.toy))
  g = ggplot(S,aes(x=cv,y=value,color=data,lty=est)) +
    facet_grid(' ~ fam') +
    geom_line() +
    coord_cartesian(ylim=c(0,20),xlim=c(0,2)) +
    scale_color_manual(values=map$data) +
    scale_linetype_manual(values=map$est) +
    labs(x='CV[x]',y=l$dur,lty='Estimate',color=l$data)
  plot.save(g,'toy','toy.est')
}

# ==============================================================================
# main

# plot.toy.distr()
# plot.toy.est()
