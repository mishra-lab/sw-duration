# ==============================================================================
# config

source('utils.r')

l = list() # labels
l$qty = 'Duration'
l$pop = 'Population'
l$poq = 'Population (data)'
fl = list() # factor labels
fl$data = list(xo='x',xa='x|a',za='z|a',zs='z|s')
fl$qty = list(x='total',z='censored')
fl$pop = list(o='source',a='active',s='sampled')
fl$poq = list(o='source (total)',a='active (2 x censored)',s='sampled (2 x censored)')
fl$fam = lapply(distrs,`[[`,'l')

map = list() # aes maps
map$est = c(mean='solid',median='31',turnover='11')
map$qty = c(x='solid',z='21');           names(map$qty) = fl$qty[names(map$qty)]
map$pop = c(o='#000',a='#099',s='#f90'); names(map$pop) = fl$pop[names(map$pop)]
map$poq = c(o='#000',a='#099',s='#f90'); names(map$poq) = fl$poq[names(map$poq)]

plot.1o = list(w1=2,h1=1.4,wo=2,ho=1)

p.sam = function(z,sp0,st95){ 1+(sp0-1)*exp(z*log(.05/(1-sp0))/st95) }

# ==============================================================================
# toy simulation model

G0 = list(fam='weibull',m=5,cv=1,dmax=50,sp0=.5,st95=1)
S0 = df(est=NA,data=NA,d=NA,value=NA)

toy.math = function(fam,m,cv,dmax,sp0,st95,dd=1e-3,zm=1,pdf=0){
  if (fam=='exp' & cv != 1){ return(S0) }
  d = seq(dd/2,dmax,dd)
  dfam = distrs[[fam]]
  args = fit.distr(dfam$d,m=m,cv=cv,u=dmax,tol=1e-5)
  if (is.na(args$a)){ return(S0) }
  pxo = sum1(distr.call(dfam$d,args,x=d))     # total durs among source
  pxa = sum1(pxo * d)                         # total durs among active
  pza = sum1(1 - distr.call(dfam$p,args,q=d)) # censored durs among active
  pzs = sum1(pza * p.sam(d,sp0,st95))         # censored durs among sampled
  S = rbind(
    df(est='mean',    data='xo',d=NA,value=xp.mean(d,pxo)),
    df(est='mean',    data='xa',d=NA,value=xp.mean(d,pxa)),
    df(est='mean',    data='za',d=NA,value=xp.mean(d,pza)*zm),
    df(est='mean',    data='zs',d=NA,value=xp.mean(d,pzs)*zm),
    df(est='median',  data='xo',d=NA,value=xp.q(d,pxo,.50)),
    df(est='median',  data='xa',d=NA,value=xp.q(d,pxa,.50)),
    df(est='median',  data='za',d=NA,value=xp.q(d,pza,.50)*zm),
    df(est='median',  data='zs',d=NA,value=xp.q(d,pzs,.50)*zm),
    df(est='turnover',data='za',d=NA,value=1/sum(pza[d<1])),
    df(est='turnover',data='zs',d=NA,value=1/sum(pzs[d<1])))
  if (pdf){ S = rbind(S,
    df(est='pdf',     data='xo',d=d, value=pxo/dd),
    df(est='pdf',     data='xa',d=d, value=pxa/dd),
    df(est='pdf',     data='za',d=d, value=pza/dd),
    df(est='pdf',     data='zs',d=d, value=pzs/dd) )}
  return(S)
}

toy.stoc = function(fam,m,cv,dmax,sp0,st95,n=1e6,zm=1,pdf=0){
  if (fam=='exp' & cv != 1){ return(S0) }
  dfam = distrs[[fam]]
  args = fit.distr(dfam$d,m=m,cv=cv,u=dmax,tol=1e-5)
  if (is.na(args$a)){ return(S0) }
  xo = distr.call(dfam$r,args,n=n)             # total durs among source
  xa = sample(xo,n,rep=1,p=xo)                 # total durs among active
  za = runif(n,0,xa)                           # censored durs among active
  zs = sample(za,n,rep=1,p=p.sam(za,sp0,st95)) # censored durs among sampled
  d = seq(0,dmax,.2) # helper for pdf
  pfun = function(x){ h = hist(x,c(d,dmax+1)); list(d=h$mids,value=h$dens) }
  S = rbind(
    df(est='mean',    data='xo',d=NA,value=mean(xo)),
    df(est='mean',    data='xa',d=NA,value=mean(xa)),
    df(est='mean',    data='za',d=NA,value=mean(za)*zm),
    df(est='mean',    data='zs',d=NA,value=mean(zs)*zm),
    df(est='median',  data='xo',d=NA,value=median(xo)),
    df(est='median',  data='xa',d=NA,value=median(xa)),
    df(est='median',  data='za',d=NA,value=median(za)*zm),
    df(est='median',  data='zs',d=NA,value=median(zs)*zm),
    df(est='turnover',data='za',d=NA,value=1/mean(za<1)),
    df(est='turnover',data='zs',d=NA,value=1/mean(zs<1)))
  if (pdf){ S = rbind(S,
    df(est='pdf',     data='xo',pfun(xo)),
    df(est='pdf',     data='xa',pfun(xa)),
    df(est='pdf',     data='za',pfun(za)),
    df(est='pdf',     data='zs',pfun(zs)) )}
  return(S)
}

clean.toy = function(S){
  S = subset(S,!is.na(est))
  S$qty = factor(substr(S$data,1,1),names(fl$qty),fl$qty)
  S$pop = factor(substr(S$data,2,2),names(fl$pop),fl$pop)
  S$poq = factor(substr(S$data,2,2),names(fl$poq),fl$poq)
  S$data = factor(S$data,names(fl$data),fl$data)
  S$di = as.numeric(S$data)
  S$fam = factor(S$fam,names(fl$fam),fl$fam)
  S$CV = str('CV[x] = ',S$cv)
  return(S)
}

plot.toy.distr = function(f=2:5){
  G = ulist(G0,cv=c(.5,1,1.5),fam=names(fl$fam)[f])
  S = clean.toy(grid.apply(G,toy.math,pdf=1))
  S$leg = paste('p(',S$data,')')
  g = ggplot(S,aes(x=d,y=value,color=pop,lty=qty,shape=qty)) +
    facet_grid('fam ~ CV') +
    geom_line() +
    geom_point(data=~subset(.x,est=='mean'),aes(y=0,x=value)) +
    coord_cartesian(xlim=c(0,25),ylim=c(0,.3)) +
    scale_shape_manual   (values=c(1,4)) +
    scale_color_manual   (values=map$pop) +
    scale_linetype_manual(values=map$qty) +
    labs(x='Duration (years)',y='Density',
      color=l$pop,lty=l$qty,shape=l$qty)
  gl = inner_legend(g,subset(S,est=='mean' & fam==fl$fam$weib & cv==1.5))
  plot.save(gl,'toy','toy.distr')
  plot.save(gl + subset(S,fam==fl$fam$weib),'toy','toy.distr.1')
  b = S$fam==fl$fam$weib & S$est!='median'
  for (s in unique(S$di)){
    plot.save(g + subset(S, b & di<=s) + facet_grid('CV'),
      'toy',str('toy.distr.s',s),size=c(4,4)) }
}

plot.toy.est = function(f=2:5){
  G = ulist(G0,cv=seq(0,2,.1),fam=names(fl$fam)[f],zm=2)
  S = clean.toy(grid.apply(G,toy.math))
  S = subset(S,data!=fl$data$xa)
  g = ggplot(S,aes(x=cv,y=value,color=poq,lty=est)) +
    facet_wrap('fam',nrow=2) +
    geom_line() +
    coord_cartesian(ylim=c(0,20),xlim=c(0,2)) +
    scale_color_manual   (values=map$poq) +
    scale_linetype_manual(values=map$est) +
    labs(x='CV[x]',y='Ê[x]',lty='Estimate',color=l$poq)
  plot.save(g,'toy','toy.est')
  plot.save(g + subset(S,fam==fl$fam$weib),'toy','toy.est.1')
}

# ==============================================================================
# main

# plot.toy.distr()
# plot.toy.est()
