# ==============================================================================
# config

source('utils.r')

l = list() # labels
l$data = 'Duration: Population'
l$dur  = 'Duration (years)'
fl = list() # factor labels
fl$data = list(
  x  = 'total: source [x]',
  xa = 'total: active [x|a]',
  za = 'trunc: active [z|a]',
  zs = 'trunc: sampled [z|s]')
fl$fam = lapply(distrs,`[[`,'l')

map = list() # aes maps
map$est  = c(mean='solid',median='31',turnover='11')
map$data = c(x='#999',xa='#3cc',za='#066',zs='#f90')
names(map$data) = fl$data[names(map$data)]

plot.1o$wo = 2 # legend is wide

p.sam = function(z,sp0,st95){ 1+(sp0-1)*exp(z*log(.05/(1-sp0))/st95) }

# ==============================================================================
# toy simulation model

G0 = list(fam='weibull',m=5,cv=1,dmax=99,sp0=.5,st95=1)
S0 = df(est=NA,data=NA,d=NA,value=NA)

toy.math = function(fam,m,cv,dmax,sp0,st95,dd=1e-3,zm=1,pdf=0){
  if (fam=='exp' & cv != 1){ return(S0) }
  d = seq(dd/2,dmax,dd)
  dfam = distrs[[fam]]
  args = fit.distr(dfam$d,m=m,cv=cv,u=dmax,tol=1e-5)
  if (is.na(args$a)){ return(S0) }
  px  = sum1(distr.call(dfam$d,args,x=d))     # total durs among source
  pxa = sum1(px * d)                          # total durs among active
  pza = sum1(1 - distr.call(dfam$p,args,q=d)) # truncated durs among active
  pzs = sum1(pza * p.sam(d,sp0,st95))         # truncated durs among sampled
  S = rbind(
    df(est='mean',     data='x', d=NA,value=xp.mean(d,px)),
    df(est='mean',     data='xa',d=NA,value=xp.mean(d,pxa)),
    df(est='mean',     data='za',d=NA,value=xp.mean(d,pza)*zm),
    df(est='mean',     data='zs',d=NA,value=xp.mean(d,pzs)*zm),
    df(est='median',   data='x', d=NA,value=xp.quant(d,px)),
    df(est='median',   data='xa',d=NA,value=xp.quant(d,pxa)),
    df(est='median',   data='za',d=NA,value=xp.quant(d,pza)*zm),
    df(est='median',   data='zs',d=NA,value=xp.quant(d,pzs)*zm),
    df(est='turnover', data='za',d=NA,value=1/sum(pza[d<1])),
    df(est='turnover', data='zs',d=NA,value=1/sum(pzs[d<1])))
  if (pdf){
    S = rbind(S,
      df(est='pdf',data='x',  d=d,value=px /dd),
      df(est='pdf',data='xa', d=d,value=pxa/dd),
      df(est='pdf',data='za', d=d,value=pza/dd),
      df(est='pdf',data='zs', d=d,value=pzs/dd) )}
  return(S)
}

toy.stoc = function(fam,m,cv,dmax,sp0,st95,n=1e6,zm=1,pdf=0){
  if (fam=='exp' & cv != 1){ return(S0) }
  dfam = distrs[[fam]]
  args = fit.distr(dfam$d,m=m,cv=cv,u=dmax,tol=1e-5)
  if (is.na(args$a)){ return(S0) }
  x  = distr.call(dfam$r,args,n=n)             # total durs among source
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
    pfun = function(x){ h = hist(x,c(d,dmax+1)); list(d=h$mids,value=h$dens) }
    d = seq(0,dmax,.2)
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
  S$CV = str('CV[x] = ',S$cv)
  return(S)
}

plot.toy.distr = function(){
  G = ulist(G0,cv=c(.5,1,1.5),fam=names(fl$fam))
  S = clean.toy(grid.apply(G,toy.math,pdf=1))
  g = ggplot(S,aes(x=d,y=value,color=data)) +
    facet_grid('fam ~ CV') +
    geom_line() +
    geom_point(data=~subset(.x,est=='mean'),aes(y=0,x=value),shape=1) +
    coord_cartesian(xlim=c(0,25),ylim=c(0,.3)) +
    scale_color_manual(values=map$data) +
    labs(x=l$dur,y='Density',color=l$data)
  plot.save(g,'toy','toy.distr')
  plot.save(g + subset(S,fam==fl$fam$weib),'toy','toy.distr.1')
}

plot.toy.est = function(f=1:6){
  G = ulist(G0,cv=seq(0,2,.1),fam=names(fl$fam)[f],zm=2)
  S = clean.toy(grid.apply(G,toy.math))
  g = ggplot(S,aes(x=cv,y=value,color=data,lty=est)) +
    facet_wrap('fam',nrow=2) +
    geom_line() +
    coord_cartesian(ylim=c(0,20),xlim=c(0,2)) +
    scale_color_manual(values=map$data) +
    scale_linetype_manual(values=map$est) +
    labs(x='CV[x]',y=l$dur,lty='Estimate',color=l$data)
  plot.save(g,'toy','toy.est')
  plot.save(g + subset(S,fam==fl$fam$weib),'toy','toy.est.1')
}

# ==============================================================================
# main

# plot.toy.distr()
# plot.toy.est()
