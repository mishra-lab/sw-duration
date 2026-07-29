# ==============================================================================
# utils

source('utils.r')

u.vec = function(du=.1,umax=50){ seq(du/2,umax,du) }
p.x   = function(x,m,cv){ het.funs$gamma$d(x,m,cv) }
p.za  = function(z,m,cv){ 1-het.funs$gamma$p(z,m,cv) }
p.sz  = function(z,ti=0){ 1-exp(z*log(.05)/ti) }
hist.pdf = function(x,b){ p = hist(x,c(-1e9,b,1e9))$dens[-1] }
i.near = function(x,ref){ which.min(abs(x-ref)) }

# ==============================================================================
# prettify

l = list() # labels
l$data = 'Duration: Pop'
l$dur = 'Duration (years)'

fl = list() # factor labels
fl$data = list(
  x  = '[x] total: source',
  xa = '[x|a] total: active',
  za = '[z|a] cens: active',
  zs = '[z|s] cens: sampled')
fl$fam = list(
  gamma   = '(a) Gamma',
  weibull = '(b) Weibull',
  lnorm   = '(c) Log-Normal')
fl$pars = list(
  ti  = 't 95% ID',
  Ex  = 'E[•]:[x]',
  CVx = 'CV[•]:[x]',
  Ez  = 'E[•]:[z|s]',
  CVz = 'CV[•]:[z|s]')

map = list()
map$data = c(x='#999',xa='#0cc',za='#066',zs='#f90')
map$est  = c(mean='solid',median='31',turnover='11')
names(map$data) = fl$data[names(map$data)]

# ==============================================================================
# toy simulation model

G0 = list(m=5,cv=0,n=1e5,fam='gamma',seed=666)

run.toy = function(m,cv,n,fam='gamma',ti=0,seed=NULL,du=.2,umax=25,pdf=FALSE){
  set.seed(seed)
  u  = u.vec(du,umax)                    # dummy durs
  x  = het.funs[[fam]]$r(n,m,cv)         # total durs among source
  xa = sample(x,n,rep=1,p=x)             # total durs among active
  za = runif(n,0,xa)                     # censored durs among active
  zs = sample(za,n,rep=1,p=p.sz(za,ti))  # censored durs among sampled
  S = rbind(
    df(est='mean',     data='x', value=mean(x)),
    df(est='mean',     data='xa',value=mean(xa)),
    df(est='mean',     data='za',value=mean(za)),
    df(est='mean',     data='zs',value=mean(zs)),
    df(est='median',   data='x', value=median(x)),
    df(est='median',   data='xa',value=median(xa)),
    df(est='median',   data='za',value=median(za)),
    df(est='median',   data='zs',value=median(zs)),
    df(est='turnover', data='za',value=1/mean(za<1)),
    df(est='turnover', data='zs',value=1/mean(zs<1)))
  if (pdf){ S = rbind(cbind(S,u=NA),
    df(est='pdf',u=u,  data='x', value=hist.pdf(x, u)),
    df(est='pdf',u=u,  data='xa',value=hist.pdf(xa,u)),
    df(est='pdf',u=u,  data='za',value=hist.pdf(za,u)),
    df(est='pdf',u=u,  data='zs',value=hist.pdf(zs,u))) }
  return(S)
}

clean.toy = function(S,z2m=FALSE){
  if (z2m){
    z = S$est!='turnover' & S$data %in% c('za','zs')
    S$value[z] = 2 * S$value[z]
  }
  S$data = factor(S$data,names(fl$data),fl$data)
  S$fam  = factor(S$fam,names(fl$fam),fl$fam)
  return(S)
}

plot.toy.est = function(){
  G = ulist(G0,cv=seq(0,3,.1),fam=names(fl$fam),ti=1,n=1e6)
  S = clean.toy(grid.apply(G,run.toy),z2m=TRUE)
  g = ggplot(S,aes(x=cv,y=value,color=data,lty=est)) +
    facet_grid(' ~ fam') +
    geom_line() +
    coord_cartesian(ylim=c(0,25)) +
    scale_color_manual(values=map$data) +
    scale_linetype_manual(values=map$est) +
    labs(x='CV[x]',y=l$dur,lty='Estimate',color=l$data)
  plot.save(g,'toy','est')
}

plot.toy.pdf = function(){
  G = ulist(G0,cv=c(.5,1,1.5),fam=names(fl$fam),ti=1,n=1e6)
  S = clean.toy(grid.apply(G,run.toy,pdf=TRUE))
  S$cv = str('(',c('i','ii','iii')[S$cv*2],') CV: ',S$cv)
  g = ggplot(S,aes(x=u,y=value,color=data)) +
    facet_grid('cv ~ fam') +
    geom_line() +
    geom_point(data=subset(S,est=='mean'),aes(y=0,x=value),shape=1) +
    coord_cartesian(xlim=c(0,25),ylim=c(0,.3)) +
    scale_color_manual(values=map$data) +
    labs(x=l$dur,y='Density',color=l$data)
  plot.save(g,'toy','pdf')
}

# ==============================================================================
# meta-analysis via stan

gen.makevars(opt=3)

gen.stan.data = function(Y,du=.1,umax=50,eps=1e-6){
  z_i = u.vec(du,umax)
  Ym  = subset(Y,meas=='mean')
  Yq  = subset(Y,meas=='q')
  q_i = sapply(Yq$value,i.near,ref=z_i)
  data = list(
    Ni = len(z_i),        # num integration steps
    Nm = nrow(Ym),        # num empiric means
    Nq = nrow(Yq),        # num empiric quantiles
    z_i = z_i,            # integration z values
    m_n = Ym$n,           # mean sample size
    m_v = Ym$value,       # mean estimate
    q_n = Yq$n,           # quantile sample sizes
    q_p = clip(Yq$p,.01), # quantile CDF estimate
    q_i = q_i,            # quantile value index
    eps = eps)            # a small number
}

get.stan.sample = function(Y,do='load',trace=TRUE){
  model = load.txt('code','meta',ext='.stan')
  args = list(data=gen.stan.data(Y),pars=names(fl$pars),
    chains=7,iter=1000,warm=500,seed=666)
  hash = hash.info(ulist(args,model=model))
  if (do=='load'){
    fit = load.rds('data','stan',hash,'fit')
  } else {
    fit = run.stan(args)
    if (do=='save'){
      save.rds (fit, 'data','stan',hash,'fit')
      save.json(args,'data','stan',hash,'info')
  }}
  if (trace){ plot.save('meta','trace',hash,g=rstan::stan_trace(fit)) }
  S = expand.grid(i=args$warm+1:args$iter,chain=factor(1:args$chains))
  S = cbind(S,as.data.frame(fit))
}

run.stan = function(args){
  model = rstan::stan_model('meta.stan',auto_write=TRUE)
  args = ulist(args,object=model,cores=args$chains)
  fit = do.call(rstan::sampling,args)
}

plot.stan.par = function(S){
  S$lp__ = NULL
  S = melt(S,id=c('chain','i'),var='par')
  S$par = factor(S$par,names(fl$pars),fl$pars)
  S[c('par','data')] = colsplit(S$par,':',0:1)
  g = ggplot(S,aes(x=data,y=value)) +
    facet_wrap('par',scales='free',space='free_x') +
    geom_violin(aes(color=chain),fill=NA,position='identity',scale='width') +
    coord_cartesian(ylim=c(0,NA)) +
    labs(x=l$data,y='Value')
}

plot.stan.pdf = function(S,du=.1,umax=50,sub=10){
  u  = u.vec(du,umax)
  S  = subset(S,i%%sub==0) # HACK
  ks = split(1:nrow(S),S[c('chain','i')])
  S = rbind.lapply(ks,function(k){
    px  = norm(p.x (u,S$Ex[k],S$CVx[k]))
    pzs = norm(p.za(u,S$Ex[k],S$CVx[k]) * p.sz(u,S$ti[k]))
    Sk = cbind(S[k,],rbind(
      df(u=u,par='PDF[•]',data='x', value=px),
      df(u=u,par='PDF[•]',data='zs',value=pzs),
      df(u=u,par='CDF[•]',data='x', value=cumsum(px)),
      df(u=u,par='CDF[•]',data='zs',value=cumsum(pzs))
    ))
  })
  S$data = factor(S$data,names(fl$data),fl$data)
  g = ggplot(S,aes(x=u,y=value,color=data,fill=data)) +
    facet_wrap('par',scales='free') +
    geom_summary(p=.95) +
    ggh4x::scale_y_facet(par=='PDF[•]',limits=c(0,.08)) +
    scale_color_manual(values=map$data) +
    scale_fill_manual(values=map$data) +
    labs(x=l$dur,y='Density',color=l$data,fill=l$data)
}

main.stan = function(popn='fsw',do='load'){
  Y = load.csv('data','Fazito2012')
  Y = subset(Y,pop==popn & method=='direct')
  S = get.stan.sample(Y,do=do)
  g = plot.stan.par(S); plot.save(g,'meta',str('par.',popn))
  g = plot.stan.pdf(S); plot.save(g,'meta',str('pdf.',popn))
}

# ==============================================================================
# main

# plot.toy.est()
# plot.toy.pdf()
# main.stan('fsw')
# main.stan('cli')
