# ==============================================================================
# config

source('utils.r')

l = list()
l$dur = 'Duration (years)'
l$prop = 'Population proportion'
fl = list()
fl$fams = lapply(distrs,`[[`,'l')
fl$pars = list(
  a   = 'α',      la   = 'log(α)',
  b   = 'β',      lb   = 'log(β)',
  Ex  = 'E[x]',   CVx  = 'CV[x]',
  Ez  = 'E[z|s]', CVz  = 'CV[z|s]',
  sp0 = 'p0',     st95 = 't95')

gen.makevars = function(opt){ cat(file='~/.R/Makevars',sep='',
  '\nCXX17=g++',
  '\nCXX17FLAGS=-march=native -mtune=native -fPIC -O',opt) }

dmax   = 50
d.vec  = function(dd){ seq(dd/2,dmax,dd) }
p.x    = function(x,fam,a,b){     distrs[[fam]]$d(x=x,a=a,b=b,u=dmax) }
p.za   = function(z,fam,a,b){ 1 - distrs[[fam]]$p(q=z,a=a,b=b,u=dmax) }
p.sz   = function(z,sp0,st95){ 1+(sp0-1)*exp(z*log(.05/(1-sp0))/st95) }
i.near = function(x,ref){ which.min(abs(x-ref)) }

# ==============================================================================
# data

prop.ci = function(Y,eps=1e-3){
  # generate phi & plo from p,n
  Y$p   = clip(Y$p,eps)
  Y$plo = qbeta(.025,Y$n*Y$p,Y$n*(1-Y$p))
  Y$phi = qbeta(.975,Y$n*Y$p,Y$n*(1-Y$p))
  return(Y)
}

prop.diff = function(Y){
  # convert cum probs -> probs & add strat
  Y = rbind(Y,dfu(Y[1,],value=Inf,p=1))
  Y = dfu(Y,p=diff(c(0,Y$p)),strat=int.cut(Y$value,Y$value[-nrow(Y)]))
  Y = prop.ci(Y)
}

# ==============================================================================
# model

model.data = function(Y,fam,gps=0,dd=.1,eps=1e-9){
  d_i = d.vec(dd)
  Ym  = subset(Y,meas=='mean')
  Yq  = subset(Y,meas=='q')
  q_i = sapply(Yq$value,i.near,ref=d_i)
  U = switch(fam, # bounds
    gamma   = list(la=c( -1,3),lb=c(-1,3)),
    weibull = list(la=c( -1,3),lb=c(-1,3)),
    lnorm   = list(la=c(-20,1),lb=c(-2,3)),
    sbeta   = list(la=c( -3,3),lb=c(0,3)),
    skumar  = list(la=c( -3,3),lb=c(0,3)))
  data = list(
    Ula = U$la,       # bounds: full* duration log(a) param
    Ulb = U$lb,       # bounds: full* duration log(b) param
    Usp0  = c(.5,1),  # bounds: proportion sampled immediately
    Ust95 = c(.02,3), # bounds: time to 95% sampled
    Ni = len(d_i),  # num integration steps
    Nm = nrow(Ym),  # num empiric means
    Nq = nrow(Yq),  # num empiric quantiles
    d_i = d_i,      # integration values
    dmax = dmax,    # maximum duration
    m_n = Ym$n,     # mean sample size
    m_v = Ym$value, # mean estimate
    q_n = Yq$n,     # quantile sample sizes
    q_p = Yq$p,     # quantile CDF estimate
    q_i = q_i,      # quantile value index
    eps = eps,      # a small number
    gps = gps,      # generate posterior samples
    fam = switch(fam,gamma=1,weibull=2,lnorm=3,sbeta=4,skumar=5))
}

get.sample = function(Y,fam,...,do='load',gps=0){
  model = load.txt('code','meta',ext='.stan')
  args = list(data=model.data(Y,...,fam=fam,gps=gps),pars=names(fl$pars),
    chains=7,iter=1000,warm=500,seed=666)
  if (gps){ args$pars = c(args$pars,'m_vs','q_ps') }
  hash = hash.info(ulist(args,model=model))
  if (do=='load'){
    fit = load.rds('data','stan',hash,'fit')
  } else {
    fit = run.stan(args)
    if (do=='save'){
      save.rds (fit, 'data','stan',hash,'fit')
      save.json(args,'data','stan',hash,'info')
      g = rstan::stan_trace(fit,inc_w=1)
      plot.save(g,'data','stan',hash,'trace',root='')
  }}
  S = expand.grid(i=args$warm+1:args$iter,chain=factor(1:args$chains))
  S = cbind(S,g=1:nrow(S),fam=fam,as.data.frame(fit))
}

run.stan = function(args){
  model = rstan::stan_model('meta.stan',auto_write=TRUE)
  args = ulist(args,object=model,cores=args$chains)
  fit = do.call(rstan::sampling,args)
}

plot.par = function(S,...){
  fl$pars$lp__ = 'log(L)'
  S = melt(S,id=c('i','chain','fam'),var='par',m=names(fl$pars))
  S$par = factor(S$par,names(fl$pars),fl$pars)
  S$fam = factor(S$fam,names(fl$fams),fl$fams)
  g = ggplot(S,aes(x='',y=value,color=fam)) +
    facet_wrap('par',scales='free',ncol=4) +
    geom_viola() +
    labs(x='Parameter',y='Value',fill='Family')
}

# ==============================================================================
# demo

plot.demo = function(S,Y,dd=.01){
  d = d.vec(dd)
  D = df(d=d,strat=int.cut(d,Y$value),value=sum1(
    p.za(d,fam=S$fam[1],a=mean(S$a),b=mean(S$b)) *
    p.sz(d,sp0=mean(S$sp0),st95=mean(S$st95)) )/dd )
  cols = str('q_ps[',1:nrow(Y),']')
  P = rbind.lapply(1:nrow(S),function(k){
    prop.diff(df(n=NA,p=unlist(S[k,cols]),value=Y$value)) })
  f = c('p[z|a]','p[L < z < U]'); f = factor(f,f,f)
  D$ff = f[1]; P$ff = f[2]; Y$ff = f[2];
  g = ggplot(D,aes(color=strat,fill=strat)) +
    facet_wrap('ff',scales='free') +
    geom_ribbon(aes(x=d,ymin=0,ymax=value),color=NA,alpha=.5) +
    geom_line(aes(x=d,y=value)) +
    geom_viola(data=P,aes(x=strat,y=p),color=NA,show.legend=0) +
    geom_estimate(data=prop.diff(Y),aes(x=strat,y=p,ymin=plo,ymax=phi),width=.5) +
    ggh4x::scale_x_facet(ff==f[1],breaks=c(0,Y$value,25),minor=0,lim=c(0,25)) +
    ggh4x::scale_x_facet(ff==f[2],type='discrete') +
    scale_clr_viridis(discrete=TRUE,option='plasma') +
    labs(x=l$dur,y=l$prop)
}

main.demo = function(do='load'){
  Y = prop.ci(load.csv('data','Baral2014'))
  Y = subset(Y,ref=='rds')
  S = rbind.lapply(names(distrs),get.sample,Y=Y,do=do,gps=1)
  g = plot.par(S); plot.save(g,'stan','demo','par',size=c(7,5))
  S = subset(S,fam=='skumar') # best fit
  g = plot.demo(S,Y); plot.save(g,'stan','demo','prop',size=c(7,3))
}

# ==============================================================================
# meta

boot.mci = function(x,w=NULL,nb=1e5,seed=666){
  set.seed(seed)
  xb = unlist(par.lapply(1:nb,function(i){ mean(sample(x,rep=1,p=w)) }))
  mci = list(value=mean(x),lo=quantile(x,.025,names=0),hi=quantile(x,.975,names=0))
}

meta.classic = function(Y){
  Yt = subset(Y,F12=='turnover')
  Ym = subset(Y,F12=='median'|F12=='mean')
  # Yr = subset(Y,F12=='retired') # TODO
  M = rbind(
    df(ns=nrow(Yt),method='turnover',meta='raw.avg',boot.mci(1/Yt$p)),
    df(ns=nrow(Yt),method='turnover',meta='wtd.avg',boot.mci(1/Yt$p,Yt$n)),
    df(ns=nrow(Ym),method='2 × m',   meta='raw.avg',boot.mci(2*Ym$value)),
    df(ns=nrow(Ym),method='2 × m',   meta='wtd.avg',boot.mci(2*Ym$value,Ym$n)))
}

plot.distr = function(S,Y,sub=100,fam='gamma',dd=.05,dmax=50){
  d = d.vec(dd,dmax)
  D = rbind.lapply(seq(1,nrow(S),sub),function(g){
    px = p.x (d,fam=fam,m=S$Eu[g],cv=S$CVu[g],b=50)
    pz = sum1(p.za(d,fam=fam,m=S$Eu[g],cv=S$CVu[g],b=50)
            * p.sz(d,sp0=S$sp0[g],st95=S$st95[g]) ) / dd
    # print(c(Ex=S$Ex[g],eEx=sum(d*px)/sum(px), # DEBUG
    #         Ez=S$Ez[g],eEz=sum(d*pz)/sum(pz)))
    Di = rbind(
      df(d=d,data='total (source) p[x]',  type='PDF',value=px),
      df(d=d,data='trunc (active) p[z|a]',type='PDF',value=pz),
      df(d=d,data='total (source) p[x]',  type='CDF',value=cumsum(px)*dd),
      df(d=d,data='trunc (active) p[z|a]',type='CDF',value=cumsum(pz)*dd))
  })
  # TODO: fix fam
  g = ggplot(D,aes(x=d,y=value,color=data,fill=data)) +
    facet_wrap('type',scales='free') +
    stat_summary(geom='ribbon',color=NA,alpha=.5,
      fun.min=function(x){ quantile(x,.025) },
      fun.max=function(x){ quantile(x,.975) }) +
    stat_summary(geom='line',fun='mean') +
    scale_clr_manual(values=c('#999','#066')) +
    ggh4x::scale_x_facet(type=='PDF',limits=c(0,25)) +
    ggh4x::scale_y_facet(type=='PDF',limits=c(0,.5)) +
    labs(x=l$dur,y=l$prop,color='Durations',fill='Durations')
}

main.meta = function(do='load'){
  Y = prop.ci(load.csv('data','Fazito2012'))
  Yi = subset(Y,kp=='fsw' & region=='Africa')
  M  = meta.classic(Yi)
  Yi = subset(Yi,K26==1) # TODO: clean-up
  S = get.sample(Yi,do=do)
  g = plot.par(S);      plot.save(g,'stan','meta.par',size=c(5,3))
  g = plot.distr(S,Yi); plot.save(g,'stan','meta.distr',size=c(7,3))
}

# ==============================================================================
# main

# gen.makevars(3)
# main.demo()
# main.meta()
