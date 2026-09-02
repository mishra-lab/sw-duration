# ==============================================================================
# config

source('utils.r')

l = list()
l$dur = 'Duration (years)'
l$prop = 'Population proportion'
l$fam = 'Family'
fl = list()
fl$fam = lapply(distrs,`[[`,'l')
fl$method = c('Complete','Turnover','2 × Med/Mean')
fl$data = list(x='Total: source [x]',za='Censored: sampled [z|s]')
fl$par = list(
  Ex  = 'E[x]',   CVx  = 'CV[x]',
  Ez  = 'E[z|s]', CVz  = 'CV[z|s]',
  a   = 'α',      la   = 'log(α)',
  b   = 'β',      lb   = 'log(β)',
  sp0 = 'ρₛ',     st95 = 'τ₉₅')

cmap = list()
cmap$fam = unlist(lapply(distrs,`[[`,'c'))
names(cmap$fam) = fl$fam

gen.makevars = function(opt){ cat(file='~/.R/Makevars',sep='',
  '\nCXX17=g++',
  '\nCXX17FLAGS=-march=native -mtune=native -fPIC -O',opt) }

dmax = 50
fams = names(distrs)
d.vec  = function(dd){ seq(dd/2,dmax,dd) }
p.x    = function(x,fam,a,b){     distrs[[fam]]$d(x=x,a=a,b=b,u=dmax) }
p.za   = function(z,fam,a,b){ 1 - distrs[[fam]]$p(q=z,a=a,b=b,u=dmax) }
p.sz   = function(z,sp0,st95){ 1+(sp0-1)*exp(z*log(.05/(1-sp0))/st95) }
i.near = function(x,ref){ which.min(abs(x-ref)) }

# ==============================================================================
# data

n.adj = function(n,m=1e3){
  # shrink effective n to simulate random effects
  n = round(1/(1/n+1/m))
}

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

model.data = function(Y,fam,gps=0,dd=.1,eps=1e-9,m=1e3){
  # TODO: support individual-level data?
  d_i = d.vec(dd)
  Ym  = subset(Y,meas=='mean')
  Yq  = subset(Y,meas=='q')
  q_i = sapply(Yq$value,i.near,ref=d_i)
  U = switch(fam, # bounds
    exp     = list(la=c(-eps,eps),lb=c(-1,3)),
    gamma   = list(la=c(-2,2),    lb=c(-3,3)),
    weibull = list(la=c(-1,2),    lb=c(-1,2)),
    lnorm   = list(la=c(-9,1),    lb=c(-2,3)),
    sbeta   = list(la=c(-3,3),    lb=c( 0,3)),
    skumar  = list(la=c(-3,3),    lb=c( 0,3)))
  data = list(
    Ula = U$la,          # bounds: full* duration log(a) param
    Ulb = U$lb,          # bounds: full* duration log(b) param
    Usp0 = c(.5,1),      # bounds: proportion sampled immediately
    Ust95 = c(.02,2),    # bounds: time to 95% sampled
    Ni = len(d_i),       # num integration steps
    Nm = nrow(Ym),       # num empiric means
    Nq = nrow(Yq),       # num empiric quantiles
    d_i = d_i,           # integration values
    dmax = dmax,         # maximum duration
    m_n = n.adj(Ym$n,m), # mean sample sizes
    m_v = Ym$value,      # mean estimates
    q_n = n.adj(Yq$n,m), # quantile sample sizes
    q_p = Yq$p,          # quantile CDF estimates
    q_i = q_i,           # quantile value index
    eps = eps,           # a small number
    gps = gps,           # generate posterior samples
    fam = which(fams==fam)) # distr family
}

get.sample = function(Y,fam,...,do='load',gps=0){
  model = load.txt('code','meta',ext='.stan')
  data = model.data(Y,...,fam=fam,gps=gps)
  args = list(data=data,pars=names(fl$par),
    chains=7,iter=1000,warm=500,seed=666)
  if (gps){ args$pars = c(args$pars,'m_vs','q_ps') }
  hash = hash.info(ulist(args,model=model))
  if (do=='load'){
    fit = load.rds('data','stan',hash,'fit')
  } else {
    fit = run.stan(args)
    if (do=='save'){
      args$data$fam = fams[args$data$fam] # HACK
      save.rds (fit, 'data','stan',hash,'fit')
      save.json(args,'data','stan',hash,'info')
      g = rstan::stan_trace(fit,inc_w=1,pars=names(fl$par)) + scale_colorfill()
      plot.save(g,'data','stan',hash,'trace',root='')
  }}
  S = expand.grid(i=args$warm+1:args$iter,chain=factor(1:args$chains))
  S = cbind(S,g=1:nrow(S),fam=fam,as.data.frame(fit))
}

run.stan = function(args){
  # TODO: diagnostic_file
  status(1,'run.stan: ',fams[args$data$fam])
  model = rstan::stan_model('meta.stan',auto_write=TRUE)
  args = ulist(args,object=model,cores=args$chains)
  fit = do.call(rstan::sampling,args)
}

plot.par = function(S,...){
  fl$par$lp__ = 'log(L)'
  S = melt(S,id=c('i','chain','fam'),var='par',m=names(fl$par))
  # print(t(sapply(split(S$value,S[c('fam','par')]),mci))) # NUM
  S$par = factor(S$par,names(fl$par),fl$par)
  S$fam = factor(S$fam,names(fl$fam),fl$fam)
  g = ggplot(S,aes(x='',y=value,color=fam,fill=fam)) +
    facet_wrap('par',scales='free',ncol=4) +
    scale_colorfill(v=cmap$fam) +
    geom_viola() +
    axis_blank('x') +
    labs(x='Parameter',y='Value',color=l$fam,fill=l$fam)
}

# ==============================================================================
# demo

plot.demo = function(S,Y,sub=10,dd=.01){
  d = d.vec(dd)
  D = rbind.lapply(seq(1,nrow(S),sub),function(g){
    pz = sum1(p.za(d,fam=S$fam[g],a=S$a[g],b=S$b[g])
            * p.sz(d,sp0=S$sp0[g],st95=S$st95[g]) ) / dd
    df(src='Model',d=d,strat=int.cut(d,Y$value),value=pz) })
  cols = str('q_ps[',1:nrow(Y),']')
  P = rbind.lapply(1:nrow(S),function(k){
    prop.diff(df(src='Model',n=NA,p=unlist(S[k,cols]),value=Y$value)) })
  Y = prop.diff(cbind(Y,src='Data'))
  f = c('(i) P[a < z < b | s]','(ii) P[z|s]'); f = factor(f,f,f)
  D$ff = f[2]; P$ff = f[1]; Y$ff = f[1];
  g = ggplot(D,aes(color=strat,fill=strat,lty=src)) +
    facet_wrap('ff',scales='free') +
    geom_mcrib(aes(x=d,y=value)) +
    geom_viola(data=P,aes(x=strat,y=p),color=NA,show.legend=FALSE) +
    geom_estimate(data=Y,aes(x=strat,y=p,ymin=plo,ymax=phi),width=.5) +
    ggh4x::scale_x_facet(ff==f[1],type='discrete') +
    ggh4x::scale_y_facet(ff==f[2],lim=c(0,.3)) +
    ggh4x::scale_x_facet(ff==f[2],breaks=c(0,Y$value,20),minor=0,lim=c(0,20)) +
    scale_colorfill('plasma') +
    scale_linetype_manual(values=c('21','solid')) +
    labs(x=l$dur,y=l$prop,color='Stratum',fill='Stratum',lty='Source')
}

main.demo = function(do='load'){
  Y = prop.ci(load.csv('data','Baral2014'))
  Y = subset(Y,ref=='rds')
  S = rbind.lapply(fams,get.sample,Y=Y,do=do,gps=1,.par=0)
  g = plot.par(S); plot.save(g,'stan','demo.pars',size=c(7,4))
  S = subset(S,fam=='gamma') # best fit
  g = plot.demo(S,Y); plot.save(g,'stan','demo.props',size=c(7,3))
}

# ==============================================================================
# meta

boot.mci = function(x,w=NULL,nb=1e5,seed=666){
  if (len(x)==0){ return(list(value=NA,lo=NA,hi=NA)) }
  set.seed(seed)
  xb = unlist(par.lapply(1:nb,function(i){ mean(sample(x,rep=1,p=w)) }))
  mci = list(value=mean(xb),lo=quantile(xb,.025,names=0),hi=quantile(xb,.975,names=0))
}

meta.classic = function(Y){
  Yt = subset(Y,F12=='turnover')
  Ym = subset(Y,F12=='median'|F12=='mean')
  # print(ulen(str(Ym$id,Ym$n))/nrow(Ym)) # DEBUG == 1
  M = rbind(
    df(ns=nrow(Yt),method='turnover',meta='raw.avg',boot.mci(1/Yt$p)),
    df(ns=nrow(Yt),method='turnover',meta='wtd.avg',boot.mci(1/Yt$p,Yt$n)),
    df(ns=nrow(Ym),method='2 × m',   meta='raw.avg',boot.mci(2*Ym$value)),
    df(ns=nrow(Ym),method='2 × m',   meta='wtd.avg',boot.mci(2*Ym$value,Ym$n)))
}

plot.distr = function(S,Y,sub=10,dd=.1,zoom=15){
  d = d.vec(dd)
  D = rbind.lapply(seq(1,nrow(S),sub),function(g){
    fam = S$fam[g]
    px = sum1(p.x (d,fam=fam,a=S$a[g],b=S$b[g])) / dd # TODO: why need sum1?
    pz = sum1(p.za(d,fam=fam,a=S$a[g],b=S$b[g])
            * p.sz(d,sp0=S$sp0[g],st95=S$st95[g]) ) / dd
    Di = cbind(fam=fam,rbind(
      df(d=d,data=fl$data$x, type='PDF',value=px),
      df(d=d,data=fl$data$za,type='PDF',value=pz),
      df(d=d,data=fl$data$x, type='CDF',value=cumsum(px)*dd),
      df(d=d,data=fl$data$za,type='CDF',value=cumsum(pz)*dd),
      df(d=0,data=fl$data$x, type='mean',value=xp.mean(d,px)),
      df(d=0,data=fl$data$za,type='mean',value=xp.mean(d,pz))))
  })
  geom_data = function(m,type,data,...){ geom_point(
    data=cbind(subset(Y,meas==m),type=type,data=data),
    inherit.aes=0,shape=16,alpha=.2,...) }
  fl$fam$data = 'Data'
  types = factor(c('CDF','PDF','mean'))
  D$type = factor(D$type,levels=types)
  D$fam = factor(D$fam,names(fl$fam),fl$fam)
  E = subset(D,type=='mean')
  D = subset(D,type!='mean' & d<=zoom)
  g = ggplot(D,aes(x=d,y=value,color=fam,fill=fam)) +
    facet_grid('type ~ data',scales='free') +
    geom_mcrib() +
    geom_viola(data=E,aes(x=value,y=''),show.legend=0,pos=dodge(w=-1)) +
    geom_data('q',   types[1],fl$data$za,map=aes(x=value,y=p,size=n.adj(n))) +
    geom_data('mean',types[3],fl$data$za,map=aes(y='',x=value,size=n.adj(n))) +
    scale_size_area(limits=c(10,1000),breaks=c(30,100,300,1000)) +
    scale_colorfill(v=cmap$fam) +
    scale_x_continuous(lim=c(0,zoom)) +
    ggh4x::scale_y_facet(type=='PDF',lim=0:1) +
    ggh4x::scale_y_facet(type=='mean',type='discrete') +
    ggh4x::force_panelsizes(rows=c(3,3,2)) +
    labs(x=l$dur,y=l$prop,color=l$fam,fill=l$fam,size='Data N*')
}

meta.forest = function(S0,pop='fsw'){
  T = subset(load.csv('data','Fazito2012'),kp==pop & region!='NorthAm') # TEMP
  S = dfu(aggregate(cbind(value=Ex)~kp+region+ns+fam,S0,mci,type=c))
  X = rbind(cbind(S,src='Ours',method=''),cbind(T,src='Fazito 2012',fam=T$method))
  X$fam = factor(X$fam,c(1:3,names(fl$fam)),c(fl$method,fl$fam))
  num.str = function(...){ gsub('NA','    •    ',sprintf(...)) }
  add_info = function(name,map,...){ list(
    geom_info(map,...),
    geom_info(ulist(map,label=name,y=''),data=df(region='',fam=NA),
      inherit.aes=0,check_overlap=TRUE,...))}
  geom_info = function(map,...){ geom_text(map=do.call(aes,map),
    size=3,family='Alegreya Sans',show.legend=FALSE,...) }
  g = ggplot(X,aes(x=value.m,xmin=value.lo,xmax=value.hi,
      y=interaction(src,method,fam),color=fam)) +
    facet_grid('region',scales='free',space='free') +
    geom_vline(xintercept=0,color='#ccc') +
    geom_estimate(size=1.5,lwd=.5,width=0,  position=dodge(w=1)) +
    add_info(hjust=1,'Source',  list(x=-6.5,label=quote(src))) +
    add_info(hjust=0,'Method / Family',list(x=-6.0,label=quote(fam))) +
    add_info(hjust=1,'Mean',    list(x=15,  label=quote(num.str('%.1f',value.m)))) +
    add_info(hjust=0,'(95% CI)',list(x=15.5,label=quote(num.str('(%.1f,%.1f)',value.lo,value.hi)))) +
    add_info(hjust=1, 'Ns',     list(x=20,label=quote(ns))) +
    scale_x_continuous(breaks=seq(0,12,2),lim=c(-10,20)) +
    scale_colorfill(v=cmap$fam,guide='none') +
    coord_cartesian(rev='y') +
    axis_blank('y',panel.grid.minor.x=blank) +
    labs(x=l$dur,y='',color=l$fam)
  plot.save(g,'stan','meta.forest',size=c(5,5))
}

main.meta = function(do='load',pop='fsw'){
  S0 = NULL
  Y0 = prop.ci(load.csv('data','Fazito2012x'))
  for (reg in c('Africa','Europe','LatAm','Asia')){
    Y = subset(Y0,kp==pop & region==reg & K26==1)
    S = rbind.lapply(fams,get.sample,Y=Y,do=do,.par=do=='load')
    g = plot.par(S);     plot.save(g,'stan',str('meta.pars.',  reg),size=c(7,4))
    g = plot.distr(S,Y); plot.save(g,'stan',str('meta.distrs.',reg),size=c(7,4))
    S0 = rbind(S0,cbind(S,kp=pop,region=reg,ns=ulen(Y$id)))
  }
  meta.forest(S0)
}

# ==============================================================================
# main

# gen.makevars(3)
# main.demo()
# main.meta()
