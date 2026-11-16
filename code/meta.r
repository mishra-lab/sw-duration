# ==============================================================================
# config

source('utils.r')
distrs$lnorm = NULL # HACK
distrs$resam = list(l='(Resample)',c='#000000')

l = list()
l$dur = 'Duration (years)'
l$fam = 'Family'
fl = list()
fl$fam = lapply(distrs,`[[`,'l')
fl$region = c(Africa='Africa',Asia='Asia',Europe='Europe',LatAm='Latin\nAm',NorAm='North\nAm')
fl$method = list(`1`='M1: former',`2`='M2: 1 / P[z < 1]',`3`='M3: mean × 2')
fl$data = list(x='Total among source (x)',za='Censored among sampled (z|s)')
fl$par = list(
  Ex  = 'E[x]',   CVx  = 'CV[x]',
  Ez  = 'E[z|s]', CVz  = 'CV[z|s]',
  a   = 'α',      la   = 'log(α)',
  b   = 'β',      lb   = 'log(β)',
  sp0 = 'ρₛ',     st95 = 'τ₉₅')

cmap = list()
cmap$sam = unlist(setNames(lapply(distrs,`[[`,'c'),fl$fam))

gen.makevars = function(opt){ cat(file='~/.R/Makevars',sep='',
  '\nCXX17=g++',
  '\nCXX17FLAGS=-march=native -mtune=native -fPIC -O',opt) }

dmax = 50
fams = names(distrs)[-6] # HACK
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
    sbeta   = list(la=c(-3,3),    lb=c( 0,3)),
    skumar  = list(la=c(-3,3),    lb=c( 0,3)),
    lnorm   = list(la=c(-9,1),    lb=c(-2,3)))
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

get.sample = function(Y,fam,reg,...,do='load',gps=0){
  model = load.txt('code','meta',ext='.stan')
  data = model.data(Y,...,fam=fam,gps=gps)
  args = list(data=data,pars=names(fl$par),
    chains=7,iter=1500,warm=500,seed=666)
  if (gps){ args$pars = c(args$pars,'m_vs','q_ps') }
  hash = hash.info(ulist(args,model=model))
  path = file.path('data','stan',str(reg,'.',fam),hash)
  if (do=='load'){
    fit = load.rds(path,'fit')
  } else {
    fit = run.stan(args)
    if (do=='save'){
      args$data$reg = reg
      args$data$fam = fam
      save.rds (fit, path,'fit')
      save.json(args,path,'info')
      g = rstan::stan_trace(fit,inc_w=1,pars=names(fl$par)) + scale_colorfill()
      plot.save(g,path,'trace',root='')
  }}
  S = expand.grid(i=1:(args$iter-args$warm),chain=factor(1:args$chains))
  S = cbind(S,g=1:nrow(S),fam=fam,region=reg,as.data.frame(fit))
  S = subset(S,!( # HACK
    (region=='Asia'   & fam=='skumar'  & chain==5) |
    (region=='Europe' & fam=='weibull' & chain==3) |
    (region=='Europe' & fam=='skumar'  & chain==1) |
    # TODO: reg=='Europe' & fam=='lnorm' ?
  0 ))
}

run.stan = function(args){
  # TODO: diagnostic_file
  status(1,'run.stan: ',fams[args$data$fam])
  model = rstan::stan_model('meta.stan',auto_write=TRUE)
  args = ulist(args,object=model,cores=args$chains)
  fit = do.call(rstan::sampling,args)
}

ffl = function(x,fli){ factor(x,names(fli),fli) }

resample = function(S,n=7000,seed=666){
  set.seed(seed)
  p = exp(S$lp__-max(S$lp__))
  R = S[sample(1:nrow(S),n,rep=1,p=p),]
  # print(table(R$fam)) # NUM
  S = rbind(cbind(S,resam=FALSE,sam=ffl(S$fam,  fl$fam)),
            cbind(R,resam=TRUE, sam=ffl('resam',fl$fam)))
}

plot.par = function(S,...){
  fl$par$lp__ = 'log(L)'
  S[S$resam,c('a','la','b','lb')] = NA
  S = melt(S,id=c('i','chain','fam','sam'),var='par',m=names(fl$par))
  S$par = ffl(S$par,fl$par)
  g = ggplot(S,aes(x='',y=value,color=sam,fill=sam)) +
    facet_wrap('par',scales='free',ncol=4) +
    scale_colorfill(v=cmap$sam) +
    geom_viola(pos=dodge(w=1,pre='single')) +
    axis_blank('x') +
    labs(x='Parameter',y='Value',color=l$fam,fill=l$fam)
}

# ==============================================================================
# demo

prop.demo = function(S,Y,dd=.01){
  d = d.vec(dd)
  s = int.cut(d,Y$value)
  P = rbind.lapply(1:nrow(S),function(k){
    px = p.x(d,S$fam[k],S$a[k],S$b[k]) * dd # HACK
    Pi = aggregate(p~strat,df(p=px,strat=s),sum)
  })
  print(aggregate(p~strat,P,mci,type=c)) # NUM
}

plot.demo = function(S,Y,sub=10,dd=.01){
  d = d.vec(dd)
  D = rbind.lapply(seq(1,nrow(S),sub),function(g){
    px = sum1(p.x (d,fam=S$fam[g],a=S$a[g],b=S$b[g])) / dd
    pz = sum1(p.za(d,fam=S$fam[g],a=S$a[g],b=S$b[g])
            * p.sz(d,sp0=S$sp0[g],st95=S$st95[g]) ) / dd
    Di = rbind(
      df(src='Model (z|s)',d=d,strat=int.cut(d,Y$value),value=pz),
      df(src='Model (x)',d=d,strat=NA,value=px) )})
  cols = str('q_ps[',1:nrow(Y),']')
  P = rbind.lapply(1:nrow(S),function(k){
    prop.diff(df(src='Model (z|s)',n=NA,p=unlist(S[k,cols]),value=Y$value)) })
  Y = prop.diff(cbind(Y,src='Data (z|s)'))
  f = c('(i) P[a < z < b | s; Nᵢ ]','(ii) f(z|s), f(x)'); f = factor(f,f,f)
  D$ff = f[2]; P$ff = f[1]; Y$ff = f[1];
  g = ggplot(D,aes(color=strat,fill=strat,lty=src)) +
    facet_wrap('ff',scales='free') +
    geom_mcrib(aes(x=d,y=value)) +
    geom_viola(data=P,aes(x=strat,y=p),color=NA,show.legend=FALSE) +
    geom_estimate(data=Y,aes(x=strat,y=p,ymin=plo,ymax=phi),width=.5) +
    ggh4x::scale_x_facet(ff==f[1],type='discrete') +
    ggh4x::scale_y_facet(ff==f[2],lim=c(0,.3)) +
    ggh4x::scale_x_facet(ff==f[2],breaks=c(0,Y$value,20),minor=0,lim=c(0,20)) +
    scale_linetype_manual(values=c('11','31','solid')) +
    scale_colorfill('plasma',na.value='#ccc',
      labels=function(x){ ifelse(is.na(x),'Overall',x) }) +
    labs(x=l$dur,y='Value',color='Stratum',fill='Stratum',lty='Source')
}

main.demo = function(do='load'){
  Y = prop.ci(load.csv('data','Baral2014'))
  Y = subset(Y,ref=='rds')
  S = resample(rbind.lapply(fams,get.sample,Y=Y,reg='demo',do=do,gps=1,.par=do=='load'))
  g = plot.par(S); plot.save(g,'stan','demo.pars',size=c(7,4))
  S = subset(S,resam)
  g = plot.demo(S,Y); plot.save(g,'stan','demo.props',size=c(7,3))
  prop.demo(S,Y)
}

# ==============================================================================
# meta

plot.distr = function(S,Y,sub=10,dd=.1,zoom=15){
  d = d.vec(dd)
  D = rbind.lapply(seq(1,nrow(S),sub),function(g){
    px = sum1(p.x (d,fam=S$fam[g],a=S$a[g],b=S$b[g])) / dd # TODO: why need sum1?
    pz = sum1(p.za(d,fam=S$fam[g],a=S$a[g],b=S$b[g])
            * p.sz(d,sp0=S$sp0[g],st95=S$st95[g]) ) / dd
    cpx = cumsum(px)*dd
    cpz = cumsum(pz)*dd
    Di = cbind(fam=S$fam[g],sam=S$sam[g],rbind(
      df(d=d,data=fl$data$x, type='PDF',value=px),
      df(d=d,data=fl$data$za,type='PDF',value=pz),
      df(d=d,data=fl$data$x, type='CDF',value=cpx),
      df(d=d,data=fl$data$za,type='CDF',value=cpz),
      df(d=0,data=fl$data$x, type='p<1',value=approx(d,cpx,1)$y),
      df(d=0,data=fl$data$za,type='p<1',value=approx(d,cpz,1)$y),
      df(d=0,data=fl$data$x, type='mean',value=xp.mean(d,px)),
      df(d=0,data=fl$data$za,type='mean',value=xp.mean(d,pz))))
  })
  # print(aggregate(value~sam+data,subset(D,type=='p<1'),mci,type=c)) # NUM
  D = subset(D,type!='p<1')
  geom_data = function(m,type,data,...){ geom_point(
    data=cbind(subset(Y,meas==m),type=type,data=data),
    inherit.aes=0,shape=16,alpha=.2,...) }
  fl$fam$data = 'Data'
  types = factor(c('CDF','PDF','mean'))
  D$type = factor(D$type,levels=types)
  E = subset(D,type=='mean')
  D = subset(D,type!='mean' & d<=zoom)
  g = ggplot(D,aes(x=d,y=value,color=sam,fill=sam)) +
    facet_grid('type ~ data',scales='free') +
    geom_mcrib() +
    geom_viola(data=E,aes(x=value,y=''),show.legend=0,pos=dodge(w=-1)) +
    geom_data('q',   types[1],fl$data$za,map=aes(x=value,y=p,size=n.adj(n))) +
    geom_data('mean',types[3],fl$data$za,map=aes(y='',x=value,size=n.adj(n))) +
    scale_size_area(limits=c(10,1000),breaks=c(30,100,300,1000)) +
    scale_colorfill(v=cmap$sam) +
    scale_x_continuous(lim=c(0,zoom),breaks=seq(0,zoom,2)) +
    ggh4x::scale_y_facet(type=='PDF',lim=0:1) +
    ggh4x::scale_y_facet(type=='mean',type='discrete') +
    ggh4x::force_panelsizes(rows=c(3,3,2)) +
    labs(x=l$dur,y='Density',color=l$fam,fill=l$fam,size='Data N\'')
}

meta.forest = function(SA,pop='fsw'){
  T = subset(load.csv('data','Fazito2012'),kp==pop)
  S = dfu(aggregate(cbind(value=Ex)~kp+region+ns+sam,SA,mci,type=c))
  X = rbind(cbind(T,src='Fazito 2012',sam=ffl(T$method,fl$method)),
            cbind(S,src='Model',method=''))
  X$region = ffl(X$region,fl$region)
  num.str = function(...){ gsub('NA','    •    ',sprintf(...)) }
  add_info = function(name,map,...){ list(
    geom_info(map,...),
    geom_info(ulist(map,label=name,y=''),data=df(region='',fam=NA),
      inherit.aes=0,check_overlap=TRUE,...))}
  geom_info = function(map,...){ geom_text(map=do.call(aes,map),
    size=3,family='Alegreya Sans',show.legend=FALSE,...) }
  g = ggplot(X,aes(x=value.m,xmin=value.lo,xmax=value.hi,
      color=sam,y=interaction(src,method,sam))) +
    facet_grid('region',scales='free',space='free') +
    geom_vline(xintercept=0,color='#ccc') +
    geom_estimate(size=1.5,lwd=.5,width=0,  position=dodge(w=1)) +
    add_info(hjust=1, 'Source',  list(x=-6.5,label=quote(src))) +
    add_info(hjust=0, 'Method / Family',list(x=-6.0,label=quote(sam))) +
    add_info(hjust=1, 'Mean',    list(x=15,  label=quote(num.str('%.2f',value.m)))) +
    add_info(hjust=.5,'(95% CI)',list(x=17.5,label=quote(num.str('(%.2f, %.2f)',value.lo,value.hi)))) +
    add_info(hjust=1, 'Ns',      list(x=21,  label=quote(ns))) +
    scale_x_continuous(breaks=seq(0,12,2),lim=c(-10,21)) +
    scale_colorfill(v=cmap$sam,guide='none',na.value='#999') +
    coord_cartesian(rev='y') +
    axis_blank('y',panel.grid.minor.x=blank) +
    theme(panel.spacing=unit(4,'pt')) +
    labs(x=l$dur,y=NULL)
  plot.save('stan','meta.forest.all',size=c(5,6),g=g)
  plot.save('stan','meta.forest.rs',size=c(5,3),g=g +
    subset(X,sam==fl$fam$resam|src=='Fazito 2012'))
}

main.meta = function(do='load',pop='fsw'){
  SA = NULL
  YA = prop.ci(load.csv('data','Fazito2012x'))
  for (reg in names(fl$region)){
    Y = subset(YA,kp==pop & region==reg & K26==1)
    S = resample(rbind.lapply(fams,get.sample,Y=Y,reg=reg,do=do,.par=do=='load'))
    g = plot.par(S);     plot.save(g,'stan',str('meta.pars.',  reg),size=c(7,4))
    g = plot.distr(S,Y); plot.save(g,'stan',str('meta.distrs.',reg),size=c(7,4))
    SA = rbind(SA,cbind(S,kp=pop,ns=ulen(Y$id)))
  }
  save.rds(SA,'data','stan','meta')
  meta.forest(S0)
}

# ==============================================================================
# main

# gen.makevars(3)
# main.demo()
# main.meta()
# TODO: csv export
