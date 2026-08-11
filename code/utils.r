library('reshape2')
options(width=180)

# ==============================================================================
# misc

len = length
str = paste0
df = data.frame
.verb = 3

status = function(lvl,...){
  if (lvl > .verb | lvl <= 0){ return() }
  pre = list(c(rep('-',80),'\n'),'',' > ')[[lvl]]
  cat(pre,...,'\n',sep='')
}

ulist = function(x=list(),xu=list(),...){
  x = c(x,xu,list(...))
  x[!duplicated(names(x),fromLast=TRUE)]
}

dfu = function(X,...){ as.data.frame(ulist(X,...)) }
sum1 = function(x){ x/sum(x) }
clip = function(x,eps){ pmin(1-eps,pmax(eps,x)) }
ulen = function(x){ len(unique(x)) }

int.cut = function(x,br,lo=0,hi=Inf){
  labels = gsub('-Inf','+',str(c(lo,br),'-',c(br,hi)))
  x.cut = cut(x,breaks=c(lo,br,hi),labels=labels)
}

# ==============================================================================
# files + i/o

proj.root = strsplit(file.path(getwd(),''),file.path('','code',''))[[1]][1]

root.path = function(...,ext='',create=FALSE){
  path = str(file.path(proj.root,...),ext)
  if (create & !dir.exists(dirname(path))){
    dir.create(dirname(path),recursive=TRUE) }
  return(path)
}

load.txt = function(...,ext='.txt'){
  fname = root.path(...,ext=ext)
  status(3,'load: ',fname)
  readLines(con=fname)
}

save.txt = function(X,...,ext='.txt'){
  fname = root.path(...,ext=ext,create=TRUE)
  status(3,'save: ',fname)
  cat(X,file=fname,sep='\n')
}

load.csv = function(...,ext='.csv'){
  fname = root.path(...,ext=ext)
  status(3,'load: ',fname)
  read.csv(file=fname,fileEncoding='Latin1')
}

save.csv = function(X,...,ext='.csv'){
  fname = root.path(...,ext=ext,create=TRUE)
  status(3,'save: ',fname)
  write.csv(X,file=fname,row.names=FALSE)
}

load.rds = function(...,ext='.rds'){
  fname = root.path(...,ext=ext)
  status(3,'load: ',fname)
  readRDS(fname)
}

save.rds = function(X,...,ext='.rds'){
  fname = root.path(...,ext=ext,create=TRUE)
  status(3,'save: ',fname)
  saveRDS(X,file=fname)
}

load.json = function(...,ext='.json'){
  fname = root.path(...,ext=ext)
  status(3,'load: ',fname)
  rjson::fromJSON(file=fname)
}

save.json = function(X,...,ext='.json',indent=2){
  fname = root.path(...,ext=ext,create=TRUE)
  status(3,'save: ',fname)
  write(rjson::toJSON(X,indent=indent),file=fname)
}

hash.info = function(info,.len=11){
  hash = substr(digest::sha1(info),1,.len)
}

# ==============================================================================
# lapply

.cores = 7

par.lapply = function(...,.par=TRUE){
  if (.par && .cores > 1){
    parallel::mclapply(...,mc.cores=.cores)
  } else {
    lapply(...)
  }
}

rbind.lapply = function(...){
  do.call(rbind,c(par.lapply(...)))
}

grid.apply = function(x,fun,args=list(),...,.rbind=TRUE,.cbind=TRUE,.par=TRUE){
  xg = expand.grid(x,stringsAsFactors=FALSE)
  gi = seq(nrow(xg))
  g.args   = lapply(gi,function(i){ ulist(as.list(xg[i,,drop=FALSE]),args,...) })
  g.fun    = ifelse(.cbind,function(...){ cbind(fun(...),...) },fun)
  g.lapply = ifelse(.rbind,rbind.lapply,par.lapply)
  g.lapply(g.args,do.call,what=g.fun,.par=.par)
}

# ==============================================================================
# stats

d.mean = function(dfun,u,...,eps=1e-7){
  f = function(x){ x * dfun(x,u=u,...) }
  m = integrate(f,lower=eps,upper=u-eps)$value
}

d.cv = function(dfun,u,...,m,eps=1e-7){
  if (missing(m)){ m = d.mean(dfun,u=u,...) }
  f = function(x){ (x-m)^2 * dfun(x,u=u,...) }
  cv = sqrt(integrate(f,lower=eps,upper=u-eps)$value) / m
}

xp.mean = function(x,p){
  m = sum(p*x) / sum(p)
}

xp.cv = function(x,p,m){
  if (missing(m)){ m = xp.mean(x,p) }
  cv = sqrt(sum(p*(x-m)^2)/sum(p)) / m
}

xp.quant = function(x,p,po=.5){
  spline(x=cumsum(sum1(p)),y=x,xout=po)$y
}

fit.n = function(p,p.025,p.975){
  err.fun = function(x){
    n = 10^x
    q2 = qbeta(c(.025,.975),n*p,n*(1-p))
    err = sum((q2-c(p.025,p.975))^2)
  }
  n = 10^optimize(err.fun,c(0,4))$minimum
}

rtru = LaplacesDemon::rtrunc
dtru = LaplacesDemon::dtrunc
ptru = LaplacesDemon::ptrunc

distrs = list(
  exp = list(
    r = function(n,a,b,u){ rtru(n=n,spec='exp',rate=1/a,a=0,b=u) },
    d = function(x,a,b,u){ dtru(x=x,spec='exp',rate=1/a,a=0,b=u) },
    p = function(q,a,b,u){ ptru(x=q,spec='exp',rate=1/a,a=0,b=u) },
    c = '#cc0033', l = 'Exponential'),
  gamma = list(
    r = function(n,a,b,u){ rtru(n=n,spec='gamma',shape=a,rate=b,a=0,b=u) },
    d = function(x,a,b,u){ dtru(x=x,spec='gamma',shape=a,rate=b,a=0,b=u) },
    p = function(q,a,b,u){ ptru(x=q,spec='gamma',shape=a,rate=b,a=0,b=u) },
    c = '#990099', l = 'Gamma'),
  weibull = list(
    r = function(n,a,b,u){ rtru(n=n,spec='weibull',shape=a,scale=b,a=0,b=u) },
    d = function(x,a,b,u){ dtru(x=x,spec='weibull',shape=a,scale=b,a=0,b=u) },
    p = function(q,a,b,u){ ptru(x=q,spec='weibull',shape=a,scale=b,a=0,b=u) },
    c = '#0099cc', l = 'Weibull'),
  lnorm = list(
    r = function(n,a,b,u){ rtru(n=n,spec='lnorm',meanlog=a,sdlog=b,a=0,b=u) },
    d = function(x,a,b,u){ dtru(x=x,spec='lnorm',meanlog=a,sdlog=b,a=0,b=u) },
    p = function(q,a,b,u){ ptru(x=q,spec='lnorm',meanlog=a,sdlog=b,a=0,b=u) },
    c = '#00cc66', l = 'Log-Normal'),
  sbeta = list(
    r = function(n,a,b,u){ rbeta(n=n,  shape1=a,shape2=b) * u },
    d = function(x,a,b,u){ dbeta(x=x/u,shape1=a,shape2=b) / u },
    p = function(q,a,b,u){ pbeta(q=q/u,shape1=a,shape2=b) },
    c = '#ffcc00', l = 'Scaled Beta'),
  skumar = list(
    r = function(n,a,b,u){ extraDistr::rkumar(n=n,  a=a,b=b) * u },
    d = function(x,a,b,u){ extraDistr::dkumar(x=x/u,a=a,b=b) / u },
    p = function(q,a,b,u){ extraDistr::pkumar(q=q/u,a=a,b=b) },
    c = '#ff6600', l = 'Scaled Kumar'))

fit.distr = function(dfun,m,cv,u,tol=1e-7,eps=1e-7){
  if (is.character(dfun)){ dfun = distrs[[dfun]]$d }
  jfun = function(x){
    mx  = d.mean(dfun,a=exp(x[1]),b=exp(x[2]),u=u,eps=eps)
    cvx = d.cv  (dfun,a=exp(x[1]),b=exp(x[2]),u=u,eps=eps,m=m)
    err = (m-mx)^2 + (cv-cvx)^2 }
  fail = function(...){ failed <<- 1 }; failed <<- 0
  opt = tryCatch(optim(c(a=0,b=0),jfun),warning=fail,error=fail)
  if (failed){ return(list(a=NA,b=NA,u=NA)) }
  if (opt$value > tol){ print(opt) }
  return(c(as.list(exp(opt$par)),u=u))
}

distr.call = function(fun,args,...){
  do.call(fun,c(args,list(...)))
}

# ==============================================================================
# ggplot2

library('ggplot2')

theme_set(theme_light())
theme_update(
  text=element_text(family='Alegreya Sans'),
  strip.background=element_rect(fill='#eeeeee'),
  strip.text.x=element_text(color='black'),
  strip.text.y=element_text(color='black'),
  legend.spacing.y=unit(-1,'mm'))

dodge = function(...,w=.5){ position_dodge(...,width=w) }

geom_viola = function(...){
  geom_violin(...,alpha=.5,scale='width')
}

geom_estimate = function(...,shape=18,width=1){
  suppressWarnings(list(
    geom_point(...,size=1,shape=shape),
    geom_errorbar(...,lwd=1/3,width=width)
))}

scale_colorfill = function(v='plasma',...,d=1,end=.85){
  require('viridis')
  if (len(v) > 1){ return(list(
    scale_color_manual(values=v,...),
    scale_fill_manual (values=v,...) ))}
  if (v %in% ls('package:viridis')){ return(list(
    scale_color_viridis(option=v,...,discrete=d,end=end),
    scale_fill_viridis (option=v,...,discrete=d,end=end) ))}
  if (v %in% names(RColorBrewer::brewer.pal.info)){ return(ifelse(d,
    scale_color_brewer   (pal=v,aes=c('color','fill'),...),
    scale_color_distiller(pal=v,aes=c('color','fill'),...) ))}
}

plot.save = function(g,...,size=NULL,ext='.pdf',root='out/fig'){
  if (is.null(size)){ size = plot.size(g) }
  if (ext=='.pdf'){ dev = cairo_pdf } else { dev = NULL }
  fname = root.path(root,...,ext=ext,create=TRUE)
  status(3,'save: ',fname)
  ggsave(plot=g,file=fname,w=size[1],h=size[2],device=dev)
}

plot.1o = list(w1=2,h1=1.5,wo=1.5,ho=1)

plot.size = function(g,...){
  s = ulist(plot.1o,...)
  layout = ggplot_build(g)$layout$layout
  size = c(w=s$wo+s$w1*max(layout$COL),
           h=s$ho+s$h1*max(layout$ROW))
}
