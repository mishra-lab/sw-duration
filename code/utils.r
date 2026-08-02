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

dfu = function(X,...){
  as.data.frame(ulist(X,...))
}

hash.info = function(info,.len=11){
  hash = substr(digest::sha1(info),1,.len)
}

sum1 = function(x){ x/sum(x) }
clip = function(x,eps){ pmin(1-eps,pmax(eps,x)) }

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

fit.n = function(p,p.025,p.975){
  # X$n.eff = mapply(fit.n,X$p.adj,X$p.025,X$p.975)
  err.fun = function(x){
    n = 10^x
    q2 = qbeta(c(.025,.975),n*p,n*(1-p))
    err = sum((q2-c(p.025,p.975))^2)
  }
  n = 10^optimize(err.fun,c(0,4))$minimum
}

fit.weibull = function(m,cv2,...){
  efun = function(k){ s = gamma(1+1/k)^2; e = ((gamma(1+2/k)-s)/s-cv2)^2 }
  k = optimize(efun,c(1e-6,1e+6))$minimum
  par = list(shape=k,scale=m/gamma(1+1/k),...)
}

het.funs = list(
  # m = mean; het = CV (sd / mean)
  gamma = list(
    l = 'Gamma',
    r = function(n,m,het){ cv2 = max(het^2,1e-9); rgamma(n,shape=1/cv2,scale=m*cv2) },
    d = function(x,m,het){ cv2 = max(het^2,1e-9); dgamma(x,shape=1/cv2,scale=m*cv2) },
    p = function(q,m,het){ cv2 = max(het^2,1e-9); pgamma(q,shape=1/cv2,scale=m*cv2) },
    q = function(p,m,het){ cv2 = max(het^2,1e-9); qgamma(p,shape=1/cv2,scale=m*cv2) }),
  weibull = list(
    l = 'Weibull',
    r = function(n,m,het){ f = fit.weibull(m,het^2); rweibull(n,shape=f$shape,scale=f$scale) },
    d = function(x,m,het){ f = fit.weibull(m,het^2); dweibull(x,shape=f$shape,scale=f$scale) },
    p = function(q,m,het){ f = fit.weibull(m,het^2); pweibull(q,shape=f$shape,scale=f$scale) },
    q = function(p,m,het){ f = fit.weibull(m,het^2); qweibull(p,shape=f$shape,scale=f$scale) }),
  lognormal = list(
    l = 'Log-Norm',
    r = function(n,m,het){ u = log(m/sqrt(1+het^2)); s = sqrt(log(1+het^2)); rlnorm(n,meanlog=u,sdlog=s) },
    d = function(x,m,het){ u = log(m/sqrt(1+het^2)); s = sqrt(log(1+het^2)); dlnorm(x,meanlog=u,sdlog=s) },
    p = function(q,m,het){ u = log(m/sqrt(1+het^2)); s = sqrt(log(1+het^2)); plnorm(q,meanlog=u,sdlog=s) },
    q = function(p,m,het){ u = log(m/sqrt(1+het^2)); s = sqrt(log(1+het^2)); qlnorm(p,meanlog=u,sdlog=s) })
  )

gen.makevars = function(opt=0){
  cat('\nCXX17=g++',
      '\nCXX17FLAGS=-march=native -mtune=native -fPIC -O',opt,
  sep='',file='~/.R/Makevars')
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

geom_summary = function(ps=NULL,alpha=.33){
  geom = lapply(ps,function(p){
    stat_summary(geom='ribbon',color=NA,alpha=alpha,
      fun.min=function(x){ quantile(x,  (1-p)/2) },
      fun.max=function(x){ quantile(x,1-(1-p)/2) })
  })
  geom = c(geom,stat_summary(geom='line',fun='median'))
}

scale_clr_viridis = function(...,end=.85){ list(
  viridis::scale_color_viridis(...,end=end),
  viridis::scale_fill_viridis(...,end=end)
)}

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
