library('reshape2')
library('LaplacesDemon')
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

fit.n = function(p,p.025,p.975){
  err.fun = function(x){
    n = 10^x
    q2 = qbeta(c(.025,.975),n*p,n*(1-p))
    err = sum((q2-c(p.025,p.975))^2)
  }
  n = 10^optimize(err.fun,c(0,4))$minimum
}

fit.mcv = function(fam,m,cv,b=+Inf,tol=1e-6){
  args = list(spec=fam,a=0,b=b)
  jfun = function(x){
    mx = do.call( extrunc,c(exp(x),args))
    vx = do.call(vartrunc,c(exp(x),args))
    err = (m-mx)^2 + (cv-sqrt(vx)/mx)^2 }
  x0 = log(switch(fam,
    gamma   = c(shape=1/cv^2,rate=m*cv^2),
    weibull = c(shape=1/cv^2,scale=1/m/cv^2),
    lnorm   = c(meanlog=log(m/sqrt(1+cv^2)),sdlog=sqrt(log(1+cv^2)))
  ))
  if (is.finite(b)|fam=='weibull'){
    opt = optim(x0,jfun)
    if (opt$value > tol){ cat('WARNING: fit.mcv\n'); print(opt) }
    return(c(exp(opt$par),args))
  } else {
    return(c(x0,args))
  }
}

mcv.fun = lapply(list(r=rtrunc,d=dtrunc,p=ptrunc,q=qtrunc),
  function(fun){ function(...){ args = fit.mcv(...)
    function(...){ do.call(fun,c(args,list(...)))
}}})

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

geom_viola = function(...){
  geom_violin(...,alpha=.5,scale='width')
}

geom_estimate = function(...,shape=18,width=1){
  suppressWarnings(list(
    geom_point(...,size=1,shape=shape),
    geom_errorbar(...,lwd=1/3,width=width)
))}

scale_clr_manual = function(...){ list(
  # TODO: combine w viridis
  scale_color_manual(...),
  scale_fill_manual(...)
)}

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
