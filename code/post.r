source('meta.r')

l$reg = 'Region'
fl$data = list(
  za='Observed durations:\ncensored among sampled (z|s)',
  x ='Target durations:\ntotal overall (x)')
cmap$reg = c(Africa='#cc3300',Asia='#ff9900',Europe='#003399',LatAm='#cc00cc',NorAm='#009933')

plot.distr = function(sub=1,dd=.1,zoom=15){
  types = list(S=' Survival function: S( • ) ',m='Mean: E[ • ]')
  Y = subset(prop.ci(load.csv('data','Fazito2012x')),kp=='fsw' & K26==1)
  S = subset(load.rds('data','stan','meta'),resam)
  d = d.vec(dd)
  D = rbind.lapply(seq(1,nrow(S),sub),function(g){
    fam = S$fam[g]
    px = sum1(p.x (d,fam=fam,a=S$a[g],b=S$b[g])) / dd # TODO: why need sum1?
    pz = sum1(p.za(d,fam=fam,a=S$a[g],b=S$b[g])
            * p.sz(d,sp0=S$sp0[g],st95=S$st95[g]) ) / dd
    Di = cbind(region=S$region[g],rbind(
      df(d=d,data=fl$data$x, type='S',value=1-cumsum(px)*dd),
      df(d=d,data=fl$data$za,type='S',value=1-cumsum(pz)*dd),
      df(d=0,data=fl$data$x, type='m',value=xp.mean(d,px)),
      df(d=0,data=fl$data$za,type='m',value=xp.mean(d,pz))))
  })
  geom_data = function(m,type,data,...){ geom_point(
    data=cbind(subset(Y,meas==m),type=type,data=data),shape=16,alpha=1/5,...) }
  D$type = factor(D$type,names(types),types)
  E = subset(D,type==types$m)
  D = subset(D,type!=types$m & d<=zoom)
  g = ggplot(D,aes(x=d,y=value,color=region,fill=region)) +
    facet_grid('type ~ data',scales='free') +
    geom_data('q',   types$S,fl$data$za,map=aes(x=value,y=1-p,size=n.adj(n))) +
    geom_data('mean',types$m,fl$data$za,map=aes(y=region,x=value,size=n.adj(n))) +
    geom_viola(data=E,aes(x=value,y=region),show.legend=0) +
    geom_mcrib() +
    scale_size_area(limits=c(10,1000),breaks=c(30,100,300,1000)) +
    scale_colorfill(v=setNames(cmap$reg,names(fl$region))) +
    scale_x_continuous(lim=c(0,zoom),breaks=seq(0,zoom,2)) +
    ggh4x::scale_y_facet(type==types$S,lim=0:1) +
    ggh4x::scale_y_facet(type==types$m,type='discrete') +
    ggh4x::force_panelsizes(rows=c(3,2)) +
    labs(x=l$dur,y=NULL,color=l$reg,fill=l$reg,size='Data N\'')
  plot.save(g,'stan','meta.distrs.p',size=c(6,3))
}

plot.distr()
