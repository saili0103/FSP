#l_2 distance
l2<-function(a,b){
  if(is.matrix(a)){
    apply(a,1,function(xx) sqrt(sum((xx-b)^2)))
  }else{
    sqrt(sum((a-b)^2))
  }
}
#kernel estimate 
kernel.est<-function(h=NULL, dat.tr, x.test){
  d<-ncol(dat.tr)-1
  y.ind<-ncol(dat.tr)
  y.tr<-dat.tr[,y.ind]
  x.tr<-dat.tr[,-y.ind]
  if(is.null(h)){
    h=sd(dat.tr[,y.ind])*nrow(x.tr)^(-1/(ncol(x.tr)+2))
  }
  yhat<-rep(0,nrow(x.test))
  for(i in 1:nrow(x.test)){
    yhat[i] <- mean(y.tr[apply(x.tr,1, function(xx) max(abs(xx-x.test[i,])/h)<=1)])
  }
  yhat[is.na(yhat)]<-mean(y.tr)
  yhat
}

#compute delta(x)
delta.hat<-function(theta,x.star,dat.local,ptr.tr.local,ptr.star){
  y.ind<-ncol(dat.local)
  x.sel<-dat.local[,-y.ind]
  y.sel<-dat.local[,y.ind]
  yhat.ptr<-ptr.tr.local
  yhat.star<-ptr.star
  ftrans<-yhat.star+pmin(abs(yhat.ptr-yhat.star),theta[1]*(l2(x.sel,x.star))^theta[2])*sign(yhat.ptr-yhat.star)
  #ftrans<-yhat.ptr
  mean(y.sel-ftrans)
}

###sampling
weighted.sample<-function(x.sig,y.sig, h.sig,range.te, n){
  n.tr <- n - nrow(x.sig)
  x.tr <- matrix(0, nrow=n.tr,ncol=d)
  accepted <- 0
  total_trials = 0
  sig.vec <- rep(0,nrow(x.sig))
  for(i in 1:length(sig.vec)){
    dmax<-apply(x.sig,1, function(xx) max(abs(xx-x.sig[i,])))
    ind.local<-which(dmax<=h.sig)
    if(length(ind.local)<=2){
      sig.vec[i]<-var(y.sig)
    }else{
      sig.vec[i] <- sum(y.sig[ind.local]^2*(1-dmax[ind.local]/h.sig))/sum(1-dmax[ind.local]/h.sig)-
        (sum(y.sig[ind.local]*(1-dmax[ind.local]/h.sig))/sum(1-dmax[ind.local]/h.sig))^2
    }
  }
  #plot(x.sig[,1],sig.vec)
  #plot(x.sig[,2],sig.vec)
  M <- max(sig.vec)/mean(sig.vec)
  if(is.vector(range.te)){
    range.te<-matrix(rep(range.te,d),nrow=d,byrow=T)
  }
  while(accepted<n.tr){
    total_trials <- total_trials + 1
    x.candidate <- sapply(1:d,function(k) runif(1,range.te[k,1],range.te[k,2]))
    dmax<-apply(x.sig,1, function(xx) max(abs(xx-x.candidate)))
    ind.local<-which(dmax<=h.sig)
    if(length(ind.local)<=2){
      sig.cand=var(y.sig)
    }else{
      sig.cand <- sum(y.sig[ind.local]^2*(1-dmax[ind.local]/h.sig))/sum(1-dmax[ind.local]/h.sig)-
        (sum(y.sig[ind.local]*(1-dmax[ind.local]/h.sig))/sum(1-dmax[ind.local]/h.sig))^2
    }
    pX.candidate <- sig.cand/mean(sig.vec)
    u <- runif(1, 0, M)
    # Accept if u <= p_{X}(x_candidate)
    if (u <= pX.candidate) {
      accepted <- accepted + 1
      x.tr[accepted,] <- x.candidate
    }
  }
  acceptance_rate <- n.tr / total_trials
  #cat("Acceptance rate:", acceptance_rate, "\n")
  
  return(list(x=x.tr, sig.vec=sig.vec))
}

#FSP estimator given theta and h
fsp.est<-function(theta, h, dat.tr, x.test,ptr.te, ptr.tr){
  fsp.re<-rep(0,nrow(x.test))
  n<-nrow(dat.tr)
  y.ind<-ncol(dat.tr)
  x.tr<-dat.tr[,-y.ind]
  delta.vec<-rep(0,nrow(x.test))
  for(i in 1:nrow(x.test)){
    x.star<-x.test[i,]
    ind.cur<-which(apply(x.tr,1, function(xx) max(abs(xx-x.star)/h)<=1))
    if(length(ind.cur)==0){
      delta.vec[i]<-mean(dat.tr[,y.ind])-ptr.te[i]
    }else if(length(ind.cur)==1){
      delta.vec[i]<-delta.hat(theta, x.star=x.star, dat.local=matrix(dat.tr[ind.cur,],nrow=1),
                              ptr.tr.local=ptr.tr[ind.cur],ptr.star=ptr.te[i])
    }else{
      delta.vec[i]<-delta.hat(theta, x.star=x.star, dat.local=dat.tr[ind.cur,], 
                              ptr.tr.local=ptr.tr[ind.cur],ptr.star=ptr.te[i])
    }
  }
  
  
  return(list(fsp.re=ptr.te + delta.vec,delta.vec=delta.vec, h=h))
}

#cross-fitting to select tuning parameters
cv.fsp<-function(dd.tr, dd.val, ptr.val, ptr.tr, theta.mat, h.mat){
  y.ind<-ncol(dd.tr)
  y.val<-dd.val[,y.ind]
  err.re<-NULL
  for(k3 in 1:nrow(theta.mat)){#range over theta
    theta<-as.numeric(theta.mat[k3,])
    for(k2 in 1:nrow(h.mat)){
      if(theta.mat[k3,1]==0){
        hh<-bw.st
        if(k3 > 1){next}
      }else{
        hh<-h.mat[k2,]
      }
      ftheta.re<-fsp.est(theta= theta, h=hh, dat.tr=dd.tr,
                         x.test=dd.val[,-y.ind],ptr.te=ptr.val, ptr.tr=ptr.tr)
      fhat.fsp<-ftheta.re$fsp.re
      err.re<-rbind(err.re,
                    c( theta, hh, mean((fhat.fsp-y.val)^2)))
    }
    
  }
  
  err.re
}
