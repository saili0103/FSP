###Table 1(right)###
source('~/FSP/FSP-main-functions.R')
#library(Rfast)
f.star<-function(theta.star,x){
  theta.star[1]*abs(x[,1])+theta.star[1]*(abs(0.3+x[,2]))^theta.star[2]#+
  # theta.star[1]*(abs(x[,1])*abs(x[,2]))^theta.star[2]#+
  # theta.star[1]*(abs(x[,1])*abs(0.5-x[,2]))/2
}
#plot(x.test[,1],f.star(theta.star,x.test))
#plot(x.test[,2],f.star(theta.star,x.test))
f.ptr<-function(theta.star,x){
  #theta.star[1]*(abs(0.4+x[,2]))^theta.star[2]
  #(theta.star[1]*abs(x[,1])+theta.star[1]*(abs(0.4+x[,2]))^theta.star[2])*0.8
  rnorm(nrow(x),0,1)
}
####generate pre-trained data
library(np)
set.seed(234)
d=2
theta.star<-c(1,0.5) # true smoothness parameters
n.te=500 #number of test samples
n=300 #number of retrieved samples
Niter=300
range.te<-c(-0.5,0.5)
options(warn = -1)
N.ptr=1000
result.mat<-NULL
for(n in c(200,300,400,500)){
  for(iter in 1: Niter){
    #generate ptr-trained samples
    x.ptr<-matrix(runif(d*N.ptr,-1,1),ncol=d)
    y.ptr<-f.ptr(theta.star,x.ptr)+rnorm(N.ptr)#pre-trained model
    dat.ptr<-cbind(x.ptr,y.ptr)
    #generate test samples
    x.test<-matrix(runif(n.te*d,range.te[1],range.te[2]),ncol=d)
    f.test<-f.star(theta.star,x.test)
    y.test<-f.test+rnorm(n.te)
    dat.test<-cbind(x.test,y.test)
    bw.ptr<-npregbw(xdat=dat.ptr[,-3],ydat=dat.ptr[,3],itmax=10^3,ckertype='uniform')$bw
    bw.ptr<-pmin(bw.ptr, rep(2,2))
    ptr.te<-kernel.est(h=bw.ptr, dat.tr=dat.ptr, x.test=x.test)
    #mean((f.test-ptr.te)^2)
    ####FSP method####
    #####Step 1: sample retrieval####  
    n.sig<-round(n/4)
    x.sig<-matrix(runif(n.sig*d,range.te[1],range.te[2]),ncol=d)#target region:
    f.sig<-f.star(theta.star,x.sig)
    y.sig<-f.sig+rnorm(n.sig)
    h.sig<-min(sd(y.sig)*n^(-1/(d+2)),(range.te[2]-range.te[1])/4)
    ##Step 1.2: generate retrieved samples
    sample.re<-weighted.sample(x.sig,y.sig,h.sig,range.te,n)
    x.tr<-sample.re$x
    f.true<-f.star(theta.star, x.tr)
    y.tr<-f.true+rnorm(n-n.sig)
    dat.tr0<-cbind(x.tr,y.tr)
    dat.sig<-cbind(x.sig,y.sig)
    dat.tr<-rbind(dat.tr0, dat.sig)
    colnames(dat.tr)<-NULL
    y.ind<-ncol(dat.tr)
    ptr.tr0<-kernel.est(h=bw.ptr,dat.tr=dat.ptr, x.test=dat.tr0[,-y.ind])#ptr-trained estimate for the retrieved data
    ptr.sig<-kernel.est(h=bw.ptr,dat.tr=dat.ptr, x.test=dat.sig[,-y.ind])
    #####Step 2: Estimation####
    theta.mat<-expand.grid(seq(2, 0, length.out=4),seq(0.1,1,length.out=4))
    bw.st <- npregbw(xdat=dat.tr[,-3],ydat=dat.tr[,3], itmax=10^3,ckertype='uniform')$bw
    bw.st<-pmin(bw.st,rep(range.te[2]-range.te[1],2))
    st.re <- kernel.est(h=bw.st, dat.tr, x.test) #single-task method
    h.mat<-cbind(bw.st[1]*(4:1),bw.st[2]*(4:1))
    cv.re<-cv.fsp(dd.tr=dat.tr0, dd.val=dat.sig, ptr.val=ptr.sig, ptr.tr=ptr.tr0,
                  theta.mat=theta.mat, h.mat=h.mat)
    theta.hat<-as.numeric(cv.re[which.min(cv.re[,5]),1:2])
    h.hat<-as.numeric(cv.re[which.min(cv.re[,5]),3:4])
    fsp.re<-fsp.est(theta=theta.hat, ptr.tr=c(ptr.tr0,ptr.sig),
                    dat.tr=dat.tr, x.test=x.test, ptr.te=ptr.te, h=h.hat)
    
    ####section 2: Using randomly retrieved samples
    x.tr2<-matrix(runif(n*d,range.te[1],range.te[2]),ncol=d)
    ftr2.true<-f.star(theta.star, x.tr2)
    y.tr2<-ftr2.true+rnorm(n)
    dat.tr2<-cbind(x.tr2, y.tr2)
    colnames(dat.tr2)<-NULL
    ptr.tr2<-kernel.est(h=bw.ptr, dat.tr=dat.ptr, x.test=dat.tr2[,-y.ind])#ptr-trained estimate for the retrieved data
    ind.val<-sample(1:n, n/4)#validation-samples
    bw.st2 <- npregbw(xdat=dat.tr2[,-3],ydat=dat.tr2[,3], itmax=10^3,ckertype='uniform')$bw
    bw.st2<-pmin(bw.st2,rep(range.te[2]-range.te[1],2))
    st.re2 <- kernel.est(h=bw.st2, dat.tr2, x.test) #single-task method
    h.mat2<-cbind(bw.st2[1]*(4:1),bw.st2[2]*(4:1))
    cv.re2<-cv.fsp(dd.tr=dat.tr2[-ind.val,], dd.val=dat.tr2[ind.val,], ptr.val=ptr.tr2[ind.val], ptr.tr=ptr.tr2[-ind.val],
                   theta.mat=theta.mat,h.mat=h.mat2)
    theta.hat2<-as.numeric(cv.re2[which.min(cv.re2[,5]),1:2])
    h.hat2<-as.numeric(cv.re2[which.min(cv.re2[,5]),3:4])
    fsp.re2<-fsp.est(theta=theta.hat2, ptr.tr=ptr.tr2,
                     dat.tr=dat.tr2,x.test=x.test,ptr.te=ptr.te, h=h.hat2)
    result.mat<-rbind(result.mat,
                      c(n, theta.hat, mean((f.test-st.re)^2),
                        mean((f.test-fsp.re$fsp.re)^2), mean((f.test-ptr.te)^2),
                        mean((f.test-st.re2)^2), mean((f.test-fsp.re2$fsp.re)^2)))
  }
  
}

write.table(result.mat,file='~/FSP/FSP-sim-noise-n.txt',row.names=F)
#colMeans(result.mat)

