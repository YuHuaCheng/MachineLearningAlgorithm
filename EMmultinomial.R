EM.multinomial<-function(data,K,tau){
  n<-nrow(data)
  data=data+0.01
  #set some initial parameters
  delta=1000
  a.previous<-matrix(-1,nrow=n,ncol=K)
  #set initial value of t and c
  c<-rep(1/K,K)
  t.mat<-data[sample(1:n,K),]
  t<-t.mat/rowSums(t.mat)
  while(delta>tau){
  #E-step
  phi<-exp(t(log(t)%*%t(data)))
  a.temp<-diag(c)%*%t(phi)
  a<-t(a.temp)/colSums(a.temp)
  a[is.nan(a)]<-0
  delta<-norm(a-a.previous,"o")
  #M-step
  c<-colSums(a)/n
  bk<-t(a)%*%data
  t<-bk/rowSums(bk)
  a.previous<-a
  }
  m<-apply(a,1,which.max)
  m
  }
