#setwd("/home/zuo_r/involution")

#source("replicator_dynamics3.R")

replicator_dynamics_xy<-function(y,x,M,beta1,beta2,N,d,c,l){
  
  #-----合理预设参数解释说明------#
  
  # y 整体内卷的比例
  # x 整体躺平的比例
  
  # #内卷相对躺平的效用
  # beta1<-2
  # #挣扎相对躺平的效用
  # beta2<-1.1
  
  # #involution的成本 (defect)
  # d<-4
  # #sit-up的成本
  # c<-1
  # #lay down的成本
  # l<-0.5
  #-----------------------#
  
  Pc<-0
  Pd<-0
  Pl<-0 
  for (Nd in 0:(N-1)){ 
    for (Nl in 0:(N-1-Nd)){
      Nc <- N-1-Nd-Nl
      pai_d <- beta1*d*M/(beta1*(Nd+1)*d +beta2*Nc*c     +Nl*l)-d 
      pai_c <- beta2*c*M/(beta1*Nd*d     +beta2*(Nc+1)*c +Nl*l)-c          
      pai_l <-       l*M/(beta1*Nd*d     +beta2*Nc*c     +(Nl+1)*l)-l  
      Pd <- Pd + choose(N-1,Nd)*choose(N-1-Nd,Nl)*y^Nd*(1-x-y)^Nc*x^Nl*pai_d
      Pc <- Pc + choose(N-1,Nd)*choose(N-1-Nd,Nl)*y^Nd*(1-x-y)^Nc*x^Nl*pai_c
      Pl <- Pl + choose(N-1,Nd)*choose(N-1-Nd,Nl)*y^Nd*(1-x-y)^Nc*x^Nl*pai_l
    }
  }
  
  R_ <- y*Pd+x*Pl+(1-x-y)*Pc  ##均值
  y. <- y*(Pd-R_)/R_
  x. <- x*(Pl-R_)/R_
  s.<-(1-x-y)*(Pc-R_)/R_
  
  result1<-data.frame(y=y,x=x,M=M,beta1=beta1,beta2=beta2,d=d,N=N,c=c,l=l,y.=y.,x.=x.,s.=s.,R_=R_)#<---补充***
  return(result1)
}

#实验1里，合理预设参数即可
#beta2=1.1
#beta1=2
#d=4
#c=1
#l=0.5











