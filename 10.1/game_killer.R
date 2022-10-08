get_star_by_M = function(M){
  helper = function(xy){
    M<-M
    beta1<-1.5
    beta2<-1.1
    N<-50
    d<-4
    c<-1
    l<-0.5
    x = xy[1]
    y = xy[2]
    l=0.5
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
    y. <- y*(Pd-R_)
    x. <- x*(Pl-R_)
    s.<-(1-x-y)*(Pc-R_)
    return (c(x.,y.))
  }
  

  out = dfsane(c(0.3,0.3),helper)$par
  return(data.frame(xstar=out[1],ystar=out[2],M=M))
}

get_star_by_M_with_beta = function(M){
  helper_with_beta = function(xy){
    beta1<-1.5
    beta2<-1.1
    M<-M
    N<-50
    d<-4
    c<-1
    l<-0.5
    x = xy[1]
    y = xy[2]
    l=0.5
    Pc<-0
    Pd<-0
    Pl<-0 
    for (Nd in 0:(N-1)){ 
      for (Nl in 0:(N-1-Nd)){
        Nc <- N-1-Nd-Nl
        B1 = beta1 - (beta1-1)/(1+exp((-10)*(Nc/N-0.5)))
        B2 = beta2 - (beta2-beta1)/(1+exp((-10)*(Nd/N-0.5)))
        pai_d <- B1*d*M/(B1*(Nd+1)*d +B2*Nc*c     +Nl*l)-d 
        pai_c <- B2*c*M/(B1*Nd*d     +B2*(Nc+1)*c +Nl*l)-c          
        pai_l <-       l*M/(B1*Nd*d     +B2*Nc*c     +(Nl+1)*l)-l  
        Pd <- Pd + choose(N-1,Nd)*choose(N-1-Nd,Nl)*y^Nd*(1-x-y)^Nc*x^Nl*pai_d
        Pc <- Pc + choose(N-1,Nd)*choose(N-1-Nd,Nl)*y^Nd*(1-x-y)^Nc*x^Nl*pai_c
        Pl <- Pl + choose(N-1,Nd)*choose(N-1-Nd,Nl)*y^Nd*(1-x-y)^Nc*x^Nl*pai_l
      }
    }
    
    R_ <- y*Pd+x*Pl+(1-x-y)*Pc  ##均值
    y. <- y*(Pd-R_)
    x. <- x*(Pl-R_)
    s.<-(1-x-y)*(Pc-R_)
    return (c(x.,y.))
  }

  out = dfsane(c(0.3,0.3),helper_with_beta)$par
  return(data.frame(xstar=out[1],ystar=out[2],M=M))
}

get_star_by_N = function(N){
  helper = function(xy){
    M<-100
    beta1<-1.5
    beta2<-1.1
    N<-N
    d<-4
    c<-1
    l<-0.5
    x = xy[1]
    y = xy[2]
    l=0.5
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
    y. <- y*(Pd-R_)
    x. <- x*(Pl-R_)
    s.<-(1-x-y)*(Pc-R_)
    return (c(x.,y.))
  }
  out = dfsane(c(0.3,0.3),helper)$par
  return(data.frame(xstar=out[1],ystar=out[2],N=N))
}

get_star_by_M_with_beta = function(N){
  helper_with_beta = function(xy){
    beta1<-1.5
    beta2<-1.1
    M<-100
    N<-N
    d<-4
    c<-1
    l<-0.5
    x = xy[1]
    y = xy[2]
    l=0.5
    Pc<-0
    Pd<-0
    Pl<-0 
    for (Nd in 0:(N-1)){ 
      for (Nl in 0:(N-1-Nd)){
        Nc <- N-1-Nd-Nl
        B1 = beta1 - (beta1-1)/(1+exp((-10)*(Nc/N-0.5)))
        B2 = beta2 - (beta2-beta1)/(1+exp((-10)*(Nd/N-0.5)))
        pai_d <- B1*d*M/(B1*(Nd+1)*d +B2*Nc*c     +Nl*l)-d 
        pai_c <- B2*c*M/(B1*Nd*d     +B2*(Nc+1)*c +Nl*l)-c          
        pai_l <-       l*M/(B1*Nd*d     +B2*Nc*c     +(Nl+1)*l)-l  
        Pd <- Pd + choose(N-1,Nd)*choose(N-1-Nd,Nl)*y^Nd*(1-x-y)^Nc*x^Nl*pai_d
        Pc <- Pc + choose(N-1,Nd)*choose(N-1-Nd,Nl)*y^Nd*(1-x-y)^Nc*x^Nl*pai_c
        Pl <- Pl + choose(N-1,Nd)*choose(N-1-Nd,Nl)*y^Nd*(1-x-y)^Nc*x^Nl*pai_l
      }
    }
    
    R_ <- y*Pd+x*Pl+(1-x-y)*Pc  ##均值
    y. <- y*(Pd-R_)
    x. <- x*(Pl-R_)
    s.<-(1-x-y)*(Pc-R_)
    return (c(x.,y.))
  }
  out = dfsane(c(0.3,0.3),helper_with_beta)$par
  return(data.frame(xstar=out[1],ystar=out[2],N=N))
}


