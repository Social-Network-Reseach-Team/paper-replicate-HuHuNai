suppressPackageStartupMessages(suppressWarnings({
  library(tidyverse)
  library(data.table)
  library(progress)
  library(parallel)
  library(foreach)
  library(doParallel)
  library(tidyverse)
  library(data.table)
  library(EvolutionaryGames)
  library(ggthemes)
  library(latex2exp)
  library(ggtext)
  library(tidyverse)
  library(data.table)
  library(EvolutionaryGames)
  library(ggthemes)
  library(latex2exp)
  library(ggtext)
  library(progress)
  library(plotly)
  library(RColorBrewer)
}))

setwd("F:/Code/R/10.1")
source("replicator_dynamics3R1.R")

# 
# #----------------------------------------------#
# #a.经济下行（GDP增速下降），群体如何演化？
# #----------------------------------------------#
# 
# #a1.参数构建
# parameters<-data.frame()
# for (y in (1:50)/50) {
#   x=seq(0,1-y,0.02)
#   temp<-data.frame(y=rep(y,length(x)),x=x)
#   parameters<-rbind(parameters,temp)
# }
# 
# M<-data.frame(M=seq(50,500,50))
# 
# parameters1<-data.frame(y=rep(parameters$y,nrow(M)),
#                           x=rep(parameters$x,nrow(M)),
#                           M=rep(M$M,nrow(parameters)))
# 
# 
# #a2.数据收集函数
# #只需要变化x,y,M即可
# get_data_by_M<-function(i){
#   
#   y<-parameters1$y[i]
#   x<-parameters1$x[i]
#   M<-parameters1$M[i]
#   
#   temp<-replicator_dynamics_xy(y,x,M,beta1=2,beta2=1.1,N=50,d=4,c=1,l=0.5)
#   
#   gc()
#   
#   data<-cbind(i,temp)
# }
# 
# #a3.并行收集数据
# t1<-Sys.time()
# n_core<-detectCores(logical = F)
# system.time({
#   cl<- makeCluster(n_core)      
#   registerDoParallel(cl)       #进行进程注册
#   clusterEvalQ(cl,{
#     N<-50
#     result<-data.frame()
#   }) 
#   result <- foreach(
#     i = 1:12750,
#     .combine=rbind   #返回结果的整合
#   ) %dopar% get_data_by_M(i)
#   stopCluster(cl)
# })
# 
# t2<-Sys.time()
# t2-t1 
# 
# fwrite(result,"result-Exp1.csv",row.names = F)

#--------------a4.数据可视化-------------------#
#fig1:经济下行（GDP增速下降），群体如何演化？

#金昊之前三维可视化z改为M
get_star_by_M<-function(M){
  x=0.3
  y=0.3
  time_start = Sys.time()
  temp<-replicator_dynamics_xy(y,x,M,beta1=1.5,beta2=1.1,N=50,d=4,c=1,l=0.5)
  x_flag = FALSE
  y_flag = FALSE

  
  while(!(x_flag & y_flag)){
    if (x>1){
      x_flag = TRUE
      x=sign(x)
    }
    if (x<0){
      x_flag = TRUE
      x=0
    }
    if (abs(y)>=1){
      y_flag = TRUE
      y=sign(y)
    }
    if (y<0){
      y_flag = TRUE
      y=0
    }
    temp <- replicator_dynamics_xy(y,x,M,beta1=1.5,beta2=1.1,N=50,d=4,c=1,l=0.5)
    x. = temp$x.
    y. = temp$y.
    if (abs(x.)>0.001 || abs(y.)>0.001){
      x = x+sign(x.)*0.002
      y = y+sign(y.)*0.002
          # if (abs(x.)>0.2){
          #   x. = sign(x.)*0.001
          # }
          # x = x+x.
          # if (abs(y.)>0.2){
          #   y. = sign(y.)*0.001
          # }
          # y = y+y.
    }
    else{
      break
    }
    if (Sys.time()-time_start>10){
       cat(M,"time out\n")
    }
  }
  cat(M," done",x,y,"\n")
  return(data.frame(xstar=x,ystar=y,M=M))
}

get_star_by_M_alt<-function(M){
  x=0.3
  y=0.3
  time_start = Sys.time()
  temp<-replicator_dynamics_xy_alt(y,x,M,beta1=1.5,beta2=1.1,N=50,d=4,c=1,l=0.5)
  x_flag = FALSE
  y_flag = FALSE
  
  
  while(!(x_flag & y_flag)){
    if (x>1){
      x_flag = TRUE
      x=sign(x)
    }
    if (x<0){
      x_flag = TRUE
      x=0
    }
    if (abs(y)>=1){
      y_flag = TRUE
      y=sign(y)
    }
    if (y<0){
      y_flag = TRUE
      y=0
    }
    temp <- replicator_dynamics_xy_alt(y,x,M,beta1=1.5,beta2=1.1,N=50,d=4,c=1,l=0.5)
    x. = temp$x.
    y. = temp$y.
    if (abs(x.)>0.001 || abs(y.)>0.001){
      if (abs(x.)>0.2){
        x. = sign(x.)*0.001
      }
      x = x+x.
      if (abs(y.)>0.2){
        y. = sign(y.)*0.001
      }
      y = y+y.
    }
    else{
      break
    }
    if (Sys.time()-time_start>10){
      cat(M,"time out\n")
    }
  }
  cat(M," done",x,y,"\n")
  return(data.frame(xstar=x,ystar=y,M=M))
}


get_star_by_N<-function(M){
  x=0.3
  y=0.3
  M=100
  time_start = Sys.time()
  temp<-replicator_dynamics_xy(y,x,M=M,beta1=1.5,beta2=1.1,N=N,d=4,c=1,l=0.5)
  x_flag = FALSE
  y_flag = FALSE
  
  
  while(!(x_flag & y_flag)){
    if (x>1){
      x_flag = TRUE
      x=sign(x)
    }
    if (x<0){
      x_flag = TRUE
      x=0
    }
    if (abs(y)>=1){
      y_flag = TRUE
      y=sign(y)
    }
    if (y<0){
      y_flag = TRUE
      y=0
    }
    temp <- replicator_dynamics_xy(y,x,M=M,beta1=1.5,beta2=1.1,N=N,d=4,c=1,l=0.5)
    x. = temp$x.
    y. = temp$y.
    if (abs(x.)>0.001 || abs(y.)>0.001){
      if (abs(x.)>0.2){
        x. = sign(x.)*0.001
      }
      x = x+x.
      if (abs(y.)>0.2){
        y. = sign(y.)*0.001
      }
      y = y+y.
    }
    else{
      break
    }
    if (Sys.time()-time_start>10){
      cat(M,"time out\n")
    }
  }
  cat(M," done",x,y,"\n")
  return(data.frame(xstar=x,ystar=y,N=N))
}

get_star_by_N_alt<-function(M){
  x=0.3
  y=0.3
  M=100
  time_start = Sys.time()
  temp<-replicator_dynamics_xy_alt(y,x,M=M,beta1=1.5,beta2=1.1,N=N,d=4,c=1,l=0.5)
  x_flag = FALSE
  y_flag = FALSE
  
  
  while(!(x_flag & y_flag)){
    if (x>1){
      x_flag = TRUE
      x=sign(x)
    }
    if (x<0){
      x_flag = TRUE
      x=0
    }
    if (abs(y)>=1){
      y_flag = TRUE
      y=sign(y)
    }
    if (y<0){
      y_flag = TRUE
      y=0
    }
    temp <- replicator_dynamics_xy_alt(y,x,M=M,beta1=1.5,beta2=1.1,N=N,d=4,c=1,l=0.5)
    x. = temp$x.
    y. = temp$y.
    if (abs(x.)>0.001 || abs(y.)>0.001){
      if (abs(x.)>0.2){
        x. = sign(x.)*0.001
      }
      x = x+x.
      if (abs(y.)>0.2){
        y. = sign(y.)*0.001
      }
      y = y+y.
    }
    else{
      break
    }
    if (Sys.time()-time_start>10){
      cat(M,"time out\n")
    }
  }
  cat(M," done",x,y,"\n")
  return(data.frame(xstar=x,ystar=y,N=N))
}

# df = data.frame()
# range = seq(50,500,10)
# for (M in range){
#   df=rbind(df,get_star_by_M(M))
# }
# df
# 
# df=df %>% mutate(zstar=1-xstar-ystar)
# 
# fig <- plot_ly(df, x = ~M, y = ~xstar, name = 'x*', type = 'scatter', mode = 'lines') 
# fig <- fig %>% add_trace(y = ~ystar, name = 'y*', mode = 'lines') 
# fig <- fig %>% add_trace(y = ~zstar, name = 'z*', mode = 'lines')
# 
# fig
# mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique(df$labelll)))
# fig <- plot_ly(df,x = ~ystar, y=~xstar, z=~M,type = "scatter3d",size = 1)
# fig

# get_star_by_M(77)
# 
# replicator_dynamics_xy(0.228507,0.2911438,77,beta1=1.5,beta2=1.1,N=50,d=4,c=1,l=0.5)
# 
# get_star_by_M(100)



