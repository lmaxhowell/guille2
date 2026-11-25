l <- as.numeric(commandArgs(trailingOnly=TRUE)[1]) #needed to queue jobs on storm

if(l==1){
  load("gf.storm.RData")
  library(parallel)
  
  Time <- ncol(ch)-1
  mv <- suff.stat(ch)
  states <- c("N","B1","LB","L_B","LB_","L_B_","S")
  
  rho.trans <- function(rho.input){
    rho <- array(0,dim=c(Time-1,Time,length(states)))
    # rho[,1,] <- 0 # sets as default anyway # age 0
    rho[,2:3,] <- rho.input[1] # age one and two
    rho[,4,] <- rho.input[2] # age 3
    rho[,5:Time,] <- rho.input[3] # age 4+
    return(rho)
  }
  
  phi.trans <- function(phi.input,age.split,sd){
    # phi input is vector of parameters to be distributed
    # age.split is the age at which the two age classes should be split if there are two age classes
    # sd is a vector to decide if there should be age dependence in N or B1 or both
    # it should be either 1,2,or 1:2 to correspond with this
    if(age.split==0){ # only state dependence
      if(length(phi.input)!=length(states)){
        print(paste("input length:",length(phi.input)))
        print(paste("correct length:",length(states)))
        stop("wrong number of parameters for phi")
      }
      phi <- untrans(phi.input,list(1:Time),list(1:(Time-1)),as.list(1:length(states)))
    }else if(age.split>0){ # state and age dependence
      if(length(sd)==1){
        if(length(phi.input)!=(length(states)+1)){
          print(paste("input length:",length(phi.input)))
          print(paste("correct length:",(length(states)+1)))
          stop("wrong number of parameters for phi")
        }
        if(sd==1){
          phi.vec <- unlist(lapply(1:length(phi.input),function(x) if(x %in% 1:2){phi.input[x]}else{ rep(phi.input[x],2)}))
        }
        else if(sd==2){
          phi.vec <- unlist(lapply(1:length(phi.input),function(x) if(x %in% 2:3){phi.input[x]}else{ rep(phi.input[x],2)}))
          # print(phi.vec)
        }
      }else if(length(sd)==2){
        if(length(phi.input)!=(length(states)+2)){
          print(paste("input length:",length(phi.input)))
          print(paste("correct length:",(length(states)+2)))
          stop("wrong number of parameters for phi")
        }
        phi.vec <- unlist(lapply(1:length(phi.input),function(x) if(x %in% 1:4){
                                                                    phi.input[x]
                                                                  }else{ 
                                                                    rep(phi.input[x],2)
                                                                  }
                                                                ))
      }
      phi <- untrans(phi.vec,list(1:age.split,(age.split+1):Time),list(1:(Time-1)),as.list(1:length(states)))
    }
    return(phi)
  }
  
  model.retrieve.phi <- function(num.s,num.a,theta){
    if(num.s=="c" && num.a=="c"){
      phi <- untrans(logistic(theta[1]),list(1:Time),list(1:(Time-1)),list(1:length(states)))
      theta.cut <- 1
    }else if(num.a %in% as.character(0:(Time-1))){
      age.split <- as.numeric(num.a)
      sd <- as.numeric(strsplit(num.s,"")[[1]]) # so if num.s is 1 or 2 is returns that vector, and if it is 12 is becomes c(1,2)
      if(age.split==0){
        theta.cut <- length(states)
      }else{
        theta.cut <- length(states) + length(sd)
      }
      
      phi <- phi.trans(logistic(theta[1:theta.cut]),age.split,sd)
      
    }
    return(list("phi"=phi,"theta.cut"=theta.cut))
  }
  
  model.retrieve.other <- function(num,removal.factor){
    # num is the code for a given parameter
    # removal.factor is whether theres restrictions on state dependence for that parameter
    if(num %in% c("c")){
      struc <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
      len <- 1
    }else if(num %in% c("d")){
      struc <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=as.list(1:length(states)))
      len <- length(states)-removal.factor
    # }else if(is.numeric(num)){
    #   struc <- list(1:num,(num+1):Time)
    #   len <- length(2)
    }
    return(list("struc"=struc,"len"=len))
  }
  
  ll.m <- function(theta,mv,i){
    # theta is the vector of parameters
    # mv is the data
    # i is a number corresponding to the row of the model df
    # rho is always the same
    
    # phi
    mcs <- df.models[i,]
    num.s <- mcs[3]
    num.a <- mcs[4]
    mrp <- model.retrieve.phi(num.s,num.a,theta) # model retrieve phi
    phi <- mrp$phi
    theta <- theta[-(1:(mrp$theta.cut))]
    # View(phi[,,1])
    
    # delta
    mro <- model.retrieve.other(df.models[i,]$delta,1) # model retrieve other
    if(mro$len==1){
      theta.del <- theta[mro$len]
    }else{
      theta.del <- c(-Inf,theta[mro$len])
    }
    del <- untrans(logistic(theta.del),mro$struc$age,mro$struc$time,mro$struc$state)
    theta <- theta[-mro$len]
    # View(del[,,1])
    
    # kappa
    mro <- model.retrieve.other(df.models[i,]$kappa,2) # model retrieve other
    if(mro$len==1){
      theta.kap <- theta[mro$len]
    }else{
      theta.kap <- c(-Inf,theta[mro$len],-Inf)
    }
    kap <- untrans(logistic(theta.kap),mro$struc$age,mro$struc$time,mro$struc$state)
    theta <- theta[-mro$len]
    # View(kap[,,1])
    
    
    # rho
    rho <- rho.trans(logistic(theta[1:3]))
    theta <- theta[-3]
    # View(rho[1,,])
    
    # gamma
    mro <- model.retrieve.other(df.models[i,]$gamma,1) # model retrieve other
    if(mro$len==1){
      theta.gam <- theta[mro$len]
    }else{
      theta.gam <- c(-Inf,theta[mro$len])
    }
    gam <- untrans(logistic(theta.gam),mro$struc$age,mro$struc$time,mro$struc$state)
    theta <- theta[-mro$len]
    # View(gam[,,1])
    
    
    # p
    p  <- untrans(logistic(c(theta[1],Inf,-Inf)),list(1:Time),list(1:(Time)),list(1,2:6,length(states)))
    # View(p[,,1])
    
    psi <- make.psi(del,kap,rho,gam)
    
    ll <- sl(phi,psi,p,mv)
    return(ll)
  }
  
  theta.length <- function(i){
    if(df.models$phi.state[i]=="c" && df.models$phi.age[i]=="c"){
      theta.phi.len <- 1
    }else{
      theta.phi.len <- length(states)*length(strsplit(df.models$phi.state[i],"")[[1]])
    }
    if(sum(df.models[i,4:7]=="c")==4){
      theta.other.len <- 4 # four constant parameters
    }else if(sum(df.models[i,4:7]=="c")==3){
      theta.other.len <- length(states) + 3 # three constant parameters and one state dependant
    }
    len <- theta.phi.len + theta.other.len + 3 # plus three for rho
    return(len)
  }
  
  replen <- 3*length(0:7)+1
  df.models <- data.frame("model.code"=1,
                          "phi.state"=rep(c("c",rep("1",length(0:7)),rep("2",length(0:7)),rep("12",length(0:7))),5), # repeat five times for all parameters to be state dependent one at a time
                          "phi.age"=rep(c("c",rep(0:7,3)),5),
                          "delta"=c(rep("c",1*replen),rep("d",1*replen),rep("c",3*replen)),
                          "kappa"=c(rep("c",2*replen),rep("d",1*replen),rep("c",2*replen)),
                          "gamma"=c(rep("c",3*replen),rep("d",1*replen),rep("c",1*replen)),
                          "p"=c(rep("c",4*replen),rep("d",1*replen)))
  df.models$model.code <- 1:nrow(df.models)
  
  tls <- sapply(1:nrow(df.models),function(i) theta.length(i))
  
  set.seed(846541147)
  seeds <- sample(1:.Machine$integer.max,10)
  cores <- 89
  op <- mclapply(1:nrow(df.models),function(i){
    lapply(seeds,function(x){
      set.seed(x)
      starter <- logit(runif(tls[i]))
      opp <- optim(starter,
                   ll.m,
                   mv=mv,
                   i=i,
                   control=list(maxit=10000,fnscale=-1),
                   # hessian=TRUE
                   )
      return(c(opp,list("start"=starter)))
    })
  },mc.cores=cores)
  
  selector <- function(i){
    reps <- op[[i]]
    vals <- sapply(1:length(seeds),function(x) reps[[x]]$value)
    mins <- which.min(vals)
    
    return(reps[[mins]])
  }
  
  selected <- lapply(1:nrow(df.models),selector)
  
  aic.sel <- sapply(1:nrow(df.models),function(x) 2*tls[x] - 2*selected[[x]]$value)
  con.sel <- sapply(1:nrow(df.models),function(x) selected[[x]]$convergence)
  
  par.name.fun <- function(i){
    if(df.models$phi.state[i]=="c" && df.models$phi.age[i]=="c"){
      phi.name <- c("phi")
    }else if(df.models$phi.state[i]!="c"){
      sd <- strsplit(df.models$phi.state[i],"")[[1]]
      phi.name <- c()
      for(x in 1:length(sd)){
        phi.name <- c(phi.name,paste0("phi",sd[x],states))
      }
    }
    otherparams <- c("delta","kappa","gamma","p")
    for(x in 4:7){
      if(df.models[i,x]=="c"){
        assign(paste0(otherparams[x-3],".name"),otherparams[x-3])
      }else if(df.models[i,x]=="d"){
        assign(paste0(otherparams[x-3],".name"),paste0(otherparams[x-3],states))
      }
    }
    rho.name <- paste0("rho",c(23,4,5))
    
    return(c(phi.name,delta.name,kappa.name,rho.name,gamma.name,p.name))
  }
  
  par.names <- lapply(1:nrow(df.models),function(i) par.name.fun(i) )
  
  df.mle <- data.frame("par"=unlist(par.names),
                       "mle"=unlist(lapply(1:nrow(df.models),function(i) logistic(selected[[i]]$par))),
                       "aic"=rep(aic.sel,tls),
                       "convergence"=rep(unlist(lapply(1:nrow(df.models),function(i) selected[[i]]$convergence)),tls),
                       "model.code"=rep(df.models$model.code,tls))
  
  df.full <- data.frame("par"=unlist(sapply(par.names,rep,length(seeds))),
                        "mle"=unlist(lapply(1,function(i) unlist(lapply(1:length(seeds),function(x) logistic(op[[i]][[x]]$par))))),
                        "aic"=rep(aic.sel,tls),
                        "convergence"=rep(unlist(lapply(1,function(i) unlist(lapply(1:length(seeds),function(x) op[[i]][[x]]$convergence)))),tls),
                        "model.code"=rep(df.models$model.code,tls*length(seeds)))
  
  df.ais <- data.frame("model.code"=1:nrow(df.models),
                       "aic"=aic.sel,
                       "convergence"=con.sel)
  
  save(df.mle,df.full,df.ais,df.models,file="ms.RData")
  save(df.mle,df.full,df.ais,df.models,op,file="ms.big.RData")
}