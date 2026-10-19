l <- as.numeric(commandArgs(trailingOnly=TRUE)[1]) #needed to queue jobs on storm
# file.names <- list.files(path = getwd(), pattern = "\\.R$") # list vector of file names
# file.names <- file.names[-grep("^make",file.names)] # any non-function scripts should start with the word make
# file.names <- file.names[-which(file.names %in% c("Main.R"))] # remove current file and Main (if different)
# lapply(file.names, source) # source all functions needed
# load("guillemot.RData")
# save(list=ls(),file="gf.storm.RData")
if(l==1){
  load("gf.storm.RData")
  library(parallel)
  cores <- 89
  n <- cores*5
  set.seed(1633762624)
  seeds <- sample(1:.Machine$integer.max,n)
  
  Time <- ncol(ch)-1
  states <- c("N","B1","LB","L_B","LB_","L_B_","S")
  struc <- list("phi"=list("age"=list(1:2,3:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "p"=list("age"=list(1:Time),"time"=list(1:6,7:Time),"state"=list(1:length(states))))
  theta <- logit(c(0.3,0.65,0.5,0.4,0.2,0.75,0.8,0.4))
  par.len <- length(theta)
  
  ni.s <- list()
  for(i in 1:length(ni2[1:15])){
    ni.s[[i]] <- rep("N",length(ni2[[i]]))
  }
  
  dat <- mclapply(seeds,function(x) dat.sim.wrap(theta,1:2,3,4,5,6,7:8,struc,ni2[1:15],ni.s,x),mc.cores=cores)
  mv <- mclapply(1:n,function(x) suff.stat(dat[[x]]),mc.cores=cores)
  
  op <- mclapply(1:n, function(x) optim(rep(0,par.len),
                                  ll.sl,
                                  phi.ind=1:2,
                                  delt.ind=3,
                                  kap.ind=4,
                                  rho.ind=5,
                                  gam.ind=6,
                                  p.ind=7:8,
                                  struc=struc,
                                  mv=mv[[x]],
                                  control=list(maxit=10000,fnscale=-1),
                                  hessian=TRUE),mc.cores=cores)
  
  hess.fun <- function(hess){
    if(det(hess)<1e-08){
      return(rep(NA,par.len))
    }else{
      return(diag(solve(-hess)))
    }
  }
  vars <- unlist(lapply(1:n,function(x) hess.fun(op[[x]]$hessian)))
  mles <- unlist(lapply(1:n,function(x) op[[x]]$par))
  cis <- ci(mles,vars)
  par.name <- factor(c("phi2","phi3","delta","kappa","rho","gamma","p7","p8"),levels=c("phi2","phi3","delta","kappa","rho","gamma","p7","p8"))
  df <- data.frame("par"=rep(par.name,n),
                   "mle"=unlist(lapply(1:n,function(x) logistic(op[[x]]$par))),
                   "convergence"=rep(sapply(1:n,function(x) op[[x]]$convergence),each=par.len),
                   "lower"=cis[,1],
                   "mle.un"=mles,
                   "upper"=cis[,3],
                   "var"=unlist(lapply(1:n,function(x) hess.fun(op[[x]]$hessian))),
                   "true"=rep(theta,n),
                   "n"=rep(1:n,each=par.len))
  df.true <- data.frame("par"=par.name,
                        "mle"=logistic(theta))
  save(df,df.true,file="sim1.RData")
  save(df,df.true,op,file="sim1.big.RData")
}

if(l==2){
  load("gf.storm.RData")
  library(parallel)
  cores <- 89
  n <- cores*5
  set.seed(1633762624)
  seeds <- sample(1:.Machine$integer.max,n)
  
  Time <- ncol(ch)-1
  states <- c("N","B1","LB","L_B","LB_","L_B_","S")
  struc <- list("phi"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "p"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))))
  theta <- logit(c(0.3,0.5,0.4,0.2,0.75,0.6))
  par.len <- length(theta)
  
  ni.s <- list()
  for(i in 1:length(ni2[1:15])){
    ni.s[[i]] <- rep("N",length(ni2[[i]]))
  }
  
  dat <- mclapply(seeds,function(x) dat.sim.wrap(theta,1,2,3,4,5,6,struc,ni2[1:15],ni.s,x),mc.cores=cores)
  mv <- mclapply(1:n,function(x) suff.stat(dat[[x]]),mc.cores=cores)
  
  op <- mclapply(1:n, function(x) optim(rep(0,par.len),
                                        ll.sl,
                                        phi.ind=1,
                                        delt.ind=2,
                                        kap.ind=3,
                                        rho.ind=4,
                                        gam.ind=5,
                                        p.ind=6,
                                        struc=struc,
                                        mv=mv[[x]],
                                        control=list(maxit=10000,fnscale=-1),
                                        hessian=TRUE),mc.cores=cores)
  
  hess.fun <- function(hess){
    if(det(hess)<1e-08){
      return(rep(NA,par.len))
    }else{
      return(diag(solve(-hess)))
    }
  }
  vars <- unlist(lapply(1:n,function(x) hess.fun(op[[x]]$hessian)))
  mles <- unlist(lapply(1:n,function(x) op[[x]]$par))
  cis <- ci(mles,vars)
  par.name <- factor(c("phi","delta","kappa","rho","gamma","p"),levels=c("phi","delta","kappa","rho","gamma","p"))
  df <- data.frame("par"=rep(par.name,n),
                   "mle"=unlist(lapply(1:n,function(x) logistic(op[[x]]$par))),
                   "convergence"=rep(sapply(1:n,function(x) op[[x]]$convergence),each=par.len),
                   "lower"=cis[,1],
                   "mle.un"=mles,
                   "upper"=cis[,3],
                   "var"=unlist(lapply(1:n,function(x) hess.fun(op[[x]]$hessian))),
                   "true"=rep(theta,n),
                   "n"=rep(1:n,each=par.len))
  df.true <- data.frame("par"=par.name,
                        "mle"=logistic(theta))
  save(df,df.true,file="sim2.RData")
  save(df,df.true,op,file="sim2.big.RData")
}

if(l==3){
  load("gf.storm.RData")
  library(parallel)
  cores <- 89
  n <- cores*5
  set.seed(1633762624)
  seeds <- sample(1:.Machine$integer.max,n)
  
  Time <- ncol(ch)-1
  states <- c("N","B1","LB","L_B","LB_","L_B_","S")
  struc <- list("phi"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "p"=list("age"=list(1:Time),"time"=list(1:6,7:Time),"state"=list(1,2:6,length(states))))
  theta <- logit(c(0.65,0.5,0.4,0.2,0.75,0.8,0.4))
  par.len <- length(theta)
  
  ni.s <- list()
  for(i in 1:length(ni2[1:15])){
    ni.s[[i]] <- rep("N",length(ni2[[i]]))
  }
  
  dat <- mclapply(seeds,function(x) dat.sim.wrap(theta,1,2,3,4,5,6:7,struc,ni2[1:15],ni.s,x),mc.cores=cores)
  mv <- mclapply(1:n,function(x) suff.stat(dat[[x]]),mc.cores=cores)
  
  op <- mclapply(1:n, function(x) optim(rep(0,par.len),
                                        ll.sl,
                                        phi.ind=1,
                                        delt.ind=2,
                                        kap.ind=3,
                                        rho.ind=4,
                                        gam.ind=5,
                                        p.ind=6:7,
                                        struc=struc,
                                        mv=mv[[x]],
                                        control=list(maxit=10000,fnscale=-1),
                                        hessian=TRUE),mc.cores=cores)
  
  hess.fun <- function(hess){
    if(det(hess)<1e-08){
      return(rep(NA,par.len))
    }else{
      return(diag(solve(-hess)))
    }
  }
  vars <- unlist(lapply(1:n,function(x) hess.fun(op[[x]]$hessian)))
  mles <- unlist(lapply(1:n,function(x) op[[x]]$par))
  cis <- ci(mles,vars)
  par.name <- factor(c("phi","delta","kappa","rho","gamma","p7","p8"),levels=c("phi","delta","kappa","rho","gamma","p7","p8"))
  df <- data.frame("par"=rep(par.name,n),
                   "mle"=unlist(lapply(1:n,function(x) logistic(op[[x]]$par))),
                   "convergence"=rep(sapply(1:n,function(x) op[[x]]$convergence),each=par.len),
                   "lower"=cis[,1],
                   "mle.un"=mles,
                   "upper"=cis[,3],
                   "var"=unlist(lapply(1:n,function(x) hess.fun(op[[x]]$hessian))),
                   "true"=rep(theta,n),
                   "n"=rep(1:n,each=par.len))
  df.true <- data.frame("par"=par.name,
                        "mle"=logistic(theta))
  save(df,df.true,file="sim3.RData")
  save(df,df.true,op,file="sim3.big.RData")
}

if(l==4){
  load("gf.storm.RData")
  library(parallel)
  cores <- 89
  n <- cores*5
  set.seed(1633762624)
  seeds <- sample(1:.Machine$integer.max,n)
  
  Time <- ncol(ch)-1
  states <- c("N","B1","LB","L_B","LB_","L_B_","S")
  struc <- list("phi"=list("age"=list(1,2:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "p"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))))
  theta <- logit(c(0.3,0.65,0.5,0.4,0.2,0.75,0.6))
  par.len <- length(theta)
  
  ni.s <- list()
  for(i in 1:length(ni2[1:15])){
    ni.s[[i]] <- rep("N",length(ni2[[i]]))
  }
  
  dat <- mclapply(seeds,function(x) dat.sim.wrap(theta,1:2,3,4,5,6,7,struc,ni2[1:15],ni.s,x),mc.cores=cores)
  mv <- mclapply(1:n,function(x) suff.stat(dat[[x]]),mc.cores=cores)
  
  op <- mclapply(1:n, function(x) optim(rep(0,par.len),
                                        ll.sl,
                                        phi.ind=1:2,
                                        delt.ind=3,
                                        kap.ind=4,
                                        rho.ind=5,
                                        gam.ind=6,
                                        p.ind=7,
                                        struc=struc,
                                        mv=mv[[x]],
                                        control=list(maxit=10000,fnscale=-1),
                                        hessian=TRUE),mc.cores=cores)
  
  hess.fun <- function(hess){
    if(det(hess)<1e-08){
      return(rep(NA,par.len))
    }else{
      return(diag(solve(-hess)))
    }
  }
  vars <- unlist(lapply(1:n,function(x) hess.fun(op[[x]]$hessian)))
  mles <- unlist(lapply(1:n,function(x) op[[x]]$par))
  cis <- ci(mles,vars)
  par.name <- factor(c("phi1","phi2","delta","kappa","rho","gamma","p"),levels=c("phi1","phi2","delta","kappa","rho","gamma","p"))
  df <- data.frame("par"=rep(par.name,n),
                   "mle"=unlist(lapply(1:n,function(x) logistic(op[[x]]$par))),
                   "convergence"=rep(sapply(1:n,function(x) op[[x]]$convergence),each=par.len),
                   "lower"=cis[,1],
                   "mle.un"=mles,
                   "upper"=cis[,3],
                   "var"=unlist(lapply(1:n,function(x) hess.fun(op[[x]]$hessian))),
                   "true"=rep(theta,n),
                   "n"=rep(1:n,each=par.len))
  df.true <- data.frame("par"=par.name,
                        "mle"=logistic(theta))
  save(df,df.true,file="sim4.RData")
  save(df,df.true,op,file="sim4.big.RData")
}

if(l==5){
  load("gf.storm.RData")
  library(parallel)
  cores <- 89
  n <- cores*5
  set.seed(1633762624)
  seeds <- sample(1:.Machine$integer.max,n)
  
  Time <- ncol(ch)-1
  states <- c("N","B1","LB","L_B","LB_","L_B_","S")
  struc <- list("phi"=list("age"=list(1:2,3:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "p"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))))
  theta <- logit(c(0.3,0.65,0.5,0.4,0.2,0.75,0.6))
  par.len <- length(theta)
  
  ni.s <- list()
  for(i in 1:length(ni2[1:15])){
    ni.s[[i]] <- rep("N",length(ni2[[i]]))
  }
  
  dat <- mclapply(seeds,function(x) dat.sim.wrap(theta,1:2,3,4,5,6,7,struc,ni2[1:15],ni.s,x),mc.cores=cores)
  mv <- mclapply(1:n,function(x) suff.stat(dat[[x]]),mc.cores=cores)
  
  op <- mclapply(1:n, function(x) optim(rep(0,par.len),
                                        ll.sl,
                                        phi.ind=1:2,
                                        delt.ind=3,
                                        kap.ind=4,
                                        rho.ind=5,
                                        gam.ind=6,
                                        p.ind=7,
                                        struc=struc,
                                        mv=mv[[x]],
                                        control=list(maxit=10000,fnscale=-1),
                                        hessian=TRUE),mc.cores=cores)
  
  hess.fun <- function(hess){
    if(det(hess)<1e-08){
      return(rep(NA,par.len))
    }else{
      return(diag(solve(-hess)))
    }
  }
  vars <- unlist(lapply(1:n,function(x) hess.fun(op[[x]]$hessian)))
  mles <- unlist(lapply(1:n,function(x) op[[x]]$par))
  cis <- ci(mles,vars)
  par.name <- factor(c("phi2","phi3","delta","kappa","rho","gamma","p"),levels=c("phi2","phi3","delta","kappa","rho","gamma","p"))
  df <- data.frame("par"=rep(par.name,n),
                   "mle"=unlist(lapply(1:n,function(x) logistic(op[[x]]$par))),
                   "convergence"=rep(sapply(1:n,function(x) op[[x]]$convergence),each=par.len),
                   "lower"=cis[,1],
                   "mle.un"=mles,
                   "upper"=cis[,3],
                   "var"=unlist(lapply(1:n,function(x) hess.fun(op[[x]]$hessian))),
                   "true"=rep(theta,n),
                   "n"=rep(1:n,each=par.len))
  df.true <- data.frame("par"=par.name,
                        "mle"=logistic(theta))
  save(df,df.true,file="sim5.RData")
  save(df,df.true,op,file="sim5.big.RData")
}
if(l==6){
  load("gf.storm.RData")
  library(parallel)
  cores <- 89
  n <- cores*5
  set.seed(1633762624)
  seeds <- sample(1:.Machine$integer.max,n)
  
  Time <- ncol(ch)-1
  states <- c("N","B1","LB","L_B","LB_","L_B_","S")
  
  theta <- logit(c(0.3,0.65,0.5,0.4,0.2,0.75,0.6))
  par.len <- length(theta)
  
  ni.s <- list()
  for(i in 1:length(ni2[1:15])){
    ni.s[[i]] <- rep("N",length(ni2[[i]]))
  }
  
  hess.fun <- function(hess){
    if(det(hess)<1e-08){
      return(rep(NA,par.len))
    }else{
      return(diag(solve(-hess)))
    }
  }
  
  phi.split <- 1:6
  phi.struc <- lapply(phi.split,function(x) list(1:x,(x+1):Time))
  par.name <- factor(c("phi1","phi2","delta","kappa","rho","gamma","p"),levels=c("phi1","phi2","delta","kappa","rho","gamma","p"))
  
  op <- lapply(phi.split,function(i){
    # par.name <- factor(c(paste0("phi",i:(i+1)),"delta","kappa","rho","gamma","p"),levels=c(paste0("phi",i:(i+1)),"delta","kappa","rho","gamma","p"))
    struc <- list("phi"=list("age"=phi.struc[[i]],"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "p"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))))
    dat <- mclapply(seeds,function(x) dat.sim.wrap(theta,1:2,3,4,5,6,7,struc,ni2[1:15],ni.s,x),mc.cores=cores)
    mv <- mclapply(1:n,function(x) suff.stat(dat[[x]]),mc.cores=cores)
  
    op <- mclapply(1:n, function(x) optim(rep(0,par.len),
                                          ll.sl,
                                          phi.ind=1:2,
                                          delt.ind=3,
                                          kap.ind=4,
                                          rho.ind=5,
                                          gam.ind=6,
                                          p.ind=7,
                                          struc=struc,
                                          mv=mv[[x]],
                                          control=list(maxit=10000,fnscale=-1),
                                          hessian=TRUE),mc.cores=cores)
    vars <- unlist(lapply(1:n,function(x) hess.fun(op[[x]]$hessian)))
    mles <- unlist(lapply(1:n,function(x) op[[x]]$par))
    cis <- ci(mles,vars)
    
    df <- data.frame("par"=rep(par.name,n),
                     "mle"=unlist(lapply(1:n,function(x) logistic(op[[x]]$par))),
                     "convergence"=rep(sapply(1:n,function(x) op[[x]]$convergence),each=par.len),
                     "lower"=cis[,1],
                     "mle.un"=mles,
                     "upper"=cis[,3],
                     "var"=unlist(lapply(1:n,function(x) hess.fun(op[[x]]$hessian))),
                     "true"=rep(theta,n),
                     "n"=rep(1:n,each=par.len),
                     "phi.split"=rep(i,par.len*n))
    df.true <- data.frame("par"=par.name,
                          "mle"=logistic(theta),
                          "phi.split"=rep(i,par.len))
    print(i)
    return(list("df"=df,"df.true"=df.true))
  })
  
  dfs <- lapply(phi.split,function(i) op[[i]]$df)
  df.trues <- lapply(phi.split,function(i) op[[i]]$df.true)
  
  df <- do.call(rbind.data.frame, dfs)
  df.true <- do.call(rbind.data.frame, df.trues)
  
  save(df,df.true,file="sim6.RData")
  save(df,df.true,op,file="sim6.big.RData")
}
if(l==7){
  load("gf.storm.RData")
  library(parallel)
  cores <- 89
  n <- cores*5
  set.seed(1633762624)
  seeds <- sample(1:.Machine$integer.max,n)
  
  Time <- ncol(ch)-1
  states <- c("N","B1","LB","L_B","LB_","L_B_","S")
  struc <- list("phi"=list("age"=list(1,2,3:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
                "p"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))))
  theta <- logit(c(0.3,0.55,0.7,0.5,0.4,0.2,0.75,0.6))
  par.len <- length(theta)
  
  ni.s <- list()
  for(i in 1:length(ni2[1:15])){
    ni.s[[i]] <- rep("N",length(ni2[[i]]))
  }
  
  dat <- mclapply(seeds,function(x) dat.sim.wrap(theta,1:3,4,5,6,7,8,struc,ni2[1:15],ni.s,x),mc.cores=cores)
  mv <- mclapply(1:n,function(x) suff.stat(dat[[x]]),mc.cores=cores)
  
  op <- mclapply(1:n, function(x) optim(rep(0,par.len),
                                        ll.sl,
                                        phi.ind=1:3,
                                        delt.ind=4,
                                        kap.ind=5,
                                        rho.ind=6,
                                        gam.ind=7,
                                        p.ind=8,
                                        struc=struc,
                                        mv=mv[[x]],
                                        control=list(maxit=10000,fnscale=-1),
                                        hessian=TRUE),mc.cores=cores)
  
  hess.fun <- function(hess){
    if(det(hess)<1e-08){
      return(rep(NA,par.len))
    }else{
      return(diag(solve(-hess)))
    }
  }
  vars <- unlist(lapply(1:n,function(x) hess.fun(op[[x]]$hessian)))
  mles <- unlist(lapply(1:n,function(x) op[[x]]$par))
  cis <- ci(mles,vars)
  par.name <- factor(c("phi1","phi2","phi3","delta","kappa","rho","gamma","p"),levels=c("phi1","phi2","phi3","delta","kappa","rho","gamma","p"))
  df <- data.frame("par"=rep(par.name,n),
                   "mle"=unlist(lapply(1:n,function(x) logistic(op[[x]]$par))),
                   "convergence"=rep(sapply(1:n,function(x) op[[x]]$convergence),each=par.len),
                   "lower"=cis[,1],
                   "mle.un"=mles,
                   "upper"=cis[,3],
                   "var"=unlist(lapply(1:n,function(x) hess.fun(op[[x]]$hessian))),
                   "true"=rep(theta,n),
                   "n"=rep(1:n,each=par.len))
  df.true <- data.frame("par"=par.name,
                        "mle"=logistic(theta))
  save(df,df.true,file="sim7.RData")
  save(df,df.true,op,file="sim7.big.RData")
}