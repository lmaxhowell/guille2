file.names <- list.files(path = getwd(), pattern = "\\.R$") # list vector of file names
file.names <- file.names[-grep("^make",file.names)] # any non-function scripts should start with the word make
file.names <- file.names[-which(file.names %in% c("Main.R"))] # remove current file and Main (if different)
lapply(file.names, source) # source all functions needed

library(parallel)
load("guillemot.RData")
Time <- 16
states <- c("N","B1","LB","L_B","LB_","L_B_","S")
struc1 <- list("phi"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
               "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
               "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
               "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
               "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
               "p"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1,2:6,length(states))))

mv <- suff.stat(ch)

set.seed(846541147)
seeds <- sample(1:.Machine$integer.max,6)
cores <- 3
op <- mclapply(seeds,function(x){
  set.seed(x)
  starter <- logit(runif(6))
  opp <- optim(starter,
               ll.sl,
               phi.ind=1,
               delt.ind=2,
               kap.ind=3,
               rho.ind=4,
               gam.ind=5,
               p.ind=6,
               struc=struc1,
               mv=mv,
               control=list(maxit=10000,fnscale=-1),
               hessian=TRUE)
  return(c(opp,list("start"=starter)))
  },mc.cores=cores)

df1 <- data.frame("par"=rep(c("phi","delta","kappa","rho","gamma","p"),length(seeds)),
                  "mle"=unlist(lapply(1:length(seeds),function(x) logistic(op[[x]]$par))),
                  "value"=rep(sapply(1:length(seeds),function(x) op[[x]]$value),each=6),
                  "start"=unlist(lapply(1:length(seeds),function(x) logistic(op[[x]]$start))))
df1.sel <- df1[df1$value==min(df1$value),]

struc2 <- list("phi"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=as.list(1:length(states))),
               "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
               "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
               "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
               "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
               "p"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1,2:6,length(states))))

set.seed(846541147)
seeds2 <- sample(1:.Machine$integer.max,9)
op2 <- mclapply(seeds,function(x){
  set.seed(x)
  starter <- logit(runif(length(states)+5))
  opp <- optim(starter,
               ll.sl,
               phi.ind=1:length(states),
               delt.ind=length(states)+1,
               kap.ind=length(states)+2,
               rho.ind=length(states)+3,
               gam.ind=length(states)+4,
               p.ind=length(states)+5,
               struc=struc2,
               mv=mv,
               control=list(maxit=10000,fnscale=-1),
               hessian=TRUE)
  return(c(opp,list("start"=starter)))
},mc.cores=cores)


df2 <- data.frame("par"=rep(c(paste0("phi",states),"delta","kappa","rho","gamma","p"),length(seeds)),
                  "mle"=unlist(lapply(1:length(seeds),function(x) logistic(op2[[x]]$par))),
                  "value"=rep(sapply(1:length(seeds),function(x) op2[[x]]$value),each=length(states)+5),
                  "start"=unlist(lapply(1:length(seeds),function(x) logistic(op2[[x]]$start))))
df2.sel <- df2[df2$value==min(df2$value),]
df2.sel <- df2[df2$value==min(unique(df2$value)[-which.min(unique(df2$value))]),]

