file.names <- list.files(path = getwd(), pattern = "\\.R$") # list vector of file names
file.names <- file.names[-grep("^make",file.names)] # any non-function scripts should start with the word make
file.names <- file.names[-which(file.names %in% c("Main.R"))] # remove current file and Main (if different)
lapply(file.names, source) # source all functions needed

load("guillemot.RData")
Time <- 16
states <- c("N","B1","LB","L_B","LB_","L_B_")
struc <- list("phi"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "p"=list("age"=list(1:Time),"time"=list(1:6,7:Time),"state"=list(1,2:length(states))))
theta <- logit(c(0.6,0.5,0.4,0.2,0.7,0.4,0.9))
ni.s <- list()
for(i in 1:length(ni2[1:15])){
  ni.s[[i]] <- rep("N",length(ni2[[i]]))
}
dat <- dat.sim.wrap(theta,1,2,3,4,5,6:7,struc,ni2[1:15],ni.s,789456132)

mv <- suff.stat(dat)
ll.sl(rep(0,7),1,2,3,4,5,6:7,struc,mv)
ll.sl(theta,1,2,3,4,5,6:7,struc,mv)
op <- optim(rep(0,length(theta)),
            ll.sl,
            phi.ind=1,
            delt.ind=2,
            kap.ind=3,
            rho.ind=4,
            gam.ind=5,
            p.ind=6:7,
            struc=struc,
            mv=mv,
            control=list(maxit=10000,fnscale=-1),
            hessian=TRUE)
op$convergence
logistic(op$par)
solve(-op$hessian)
logistic(ci(op$par,diag(solve(-op$hessian))))

phi <- untrans(c(0.3,0.6,0.9),list(1,2,3:Time),list(1:(Time-1)),list(1:length(states)))
del <- untrans(c(0.4),list(1:Time),list(1:(Time-1)),list(1:length(states)))
kap <- untrans(c(0.2),list(1:Time),list(1:(Time-1)),list(1:length(states)))
rho <- untrans(c(0.1,0.4,0.7),list(1:3,4,5:Time),list(1:(Time-1)),list(1:length(states)))
gam <- untrans(c(0.75),list(1:Time),list(1:(Time-1)),list(1:length(states)))
psi <- make.psi(del,kap,rho,gam)
p  <- untrans(c(0.8,0),list(1:Time),list(1:(Time)),list(1,2:length(states)))

Q.prob(1,1,4,1,phi,psi,p)
phi[1,1,1]*psi[1,1,1,1]*(1-p[2,2,1])*phi[2,2,1]*psi[1,1,2,2]*(1-p[3,3,1])*phi[3,3,1]*psi[1,1,3,3]*p[4,4,1]

Y.fun(1,1,4,phi,psi,p)
phi[1,1,1]*psi[1,1,1,1]*(1-p[2,2,1])*phi[2,2,1]*psi[1,1,2,2]*(1-p[3,3,1])

ch[14,]
mv2 <- suff.stat(ch[14,])
struc2 <- list("phi"=list("age"=list(1,2,3:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "rho"=list("age"=list(1:3,4,5:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "p"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1,2:6,length(states))))
theta2 <- logit(c(0.3,0.6,0.9,0.4,0.2,0.1,0.4,0.7,0.75,0.8))
phi <- untrans(logistic(theta2[1:3]),struc2$phi$age,struc2$phi$time,struc2$phi$state)
delt <- untrans(logistic(theta2[4]),struc2$delt$age,struc2$delt$time,struc2$delt$state)
kap <- untrans(logistic(theta2[5]),struc2$kap$age,struc2$kap$time,struc2$kap$state)
rho <- untrans(logistic(theta2[6:8]),struc2$rho$age,struc2$rho$time,struc2$rho$state)
gam <- untrans(logistic(theta2[9]),struc2$gam$age,struc2$gam$time,struc2$gam$state)
psi <- make.psi(delt,kap,rho,gam)
p <- untrans(c(logistic(theta2[10]),1),struc2$p$age,struc2$p$time,struc2$p$state)

ll <- ll.sl(theta2,1:3,4,5,6:8,9,10,struc2,mv2) # changing the function here to return probs matrices as well

ll$probs.v[3,10,10]
mv2$v[3,10,10] # r,t,a
log(phi[10,10,3]*psi[3,7,10,10]*(1-phi[11,11,7]) + (1-phi[10,10,3]))

ll$probs.mN[1,1,1,4] # r,t,a,t'
phi[1,1,1]*psi[1,1,1,1]*(1-p[2,2,1])*phi[2,2,1]*psi[1,1,2,2]*(1-p[3,3,1])*phi[3,3,1]*psi[1,1,3,3]*p[4,4,1]
log(phi[1,1,1]*psi[1,1,1,1]*(1-p[2,2,1])*phi[2,2,1]*psi[1,1,2,2]*(1-p[3,3,1])*phi[3,3,1]*psi[1,1,3,3]*p[4,4,1])
log(Q.prob(1,1,4,1,phi,psi,p))


ll$probs.mN[1,4,4,5]
mv2$mN[1,4,4,5] # r,t,a,t'
phi[5-1,4+5-4-1,1]*psi[1,1,5-1,4+5-4-1]*p[5,4+5-4,1]
log(phi[5-1,4+5-4-1,1]*psi[1,1,5-1,4+5-4-1]*p[5,4+5-4,1])
log(Q.prob(4,4,5,1,phi,psi,p))

ll$probs.mN[2,6,6,7]
mv2$mN[2,6,6,7] # r,t,a,t'
phi[7-1,6+7-6-1,2]*psi[1,2,7-1,6+7-6-1]*p[7,6+7-6,2]
log(phi[7-1,6+7-6-1,2]*psi[1,2,7-1,6+7-6-1]*p[7,6+7-6,2])
log(Q.prob(6,6,7,2,phi,psi,p))

ll$probs.mB[2,5,7,7]  # r,s,t,a
phi[7,7,2]*psi[2,5,7,7]
log(phi[7,7,2]*psi[2,5,7,7])




# try simulation with survival==1 and survival==0
# big population constant model to see if it can do that
# what changes if we model the skip state differently (unobservable state)


struc3 <- list("phi"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "p"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1,2:6,length(states))))
theta3 <- logit(c(0.6,0.5,0.4,0.2,0.7,0.6))
dat3 <- dat.sim.wrap(theta3,1,2,3,4,5,6,struc3,ni2[1:15],ni.s,789456132)

mv3 <- suff.stat(dat3)
ll.sl(rep(0,6),1,2,3,4,5,6,struc3,mv3)
ll.sl(theta3,1,2,3,4,5,6,struc3,mv3)
op3 <- optim(rep(0,length(theta3)),
             ll.sl,
             phi.ind=1,
             delt.ind=2,
             kap.ind=3,
             rho.ind=4,
             gam.ind=5,
             p.ind=6,
             struc=struc3,
             mv=mv3,
             control=list(maxit=10000,fnscale=-1),
             hessian=TRUE)
op3$convergence
logistic(op3$par)

df.t <- data.frame("t"=rep(NA,length(2:15)),"chi"=rep(NA,length(2:15)),"time"=rep(NA,length(2:15)))
for(t in 15:2){
  i <- 15-t+2-1
  t.start <- Sys.time()
  val <- Chi(2,t,2,phi,psi,p)
  t.end <- Sys.time()
  t.diff <- difftime(t.end,t.start,units="secs")
  df.t[i,] <- list(15-t+2,val,t.diff)
  print(t)
  if(t.diff>90){
    for(q in t:2){
      df.t[15-q+2,] <- list(15-q+2,val,t.diff)
    }
    break()
  }
}
ggplot(df.t,aes(t,time)) + geom_line()
ggplot(df.t,aes(t,chi)) + geom_line()

df.t2 <- data.frame("r"=rep(states[2:6],each=length(2:15)),"t"=rep(2:15,length(2:6)),"chi"=rep(NA,length(2:6)*length(2:15)),"time"=rep(NA,length(2:6)*length(2:15)))
for(r in 2:length(states)){
  i <- seq(1,84,14)[r-1]
  for(t in 15:2){
    t.start <- Sys.time()
    val <- Chi(r,t,2,phi,psi,p)
    t.end <- Sys.time()
    t.diff <- difftime(t.end,t.start,units="secs")
    df.t2[i+15-t,3:4] <- list(val,t.diff)
    print(c(r,t,i))
    if(t.diff>90){
      for(q in t:2){
        df.t2[i+15-q,3:4] <- list(val,t.diff)
      }
      break()
    }
  }
}
ggplot(df.t2,aes(t,time,col=r)) + geom_line()
ggplot(df.t2,aes(t,chi,col=r)) + geom_line()
 