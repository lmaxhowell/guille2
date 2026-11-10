file.names <- list.files(path = getwd(), pattern = "\\.R$") # list vector of file names
file.names <- file.names[-grep("^make",file.names)] # any non-function scripts should start with the word make
file.names <- file.names[-which(file.names %in% c("Main.R"))] # remove current file and Main (if different)
lapply(file.names, source) # source all functions needed

load("guillemot.RData")
Time <- 16
states <- c("N","B1","LB","L_B","LB_","L_B_","S")
struc <- list("phi"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "p"=list("age"=list(1:Time),"time"=list(1:6,7:Time),"state"=list(1,2:6,length(states))))
theta <- logit(c(0.6,0.5,0.4,0.2,0.7,0.8,0.4))
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

############################################################
# changing things to include skip as a state but fully unobservable
rm(list=ls())
file.names <- list.files(path = getwd(), pattern = "\\.R$") # list vector of file names
file.names <- file.names[-grep("^make",file.names)] # any non-function scripts should start with the word make
file.names <- file.names[-which(file.names %in% c("Main.R"))] # remove current file and Main (if different)
lapply(file.names, source) # source all functions needed

load("guillemot.RData")
Time <- 16
states <- c("N","B1","LB","L_B","LB_","L_B_","S")
struc <- list("phi"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "p"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1,2:6,length(states))))
theta <- logit(c(0.6,0.3,0.4,0.2,0.7,0.6))
ni.s <- list()
for(i in 1:length(ni2[1:15])){
  ni.s[[i]] <- rep("N",length(ni2[[i]]))
}

dat <- dat.sim.wrap(theta,1,2,3,4,5,6,struc,ni2[1:15],ni.s,789456132)

# tb <- table(unlist(dat))
# unname((tb["L_B"] + tb["L_B_"])/(tb["L_B"] + tb["L_B_"] + tb["LB"] + tb["LB_"])) # delta
# unname((tb["L_B"] + tb["LB"])/(tb["L_B"] + tb["L_B_"] + tb["LB"] + tb["LB_"])) # gamma
dg.rates(dat)

mv <- suff.stat(dat)
ll.sl(rep(0,6),1,2,3,4,5,6,struc,mv)
ll.sl(theta,1,2,3,4,5,6,struc,mv)
op <- optim(rep(0,length(theta)),
             ll.sl,
             phi.ind=1,
             delt.ind=2,
             kap.ind=3,
             rho.ind=4,
             gam.ind=5,
             p.ind=6,
             struc=struc,
             mv=mv,
             control=list(maxit=10000,fnscale=-1),
             hessian=TRUE)
op$convergence
logistic(op$par)

phi <- untrans(logistic(theta[1]),struc$phi$age,struc$phi$time,struc$phi$state)
delt <- untrans(logistic(theta[2]),struc$delt$age,struc$delt$time,struc$delt$state)
kap <- untrans(logistic(theta[3]),struc$kap$age,struc$kap$time,struc$kap$state)
rho <- untrans(logistic(theta[4]),struc$rho$age,struc$rho$time,struc$rho$state)
gam <- untrans(logistic(theta[5]),struc$gam$age,struc$gam$time,struc$gam$state)
p <- untrans(logistic(theta[6]),struc$p$age,struc$p$time,struc$p$state)
psi <- make.psi(delt,kap,rho,gam)

mv$m[1,1,1,1,3]
phi[1,1,1]*psi[1,1,1,1]*(1-p[2,2,1])*phi[2,2,1]*psi[1,1,2,2]*p[3,3,1]
Q.prob(1,1,1,1,3,phi,psi,p)

mv$m[3,4,9,9,11]
sum(sapply(1:length(states),function(x) {
  prob <- phi[9,9,3]*psi[3,x,9,9]*(1-p[10,10,x])*phi[10,10,x]*psi[x,4,10,10]*p[11,11,x]
  return(prob)}))
Q.prob(3,4,9,9,11,phi,psi,p)

mv$m[5,3,7,4,9]
sum(sapply(1:length(states),function(x) {
  prob <- phi[7,4,5]*psi[5,x,7,4]*(1-p[8,5,x])*phi[8,5,x]*psi[x,3,8,5]*p[9,6,x]
  return(prob)}))
Q.prob(5,3,7,4,9,phi,psi,p)

mv$m[1,1,11,1,14]
sum(sapply(1:length(states),function(x) {
  prob <- phi[11,1,1]*psi[1,x,11,1]*(1-p[12,2,x])*(sum(sapply(1:length(states),function(y) {
    prob <- phi[12,2,x]*psi[x,y,12,2]*(1-p[13,3,y])*phi[13,3,y]*psi[y,1,13,3]*p[14,4,1]
    # print(prob)
    return(prob)} )))
  # print(prob)
  return(prob)}))
Q.prob(1,1,11,1,14,phi,psi,p)

mv$v[1,14,2]
(1-phi[14,2,1])+phi[14,2,1]*psi[1,1,14,2]*(1-p[15,3,1])*((1-phi[15,3,1])+phi[15,3,1]*psi[1,1,15,3]*(1-p[16,4,1]))
Chi(1,14,2,phi,psi,p)

mv$v[2,13,3]
(1-phi[13,3,2]) + phi[13,3,2]*sum(sapply(2:7,function(x) 
  psi[2,x,13,3]*(1-p[14,4,x])*((1-phi[14,4,x])+phi[14,4,x]*sum(sapply(2:7,function(y) 
    psi[x,y,14,4]*(1-p[15,5,y])*((1-phi[15,5,y])+phi[15,5,y]*sum(sapply(2:7,function(z) 
      psi[y,z,15,5]*(1-p[16,6,z]) ))))))))
(1-phi[13,3,2]) + phi[13,3,2]*psi[2,7,13,3]*(1-phi[14,4,7])
Chi(2,13,3,phi,psi,p)

# all prob for individuals released at time 1 in state 1 at age 1
sum(sapply(1:length(states),function(x) sapply(14:Time,function(y) Q.prob(1,x,13,1,y,phi,psi,p)))) + Chi(1,13,1,phi,psi,p)
summingtime <- 0
for(x in 1:length(states)){
  for(y in 14:Time){
    summingtime <- summingtime + Q.prob(1,x,13,1,y,phi,psi,p)
    print(y)
  }
  print(x)
}
summingtime
Chi(1,13,1,phi,psi,p)
summingtime + Chi(1,13,1,phi,psi,p)


nrows <- 5:20
df <- data.frame("nrows"=nrows,"ll"=rep(0,length(nrows)),"time"=rep(0,length(nrows)))
for(i in 1:length(nrows)){
  mv <- suff.stat(dat[1:nrows[i],])
  time.start <- Sys.time()
  ll <- ll.sl(theta,1,2,3,4,5,6,struc,mv)
  time.end <- Sys.time()
  time.diff <- difftime(time.end,time.start,units="secs")
  df[i,] <- list(nrows[i],ll,time.diff) 
}

testvals <- 14:6
df <- data.frame("t"=testvals,"value"=rep(0,length(testvals)),"time"=rep(0,length(testvals)))
for(i in 1:length(testvals)){
  time.start <- Sys.time()
  val <- Chi(2,testvals[i],3,phi,psi,p)
  time.end <- Sys.time()
  time.diff <- difftime(time.end,time.start,units="secs")
  df[i,] <- list(testvals[i],val,time.diff)
  print(i)
}
df <- cbind(df,"logval"=log(df$value),"logdiff"=log(df$value)-min(log(df$value)))
ggplot(df,aes(time,value)) + geom_point()
ggplot(df,aes(t,time,col=value)) + geom_point()
ggplot(df,aes(t,value,col=time)) + geom_point()

testvals <- 14:6
df2 <- data.frame("t"=rep(testvals,each=length(states)),
                  "r"=rep(states,length(testvals)),
                  "value"=rep(0,length(testvals)*length(states)),
                  "time"=rep(0,length(testvals)*length(states)))
cr <- 1
for(i in 1:length(testvals)){
  for(r in 1:length(states)){
    time.start <- Sys.time()
    val <- Chi(r,testvals[i],3,phi,psi,p)
    time.end <- Sys.time()
    time.diff <- difftime(time.end,time.start,units="secs")
    df2[cr,] <- list(testvals[i],states[r],val,time.diff)
    cr <- cr + 1
    cat("\r",c(i,r))
  }
}
ggplot(df2,aes(t,time,col=value)) + geom_point() + facet_wrap(~r)
ggplot(df2,aes(r,time)) + geom_boxplot()

pll <- 15
ax <- seq(0.01,0.99,length.out=pll)
ll <- array(0,dim=rep(pll,6))
# ll <- array(0,dim=rep(pll,4))
timer(
for(a in 1:pll){
  for(b in 1:pll){
    for(cc in 1:pll){
      for(d in 1:pll){
        for(e in 1:pll){
          for(f in 1:pll){
            ll[a,b,cc,d,e,f] <- ll.sl(logit(c(ax[a],ax[b],ax[cc],ax[d],ax[e],ax[f])),1,2,3,4,5,6,struc,mv)
            # ll[a,b,e,f] <- ll.sl(logit(c(ax[a],0.3,ax[cc],ax[d],0.7,ax[f])),1,2,3,4,5,6,struc,mv)
            cat("\r",c(a,b,cc,d,e,f))
            # cat("\r",c(a,cc,d,f))
          }
        }
      }
    }
  }
}) #Time difference of 2.347487 days

library(reshape2)
df3 <- melt(ll,varnames = c("phi","delta","kappa","rho","gamma","p"))
for(i in 1:pll){
  df3[df3==i] <- ax[i]
  cat("\r",i)
}
df4 <- df3
for(i in 1:6){
  df4[,i] <- factor(df4[,i],levels=ax)
}
ggplot(df3,aes(phi,delta,z=value,col=gamma,linetype=p)) + geom_contour() + facet_grid(vars(kappa),vars(rho))




