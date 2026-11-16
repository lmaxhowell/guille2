file.names <- list.files(path = getwd(), pattern = "\\.R$") # list vector of file names
file.names <- file.names[-grep("^make",file.names)] # any non-function scripts should start with the word make
file.names <- file.names[-which(file.names %in% c("Main.R"))] # remove current file and Main (if different)
lapply(file.names, source) # source all functions needed


# constant model
Time <- 16
states <- c("N","B1","LB","L_B","LB_","L_B_","S")

# true values going phi,delta,kappa,rho,gamma,p
theta <- logit(c(0.45,0.5,0.4,0.2,0.75,0.6))
par.len <- length(theta)

# the number of individuals released at each time point. it can very year to year
released <- rep(10,Time-1) # must be a vector of length Time-1
# released <- c(10,10,10,10,10,10,10,10,10,10,10,10,10,10,10) # would also work

# initial age for each individual
# doing all juveniles here but doesn't have to be
ni <- lapply(1:(Time-1),function(x) rep(1,released[x])) 

# initial state for each individual
ni.s <- lapply(1:(Time-1),function(x) rep("N",released[x])) 

# what is the structure of the parameters. this corresponds to constant model
# can see that p has to vary due to constraints but is constant for our definitions
# theoreically dont need p$state here as this is constrained within the likelihood functions but
# is here for constistancy
struc <- list("phi"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "p"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1,2:6,length(states))))

seed <- sample(1:.Machine$integer.max,1) # I use this to avoid having to come up with seeds myself
# and in the case of setting multiple seeds for a simulation I can set one seed and then generate seeds using that
# but equally
# seed <- 1216383180 # would work

# generate capture histories
dat <- dat.sim.wrap(theta,1,2,3,4,5,6,struc,ni,ni.s,seed)

# turn that into sufficient statistics
mv <- suff.stat(dat)

# optimise for mle
op <- optim(rep(0,par.len),
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
