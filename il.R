make.psi <- function(delta,kap,rho,gam){
  states <- c("N","B1","LB","L_B","LB_","L_B_","S")
  # all five input arrays should be indexed by [time,age,state]
  Time <- dim(delta)[1]
  Ages <- dim(delta)[2]
  # print(sapply(list(delta,kap,rho,gam,epsilon),ncol))
  if(nrow(kap)!=Time | nrow(rho)!=Time | nrow(gam)!=Time){
    stop("delta,kappa,rho and gamma must all have the same number of rows, corresponding to time")
  }
  psi <- array(0,dim=c(length(states),length(states),Time,Ages))
  for(t in 1:Time){
    for(a in 1:Ages){
      # psi[1,,t,a] <- c((1-epsilon[t,a,1])*(1-rho[t,a,1]),epsilon[t,a,2],(1-epsilon[t,a,3])*rho[t,a,3],0,0,0,0,0)
      # psi[2,,t,a] <- c((1-epsilon[t,a,1])*(1-rho[t,a,1]),epsilon[t,a,2],(1-epsilon[t,a,3])*rho[t,a,3],0,0,0,0,0)
      psi[1,,t,a] <- c((1-rho[t,a,1]),rho[t,a,1],0,0,0,0,0)
      row.names(psi) <- states
      colnames(psi) <- states
      # this is if these rows depend on the ROWS theyre in
      for(i in 2:6){
        psi[i,,t,a] <- c(0,0,
                           (1-kap[t,a,i])*(1-delta[t,a,i])*gam[t,a,i],
                           (1-kap[t,a,i])*delta[t,a,i]*gam[t,a,i],
                           (1-kap[t,a,i])*(1-delta[t,a,i])*(1-gam[t,a,i]),
                           (1-kap[t,a,i])*delta[t,a,i]*(1-gam[t,a,i]),
                           kap[t,a,i])
      }
      psi[7,,t,a] <- c(0,0,
                       (1-delta[t,a,7])*gam[t,a,7],
                       delta[t,a,7]*gam[t,a,7],
                       (1-delta[t,a,7])*(1-gam[t,a,7]),
                       delta[t,a,7]*(1-gam[t,a,7]),
                       0)
      # # this is if these rows depend on the columns theyre in
      # for(i in 4:8){
      #   psi[i-1,,t,a] <- c(0,0,0,
      #                      (1-kap[t,a,4])*(1-delta[t,a,4])*gam[t,a,4],
      #                      (1-kap[t,a,5])*delta[t,a,5]*gam[t,a,5],
      #                      (1-kap[t,a,6])*(1-delta[t,a,6])*(1-gam[t,a,6]),
      #                      (1-kap[t,a,7])*delta[t,a,7]*(1-gam[t,a,7]),
      #                      kap[t,a,8])
      # }
      # psi[8,,t,a] <- c(0,0,0,
      #                  (1-delta[t,a,4])*gam[t,a,4],
      #                  delta[t,a,5]*gam[t,a,5],
      #                  (1-delta[t,a,6])*(1-gam[t,a,6]),
      #                  delta[t,a,7]*(1-gam[t,a,7]),
      #                  0)
    }
  }
  return(psi)
}

find.transitions <- function(ch){ # need to add in age as component
  age <- as.numeric(ch[length(ch)])
  ch <- ch[1:(length(ch)-1)]
  Time <- length(ch)
  where <- which(ch!="0")
  df <- as.data.frame(array(NA,dim=c(max(length(where)-1,1),5))) # want it to be a character df but full of nothing
  # df <- as.data.frame(array(NA,dim=c(length(where),5))) # want it to be a character df but full of nothing
  colnames(df) <- c("r","s","t_r","t_s","age") # state r, state s, time at state r, time at states
  if(nrow(df)==1){ # if only one observation
    if(length(where)==1){ # if only one state observed
      df[1,] <- list(ch[ch!="0"],"0",where,where+1,age)
      return(df)
    }else if(length(where)==2 && where[2]==Time){ # the case where there is a state at the final time
      df[1,] <- list(ch[where[1]],ch[where[2]],where[1],where[2],age)
      return(df)
    }else if(length(where)==2 && where[2]!=Time){
      df[1,] <- list(ch[where[1]],ch[where[2]],where[1],where[2],age)
      lw <- where[2]
      df <- rbind(df,list(ch[lw],"0",lw,lw+1,age+length(where)-1))
      return(df)
    }
  }else{ # more than one observation
    for(t in 1:(length(where)-1)){ # go through all the states
      df[t,] <- list(ch[where[t]],ch[where[t+1]],where[t],where[t+1],age+t-1)
      # the age+t-1 ONLY WORKS if the capture history has NO gaps once the bird is seen
    } # end for t in 1:length(where)-1
    lw <- where[length(where)] # the last element of where
    if(lw!=Time){ # if the last place that a state is observed ISNT the last time
      df <- rbind(df,
                  list(ch[lw],"0",lw,lw+1,age+length(where)-1))
    }
    return(df)
  } # end else more than one observation
}

Pr_rs <- function(r,s,t,a,phi,psi){
  # print(c("r",r))
  # print(c("s",s))
  # print(c("t",t))
  # print(c("a",a))
  prob <- phi[t,a,r]*psi[r,s,t,a]
  return(prob)
}

Chi <- function(r,t,a,phi,psi,p){
  # print(c(r,t,a))
  Time <- dim(phi)[1]+1
  if(t==Time){
    return(1)
  }else{
    if(r==1){
      prob <- 1-phi[t,a,r] + phi[t,a,r]*psi[r,1,t,a]*(1-p[t+1,a+1,1])*Chi(1,t+1,a+1,phi,psi,p)
    }else{
      
    }
    
    if(length(prob)==0){
      print(parent.frame(2)$ch)
    }
    return(prob)
  }
}

Y.fun <- function(t,a,tprime,phi,psi,p){
  prob <- prod(sapply(t:(tprime-2),function(x) phi[x,a+x-t,1]*psi[1,1,x,a+x-t]*(1-p[x+1,a+x-t+1,1]) ))
  return(prob)
}

Q.prob <- function(t,a,tprime,r,phi,psi,p){
  if((tprime-t)==1){
    prob <- phi[tprime-1,a+tprime-t-1,1]*psi[1,r,tprime-1,a+tprime-t-1]*p[tprime,a+tprime-t,r]
  }else if((tprime-t)>1){
    prob <- Y.fun(t,a,tprime,phi,psi,p)*phi[tprime-1,a+tprime-t-1,1]*psi[1,r,tprime-1,a+tprime-t-1]*p[tprime,a+tprime-t,r]
  }
  return(prob)
}

Pr_r0 <- function(r,t,a,phi,psi){
  skip <- which(row.names(psi)=="S")
  if(r %in% 3:7){
    if(t<(Time-1)){
      prob <- phi[t,a,r]*psi[r,skip,t,a]*(1-phi[t+1,a+1,skip]) + (1-phi[t,a,r])
    }else if(t==(Time-1)){
      # if we are in a breeding state at Time-1 and then a zero at Time
      # then they dont have to have died they could still be in the skipping state
      # with no information about them being alive
      prob <- phi[t,a,r]*psi[r,skip,t,a] + (1-phi[t,a,r])
    }
    # prob <- phi[t,a,r]*psi[r,skip,t,a]*(1-phi[t+1,a+1,skip]) + (1-phi[t,a,r])
  # }else if(r==1){
  }else if(r %in% 1:2){
    prob <- Chi(r,t,a,phi,psi)
  }else{
    # stop(print("I shouldnt be here",r,t,a))
  }
  if(!exists("prob")){
    print(c(r,t,a,as.vector(parent.frame()$ch),as.vector(parent.frame()$transitions)))
  }
  return(prob)
}

il <- function(ch,phi,psi){ # il is "individual likelihood"
  # phi is indexed by [time,age,state]
  # psi is indexed by [state_r,state_s,time,age]
  Indicator <- function(prob){
    # a function that checks if the probability pu into it will
    # "cause problems" in the likelihood - i.e. if we log
    # this probability will it be infinite and cause the whole ll
    # to be infinite as a result
    # or has a small amount of numerical error crept in
    # and made a probability negative?
    # in which case just set it to zero
    # will add a warning in in this case
    # this function just means we can call log(prob)
    # without worrying about the errors this can cause
    if(length(prob)==0){
      print(c(parent.frame()$i,parent.frame()$t,parent.frame()$Time))
      print(c("r",parent.frame()$current_state_index))
      print(c("s",parent.frame()$next_state_index))
      print(c("t",parent.frame()$current_time))
      print(c("a",parent.frame()$current_age))
      print(parent.frame()$ch)
      View(parent.frame()$transitions)
    }
    if(!exists("prob") | is.na(prob)){
      print(parent.frame()$ch)
      print(parent.frame()$transitions)
    }
    if(prob<0){
      return(0)
    }else{
      return(ifelse(is.finite(log(prob)),log(prob),0)) 
    }
  }
  library(plyr, include.only = c("count"))
  # just need the count function from this package
  # to figure out how many unique chs there are
  uch <- count(ch) # unique number of capture histories
  # what are the states
  states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
  # number of columns of the chs minus is Time to account for the age column
  Time <- ncol(ch)-1
  # how many unique capture histories are there
  Individuals <- nrow(uch)
  ll_i <- rep(0,Individuals) # want to create an individual likelihood for each capture history
  # L1 <- rep(0,Individuals)
  # L2 <- rep(0,Individuals)
  for(i in 1:Individuals){
    ch <- uch[i,1:(Time+1)] # the current capture history we are looking at
    # can rename ch to be "the current ch" because we no longer need ch, just uch
    transitions <- find.transitions(ch)
    # print(transitions)
    # need to deal with all the capture histories that
    # only have one state in them first
    if(nrow(transitions)==1){
      current_state <- transitions[1,1]
      current_time <- transitions[1,3]
      if(transitions[1,4]==(Time+1)){ # if this occurs at the final time...do nothing
        next
      }else if(transitions[1,2]!="0"){
        current_state_index <- which(states==current_state) # what number is the current state
        current_age <- transitions[1,5]
        next_state_index <- which(states==transitions[1,2]) # what number is the next state
        # L1[i] <- L1[i] + Indicator(Pr_rs(current_state_index,next_state_index,current_time,current_age,phi,psi))
        ll_i[i] <- ll_i[i] + Indicator(Pr_rs(current_state_index,next_state_index,current_time,current_age,phi,psi))
      }else{ # else we are in a Pr_r0 situation
        current_state_index <- which(states==current_state) # what number is the current state
        current_age <- transitions[1,5]
        # print(c(current_state_index,current_time,current_age,dim(phi),dim(psi)))
        # print(Indicator(Pr_r02(current_state_index,current_time,current_age,phi,psi)))
        # L2[i] <- L2[i] + Indicator(Pr_r0(current_state_index,current_time,current_age,phi,psi))
        ll_i[i] <- ll_i[i] + Indicator(Pr_r0(current_state_index,current_time,current_age,phi,psi))
      }
    }else{ # therefore ch has multiple observations
      for(t in 1:nrow(transitions)){
        if(transitions[t,4]==(Time+1)){ # in the situation where the final time point is an observed state and not a zero
          next # then we just skip it
        }else if((transitions[t,4]-transitions[t,3])==1){ # aka are the states next to each other and not a state in the final time
          current_state_index <- which(states==transitions[t,1]) # what number is the current state
          current_time <- transitions[t,3]
          current_age <- transitions[t,5]
          # need to deal with if the transition is the final one, a state to zero
          if(transitions[t,2]=="0"){
            # L2[i] <- L2[i] + Indicator(Pr_r0(current_state_index,current_time,current_age,phi,psi))
            # print(transitions[nrow(transitions),])
            ll_i[i] <- ll_i[i] + Indicator(Pr_r0(current_state_index,current_time,current_age,phi,psi))
          }else{
            next_state_index <- which(states==transitions[t,2]) # what number is the next state
            # print(c(Pr_rs2(current_state_index,next_state_index,current_time,current_age,phi,psi),current_state_index,next_state_index,current_time,current_age))
            # L1[i] <- L1[i] + Indicator(Pr_rs(current_state_index,next_state_index,current_time,current_age,phi,psi))
            ll_i[i] <- ll_i[i] + Indicator(Pr_rs(current_state_index,next_state_index,current_time,current_age,phi,psi))
          } # end else
        } # all the states should now be next to each other so commented out the next bit
        # }else if((transitions[t,4]-transitions[t,3])==2 & transitions[t,1]!="N"){ # we have a skipped a time and it was previously breeding
        #   current_state_index <- which(states==transitions[t,1]) # what number is the current state
        #   next_state_index <- which(states==transitions[t,2]) # what number is the next state
        #   current_time <- transitions[t,3]
        #   ll_i[i] <- ll_i[i] + Indicator(Pr_r0s(current_state_index,next_state_index,current_time,phi,psi))
        # }else if(((transitions[t,4]-transitions[t,3])==2 & transitions[t,1]=="N") | 
        #          ((transitions[t,4]-transitions[t,3])>2)){
        #   # we have a skipped a time and it was NOT previously breeding OR
        #   # we have skipped multiple times (and so theoretically should not be breeding) OR
        #   current_state_index <- which(states==transitions[t,1]) # what number is the current state
        #   current_time <- transitions[t,3]
        #   ll_i[i] <- ll_i[i] + Indicator(Pr_r0(current_state_index,current_time,phi,psi))
        # }
        # print(transitions[t,])
        # print(c(t,ll_i[i]))
      } # end for t in transitions
    } # end else ch has multiple observations
    # print(i)
  } # end i
  # LL1 <- sum((uch$freq)*(L1))
  # LL2 <- sum((uch$freq)*(L2))
  ll <- sum((uch$freq)*(ll_i)) # multiply the log-likelihood for each unique capture history
  # by how many times that capture history occurs and sum it up
  # to make the multinomial loglikelihood
  return(ll)
  # return(c(ll,LL1,LL2,LL1+LL2))
}

ll.il <- function(theta,phi.ind,delt.ind,kap.ind,rho.ind,gam.ind,p.ind,struc,ch){
  states <- c("N","B1","LB","L_B","LB_","L_B_","S")
  
  phi <- untrans(logistic(theta[phi.ind]),struc$phi$age,struc$phi$time,struc$phi$state)
  delt <- untrans(logistic(theta[delt.ind]),struc$delt$age,struc$delt$time,struc$delt$state)
  kap <- untrans(logistic(theta[kap.ind]),struc$kap$age,struc$kap$time,struc$kap$state)
  rho <- untrans(logistic(theta[rho.ind]),struc$rho$age,struc$rho$time,struc$rho$state)
  gam <- untrans(logistic(theta[gam.ind]),struc$gam$age,struc$gam$time,struc$gam$state)
  p <- untrans(logistic(theta[p.ind]),struc$p$age,struc$p$time,struc$p$state)
  
  psi <- make.psi(delt,kap,rho,gam,eps,p)
  # print(psi)
  
  ll <- il(ch,phi,psi)
  
  return(ll)
}

source("suff.stat.R")

sl <- function(phi,psi,p,mv){ # sl is sufficient statistic likelihood
  Indicator <- function(prob){
    # a function that checks if the probability pu into it will
    # "cause problems" in the likelihood - i.e. if we log
    # this probability will it be infinite and cause the whole ll
    # to be infinite as a result
    # or has a small amount of numerical error crept in
    # and made a probability negative?
    # in which case just set it to zero
    # will add a warning in in this case
    # this function just means we can call log(prob)
    # without worrying about the errors this can cause
    if(length(prob)==0){
      print(c(parent.frame()$i,parent.frame()$t,parent.frame()$Time))
      print(c("r",parent.frame()$current_state_index))
      print(c("s",parent.frame()$next_state_index))
      print(c("t",parent.frame()$current_time))
      print(c("a",parent.frame()$current_age))
      print(parent.frame()$ch)
      View(parent.frame()$transitions)
    }
    if(!exists("prob") | is.na(prob)){
      print(parent.frame()$ch)
      print(parent.frame()$transitions)
    }
    if(prob<0){
      return(0)
    }else{
      return(ifelse(is.finite(log(prob)),log(prob),0)) 
    }
  }
  states <- c("N","B1","LB","L_B","LB_","L_B_","S")
  mB <- mv$mB
  mN <- mv$mN
  v <- mv$v
  Time <- dim(mv$mB)[3]
  Age <- Time
  
  if(is.null(unlist(mv$violations))==FALSE){
    warning("some transitions not allowed")
  }
  
  probs.mB <- array(0,dim=dim(mB)) # r,s,t,a
  probs.mN <- array(0,dim=dim(mN)) # r,t,a,t'
  probs.v <- array(0,dim=dim(v)) # r,t,a
  
  skip <- which(row.names(psi)=="S")
  for(t in 1:Time){
    for(a in 1:Age){
      for(r in 1:length(states)){
        for(s in 1:length(states)){
          if(mB[r,s,t,a]>0){
            probs.mB[r,s,t,a] <- log(phi[t,a,r]*psi[r,s,t,a]*p[t+1,a+1,s])
            if(is.infinite((probs.mB[r,s,t,a]))){
              print(c(r,s,t,a))
            }
            # probs.mB[r,s,t,a] <- log(phi[t,a,r]*psi[r,s,t,a])
          }
        }
        if(r==1){
          if(v[r,t,a]>0){
            # print(c(r,t,a))
            probs.v[r,t,a] <- log(Chi(r,t,a,phi,psi,p))
          }
        }else if(r %in% 2:6){
          if(t<(Time-1)){
            if(v[r,t,a]>0){
              # print(c(r,t,a))
              # print(c(r,t,a))
              # print(c(phi[t,a,r],psi[r,skip,t,a],(1-phi[t+1,a+1,skip]),(1-phi[t,a,r])))
              probs.v[r,t,a] <- log(phi[t,a,r]*psi[r,skip,t,a]*(1-phi[t+1,a+1,skip]) + (1-phi[t,a,r]))
              # print(probs.v[r,t,a])
            }
          }else if(t==(Time-1)){
            # if we are in a breeding state at Time-1 and then a zero at Time
            # then they dont have to have died they could still be in the skipping state
            # with no information about them being alive
            if(v[r,t,a]>0){
              # print(c(r,t,a))
              # probs.v[r,t,a] <- log(phi[t,a,r]*psi[r,skip,t,a]*p[t+1,a+1,skip] + (1-phi[t,a,r]))
              probs.v[r,t,a] <- log(phi[t,a,r]*psi[r,skip,t,a] + (1-phi[t,a,r]))
            }
          } # end elif t==Time-1
        } # end if r in 2:6
      } # end r in States
      for(tprime in t:Time){
        for(r in 1:2){
          if(mN[r,t,a,tprime]>0){
            # print(c(r,t,a,tprime))
            probs.mN[r,t,a,tprime] <- log(Q.prob(t,a,tprime,r,phi,psi,p))
          }
        } # end r  in 1:2
      } # end tprime in t:Time
    } # end a in Age
  } # end t in Time
  # print(c(dim(probs.mB),0,dim(mB)))
  # print(c(dim(probs.mN),0,dim(mN)))
  # print(c(dim(probs.v),0,dim(v)))
  ll <- probs.mB%*%mB + probs.mN%*%mN + probs.v%*%v
  return(as.numeric(ll))
  # return(list("probs.mB"=probs.mB,"probs.mN"=probs.mN,"probs.v"=probs.v,"ll"=as.numeric(ll)))
}

ll.sl <- function(theta,phi.ind,delt.ind,kap.ind,rho.ind,gam.ind,p.ind,struc,mv){
  states <- c("N","B1","LB","L_B","LB_","L_B_","S")
  
  phi <- untrans(logistic(theta[phi.ind]),struc$phi$age,struc$phi$time,struc$phi$state)
  delt <- untrans(logistic(theta[delt.ind]),struc$delt$age,struc$delt$time,struc$delt$state)
  kap <- untrans(logistic(theta[kap.ind]),struc$kap$age,struc$kap$time,struc$kap$state)
  rho <- untrans(logistic(theta[rho.ind]),struc$rho$age,struc$rho$time,struc$rho$state)
  gam <- untrans(logistic(theta[gam.ind]),struc$gam$age,struc$gam$time,struc$gam$state)
  # if(length(struc$p$age)==1 && length(struc$p$time)==1){ # both are length one
  #   onesandzeros <- c(1,0)
  # }else if(length(struc$p$age)==1 || length(struc$p$time)==1){ # one is length one (can use or rather than xor as it is below the and statement)
  #   rl <- max(length(struc$p$age),length(struc$p$time))
  #   onesandzeros <- c(rep(1,rl),rep(0,rl))
  # }else{ # else they are both longer than length one
  #   onesandzeros <- c(rep(1,length(struc$p$age)+length(struc$p$time)),rep(0,length(struc$p$age)+length(struc$p$time)))
  # }
  
  onesandzeros <- c(rep(1,length(p.ind)),rep(0,length(p.ind)))
  p <- untrans(c(logistic(theta[p.ind]),onesandzeros),struc$p$age,struc$p$time,list(1,2:6,length(states)))
  # print(p[,1,])
  # print(p[,2,])
  # print(p[,3,])
  print(phi[1,,])
  print(phi[2,,])
  # print(p.ind)
  # print(theta[p.ind])
  # if(all(unlist(Map(`%in%`,struc$p$state,list(1,2:length(states)))))){
  #   stop("p must be 1 for non-breeding states")
  # }
  
  psi <- make.psi(delt,kap,rho,gam)
  # print(psi[,,1,1])
  # print(p[,,1])
  # print(phi[,,1])
  
  # ll <- sl(phi,psi,p,mv)$ll
  ll <- sl(phi,psi,p,mv)
  
  return(ll)
}

phi.trans <- function(alpha,beta,beta.struc){
  # inputs should be vectors
  alpha <- c(0,alpha) # to constrain alpha_0=0
  Ages <- length(alpha)
  Time_1 <- length(unlist(beta.struc)) # cause this is Time-1
  phi <- array(0,dim=c(Time_1,Time_1+1,length(states)))
  for(a in 1:(Time_1+1)){
    for(t in 1:Time_1){
      # print(c(t,a,min(a,Ages)))
      for(i in 1:length(beta.struc)){
        if(t %in% beta.struc[[i]]){
          break
        }
      }
      phi[t,a,] <- alpha[min(a,Ages)]+beta[i]
    }
  }
  return(phi)
}

ll.sl2 <- function(theta,ageclasses,timeclasses,beta.struc,mv){
  states <- c("N","B1","LB","L_B","LB_","L_B_","S")
  
  phi.trans <- function(alpha,beta,beta.struc){
    # inputs should be vectors
    alpha <- c(0,alpha) # to constrain alpha_0=0
    Ages <- length(alpha)
    Time_1 <- length(unlist(beta.struc)) # cause this is Time-1
    phi <- array(0,dim=c(Time_1,Time_1+1,length(states)))
    for(a in 1:(Time_1+1)){
      for(t in 1:Time_1){
        # print(c(t,a,min(a,Ages)))
        for(i in 1:length(beta.struc)){
          if(t %in% beta.struc[[i]]){
            break
          }
        }
        phi[t,a,] <- alpha[min(a,Ages)]+beta[i]
      }
    }
    return(phi)
  }
  
  rho.trans <- function(rho.input){
    rho <- array(0,dim=c(Time-1,Time,length(states)))
    # rho[,1,] <- 0 # sets as default anyway # age 0
    rho[,2:3,] <- rho.input[1] # age one and two
    rho[,4,] <- rho.input[2] # age 3
    rho[,5:Time,] <- rho.input[3] # age 4+
    return(rho)
  }
  
  alpha <- theta[1:(ageclasses-1)]
  beta <- theta[ageclasses:(timeclasses+ageclasses-1)]
  rho.input <- theta[(timeclasses+ageclasses):(timeclasses+ageclasses+2)]
  phi <- logistic(phi.trans(alpha,beta,beta.struc))
  
  rho <- logistic(rho.trans(rho.input))
  
  Time <- dim(mv$mB)[3]
  
  struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
  
  theta <- theta[-c(1:(timeclasses+ageclasses+2))]
  
  delt <- untrans(logistic(theta[1]),struc.const$age,struc.const$time,struc.const$state)
  kap <- untrans(logistic(theta[2]),struc.const$age,struc.const$time,struc.const$state)
  gam <- untrans(logistic(theta[3]),struc.const$age,struc.const$time,struc.const$state)
  p <- untrans(logistic(theta[4]),struc.const$age,list(1:(Time)),struc.const$state)
  
  psi <- make.psi(delt,kap,rho,gam)
  
  ll <- sl(phi,psi,p,mv)
  
  return(ll)
}