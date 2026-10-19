load("new_data_AUK.RData")
states <- c("N","B1","LB","L_B","LB_","L_B_","S")
Time <- ncol(new_data_AUK)-1

ch <- new_data_AUK[,1:16]
colnames(ch) <- c(2010:2025)
states_code <- c(1,"B1","LSB","DSB","LUB","DUB",13)

for(s in 1:length(states_code)){
  ch[which(ch==states_code[s],arr.ind=TRUE)] <- states[s]
}

for(i in 1:nrow(ch)){
  # want to add Nbar as a state between any Ns or until N becomes B1
  # and want to add S as a state between the breeding states
  wn <- which(ch[i,]=="N")
  wb1 <- which(ch[i,]=="B1")
  wnlss <- wn[wn<wb1] # which non breeding states are before the first breeding state (as opposed to the non breeding due to skipping)
  # if(length(wnlss)==0){
  #   if(length(wn)>=2){
  #     for(j in 1:(Time-1)){
  #       if(j>min(wn) & j<max(wn) & ch[i,j]=="0"){
  #         ch[i,j] <- "E"
  #       }
  #     }
  #   }
  # }else{
  #   if(length(wnlss)>0){
  #     for(j in 1:(Time-1)){
  #       if(j>min(wn) & j<wb1 & ch[i,j]=="0"){
  #         ch[i,j] <- "E"
  #       }
  #     }
  #   }
  # }
  wbb <- which(ch[i,] %in% c("B1","LB","L_B","LB_","L_B_"))
  if(length(wbb)>1){
    for(j in 1:(length(wbb)-1)){
      if((wbb[j+1]-wbb[j])>1){
        ch[i,(wbb[j]+1):(wbb[j+1]-1)] <- "S"
      }
    }
  }
}

ch <- cbind(ch,"age"=rep(1,nrow(ch)))

ni <- rep(0,16)
for(i in 1:nrow(ch)){
  ni[which(ch[i,1:16]!=0)[1]] <- ni[which(ch[i,1:16]!=0)[1]]+1
}

ni2 <- list()
for(i in 1:length(ni)){
  ni2[[i]] <- rep(1,ni[i])
}

ch2 <- ch
ch <- ch[-c(250,289,367,391,602,69,483,124,386,415,169,128),]

save(ch,ch2,ni,ni2,file="guillemot.RData")