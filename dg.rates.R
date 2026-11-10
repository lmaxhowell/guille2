dg.rates <- function(ch){
  tb <- table(unlist(ch))
  delta <- unname((tb["L_B"] + tb["L_B_"])/(tb["L_B"] + tb["L_B_"] + tb["LB"] + tb["LB_"])) # delta
  gamma <- unname((tb["L_B"] + tb["LB"])/(tb["L_B"] + tb["L_B_"] + tb["LB"] + tb["LB_"])) # gamma
  return(list("delta"=delta,"gamma"=gamma))
}