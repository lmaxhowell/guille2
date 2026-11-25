library(ggplot2)

# every datafile referenced here contains two dataframes.
# df is the results from the simulation
# df.true is the true values used for the datasets
# every simulation generates 445 datasets and optimises for the mles

# df has the following columns 
# par: the parameter name e.g. phi,p6
# mle: the optimised mle on the probability scale i.e. 0<=mle<=1
# convergence: whether the corresponding optimisation converged. In this case as optim is used, 0 means yes and anything else is a variation on no
# lower: the lower value of the 95% confidence interval, using the hessian, on the real scale. i.e. -infinty<lower<infinity
# mle.un: the optimised mle on the real scale i.e. -infinty<mle.un<infinity
# upper: the upper value of the 95% confidence interval, using the hessian, on the real scale. i.e. -infinty<upper<infinity
# var: the variance used in the 95% coming from the hessian
# true: the true value on the real scale i.e. -infinty<true<infinity
# n: a numbering of the datasets

# df.true has the following columns
# par: the parameter name e.g. phi,p6
# mle: the true value used to generate the dataset on the probability scale i.e. 0<=mle<=1
# note that mle is used here so it lines up with df to make plotting easier but
# it is NOT an mle it is the true value

# sim1:
#   survival age classes:2
#   split between years 1:2, and 3+
#   p time classes:2
#   split between years 1:6, and 7:16
#   everything else constant
load("sim1.RData")
ggplot(df,aes(par,mle)) + geom_boxplot() + geom_point(data=df.true,aes(par,mle),col="darkorange")

# sim2:
#   constant model
load("sim2.RData")
ggplot(df,aes(par,mle)) + geom_boxplot() + geom_point(data=df.true,aes(par,mle),col="darkorange")


# sim3:
#   survival constant
#   p time classes:2
#   split between years 1:6, and 7:16
#   everything else constant
load("sim3.RData")
ggplot(df,aes(par,mle)) + geom_boxplot() + geom_point(data=df.true,aes(par,mle),col="darkorange")


# sim4:
#   survival age classes:2
#   split between years 1, and 2+
#   p constant
#   everything else constant
load("sim4.RData")
ggplot(df,aes(par,mle)) + geom_boxplot() + geom_point(data=df.true,aes(par,mle),col="darkorange")


# sim5:
#   survival age classes:2
#   split between years 1:2, and 3+
#   p constant
#   everything else constant
load("sim5.RData")
ggplot(df,aes(par,mle)) + geom_boxplot() + geom_point(data=df.true,aes(par,mle),col="darkorange")


# sim6:
#   survival age classes:2
#   interates over where the two classes are split
#   p constant
#   everything else constant
load("sim6.RData")
ggplot(df,aes(par,mle)) + geom_boxplot() + geom_point(data=df.true,aes(par,mle),col="darkorange") + facet_wrap(~phi.split)


# sim7:
#   survival age classes:3
#   split between years 1,2, and 3+
#   p constant
#   everything else constant
load("sim7.RData")
ggplot(df,aes(par,mle)) + geom_boxplot() + geom_point(data=df.true,aes(par,mle),col="darkorange")

# sim8:
#   survival age classes:1
#   survival time classes:15
#   p constant
#   everything else constant
load("sim8.RData")
ggplot(df,aes(par,mle)) + geom_boxplot() + geom_point(data=df.true,aes(par,mle),col="darkorange")
