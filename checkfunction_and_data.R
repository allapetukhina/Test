# data
rm(list=ls(all=TRUE))

# please change your working directory
#setwd("C:/...")

library(xts)
load("data.RData")

window_length = 360

mean_sd_window_comp = function(dat, w_l) {
  
  if (length(na.omit(dat)) - w_l < 0) {
    return(NA)
  }
  
  startdate   = index(head(dat[! is.na(dat)], 1))
  dat         = dat[paste0(startdate,"::")]
  
  enddates    = index(dat)[-(1:w_l)]
  mean_w      = xts(rep(NA, length(dat) - w_l), order.by=enddates)
  sd_w        = xts(rep(NA, length(dat) - w_l), order.by=enddates)
  skew_w      = xts(rep(NA, length(dat) - w_l), order.by=enddates)
  kurt_w        = xts(rep(NA, length(dat) - w_l), order.by=enddates)
  
  for (enddate in enddates) {
    enddate         = as.Date(enddate)  # R for loop breaks class
    startdate       = enddate - w_l
    daterange       = paste0(startdate, "::", enddate)
    mean_w[enddate] = mean(dat[daterange], na.rm=TRUE)
    sd_w[enddate]   = sd(dat[daterange], na.rm=TRUE)
    skew_w[enddate] = skewness(dat[daterange], na.rm = TRUE)
    kurt_w[enddate] = kurtosis(dat[daterange], na.rm = TRUE)
  }
  return(cbind(mean_w, sd_w,  skew_w,  kurt_w))
}


  dat = crypto_returns_selec_xts#[, 1]
  w_l = window_length
  startdate   = index(head(dat[! is.na(dat[,1])], 1))
  dat         = dat[paste0(startdate,"::")]
  enddates    = index(dat)[-(1:w_l)]
  data = list()
  for (enddate in enddates) {
    enddate         = as.Date(enddate)  # R for loop breaks class
    startdate       = enddate - w_l
    daterange       = paste0(startdate, "::", enddate)
    enddate_str = as.character(enddate)
    data[[enddate_str]]  = dat[daterange]
  }


mean_sd_list = list()
for (crypto in max_cryptos) {
  mean_sd_list[[crypto]] = mean_sd_window_comp(crypto_returns_selec_xts[, 
                                                                        crypto], window_length)
}

cc_list = list()
for (crypto in max_cryptos) {
  cc_list[[crypto]] = window_comp(crypto_returns_selec_xts[, 
                                                           crypto], window_length)
}

#Calculation of weights for rolling windows
for(enddate in enddates[1]){
  enddate         = as.Date(enddate)
  enddate_str = as.character(enddate)
data_w =  data[[enddate_str]] 
# vector of expected returns
mu_na<- colMeans(data_w, na.rm = T)
mu_na = as.matrix(mu_na)
mu = mu_na[! is.na(mu_na[,1])]

# the dimension

n<- length(mu)
n

# matrix of covariance
#Sigma<-var(data_w,use="complete.obs")  
Sigma<- cov(data_w[,! is.na(mu_na[,1])],use="pairwise.complete.obs")
Sigma

CheckPositive<-function(matrix)
{
  if (qr(matrix)$rank <n) stop("positive definite assumption not satisified: the matrix is not of the full rank")
  else print("go for it")
} 

CheckPositive(Sigma)


# alpha and beta
# alpha<- pnorm(-1)
# beta<- pnorm(-2)


############ actual code

qq<- as.matrix(Sigma[1,])/sqrt(Sigma[1,1])
qq

###########
# check linear independence

CheckIndependence<- function(mu)
{
  if (qr(cbind(mu,qq,rep(1,n)))$rank <3) stop("linear independence assumption not satisified")
  else print("go for it")
}

CheckIndependence(mu)

############
# other calculations


Q<- Sigma- qq%*%t(qq)
Q

# function calculating CoVaR for a given x
CoVaR <- function(x) {as.numeric(- x%*%mu + a*x%*%qq + b*sqrt(t(x)%*%Q%*%x))}
#CoVaR(data_w)
#CoVaR(c(1/6,1/3,1/2))

Qhat<- Q[2:n,2:n]
Qhat

QhatInv<- solve(Qhat)

muhat<- as.matrix((mu - mu[1]*rep(1,n))[2:n])
muhat

qqhat<- as.matrix((qq - qq[1]*rep(1,n))[2:n])
qqhat


alphaC<- as.numeric(t(muhat)%*%QhatInv%*%muhat)
betaC<- as.numeric(t(muhat)%*%QhatInv%*%qqhat)
gammaC<- as.numeric(t(qqhat)%*%QhatInv%*%qqhat)

alphaC
betaC
gammaC

# have everything checked for given significance levels
Check<- function(alpha,beta) {
  if (alpha>=1/2 | beta>=1/2) stop("significance levels assumption not satisfied!")
  a<- -qnorm(alpha)
  b<- - qnorm(beta)
  delta<- b^2*alphaC-a^2*(alphaC*gammaC-betaC^2)
  expr <- a*betaC- alphaC
  print(c("square root delta:", sqrt(delta)))
  print(c("a*betaC- alphaC:", expr))
  if (delta<= 0) print("CoVaR unbounded below")
  if (expr<= - sqrt(delta)) print("no efficient portoflios")
  if (- sqrt(delta)< expr & expr<= sqrt(delta)) print("critical portfolio efficient iff E>mu1")
  if (expr>sqrt(delta)) print("all critical portfolios efficient")
}

#Check(1/2,1/10)
#Check(1/4,1/10)


# have everything checked for the neqative probit function of the significance levels
# ( the negative quantiles)
Check2<- function(a,b) {
  alpha<- pnorm(-a)
  beta<- pnorm(-b)
  if (alpha>=1/2 | beta>=1/2) stop("significance levels assumption not satisfied!")
  delta<- b^2*alphaC-a^2*(alphaC*gammaC-betaC^2)
  expr <- a*betaC- alphaC
  print(c("square root delta:", sqrt(delta)))
  print(c("a*betaC- alphaC:", expr))
  if (delta<= 0) {  print("CoVaR unbounded below")}
  else if (expr<= - sqrt(delta)) {check = 1 
  print("no efficient portoflios")}
  else if (- sqrt(delta)< expr & expr<= sqrt(delta)) {check = 2 
  print("critical portfolio efficient iff E>mu1")}
  else if (expr>sqrt(delta)) {check = 3
  print("all critical portfolios efficient")}
  return(check)
}

# Check2(2,3/10)
# Check2(2,sqrt(3/23))
# Check2(2,1/2)
# Check2(4,3/sqrt(23))
# Check2(4,3/5)

# choose your quantiles 
alpha = 0.05
beta = 0.05

a<- -qnorm(alpha)
#a
b<- -qnorm(beta)
#b
delta<- b^2*alphaC-a^2*(alphaC*gammaC-betaC^2)
delta

Check(alpha, beta)
if(Check2(a, b)==1){
  R_T = seq(mu[1]-1,max(mu),length.out = 10)
}
else if(Check2(a, b)==2){
  R_T = seq(mu[1], max(mu), length.out = 10)
}
else if {
  R_T = seq(mu[1]-1,max(mu), length.out = 10)
}

R_T = seq(-1,1, length.out =  10)
R_T = as.matrix(R_T)

xhat <- function(E) { (E-mu[1])/alphaC *(QhatInv%*%muhat) + abs(E - mu[1])*a/(alphaC*sqrt(delta))*QhatInv%*%(betaC*muhat - alphaC*qqhat)}

apply(R_T, 1, xhat)

apply(R_T, 1, image)

# those are the weights in the portfolio; xhat only helps to create this
x<- function(E) {rbind(1 - sum(xhat(E)),xhat(E))}
x(R_T)

# this is the image of the critical portfolios; part of it can be the efficient frontier
image<- function(E) {-mu[1]+a*sqrt(Sigma[1,1])+(E-mu[1])*(a*betaC/alphaC-1) + abs(E-mu[1])/alphaC*sqrt(delta)}

plot(x=apply(R_T,1,image),y=R_T)

plot(R_T, image(R_T), type = "p")
}

