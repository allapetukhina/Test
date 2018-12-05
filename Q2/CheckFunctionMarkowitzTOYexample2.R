# data

# vector of expected returns
mu<- as.matrix(c(1,2,3))
mu

# the dimension

n<- length(mu)
n

# matrix of covariance

Sigma<- rbind(c(1,1,2),c(1,9,0),c(2,0,16))
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

#a<- -qnorm(alpha)
#a
#b<- -qnorm(beta)
#b

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
CoVaR(c(1/6,1/3,1/2))

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
Check(1/2,1/10)
Check(1/4,1/10)


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
if (delta<= 0) print("CoVaR unbounded below")
if (expr<= - sqrt(delta)) print("no efficient portoflios")
if (- sqrt(delta)< expr & expr<= sqrt(delta)) print("critical portfolio efficient iff E>mu1")
if (expr>sqrt(delta)) print("all critical portfolios efficient")
}

Check2(2,3/10)
Check2(2,sqrt(3/23))
Check2(2,1/2)
Check2(4,3/sqrt(23))
Check2(4,3/5)

# choose your quantiles 
a
b
delta<- b^2*alphaC-a^2*(alphaC*gammaC-betaC^2)
delta



xhat <- function(E) { (E-mu[1])/alphaC *(QhatInv%*%muhat) + abs(E - mu[1])*a/(alphaC*sqrt(delta))*QhatInv%*%(betaC*muhat - alphaC*qqhat)}
xhat(4)


# those are the weights in the portfolio; xhat only helps to create this
x<- function(E) {rbind(1 - sum(xhat(E)),xhat(E))}
x(4)

# this is the image of the critical portfolios; part of it can be the efficient frontier
image<- function(E) {-mu[1]+a*sqrt(Sigma[1,1])+(E-mu[1])*(a*betaC/alphaC-1) + abs(E-mu[1])/alphaC*sqrt(delta)}

