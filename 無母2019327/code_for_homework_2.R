
library("np")
data("cps71")
attach(cps71)

y <- logwage
x <- age
x2 <- age^2

n <- length(y)

########## OLS ##########

ls <- lm(y~x+x2)
summary(ls)
fit <- fitted(ls)

##########

plot(x,fit,type="l",col="black",xlab="Age",ylab="Log(Wage)",main="OLS",ylim=c(11.8,14.2),lwd=2)
points(x,y)
#################### codes for questions 1 and 5 ####################

#################### step 1 ####################

########## ROT ##########

bw.x <- 1.059*sd(x)*n^(-1/5)
bw.x

########## LSCV-LCLS ##########

h <- seq(1.5,2,0.001)

lscv <- matrix(0,nrow=length(h),ncol=1)

for (jj in 1:length(h)){

mhat <- matrix(0,nrow=n,ncol=1)

for (j in 1:n){
      
dx <- (x-x[j])/h[jj]
KK <- solve(sqrt(2*pi))*exp(-0.5*dx^2)

KK[j] <- 0	#leave one out
mhat[j] <- sum(y*KK)/sum(KK)

}

lscv[jj] <- sum((y-mhat)^2)

}

min.lscv <- which(lscv==min(lscv))
bw.x <- 1.5+0.001*min.lscv
bw.x

########## LSCV-LLLS ##########

h <- seq(3,3.5,0.001)

lscv <- matrix(0,nrow=length(h),ncol=1)

for (jj in 1:length(h)){

ones <- rep(1,n)
deltahat <- matrix(0,nrow=n,ncol=2) #一個mhat一個gradiant

for (j in 1:n){
      
dx <- (x-x[j])/h[jj]
KK <- solve(sqrt(2*pi))*exp(-0.5*dx^2)

KK[j] <- 0 #leave one out

xx <- cbind(ones,(dx*h[jj]))#X矩陣
xk <- t(xx*KK) #照片左半部，一種處理技巧

deltahat[j,] <- solve(xk%*%xx)%*%(xk%*%y)

}

mhat <- deltahat[,1]

lscv[jj] <- sum((y-mhat)^2)

}

min.lscv <- which(lscv==min(lscv))
bw.x <- 3+0.001*min.lscv
bw.x

#################### step 2 ####################

########## LCLS ##########

mhat <- matrix(0,nrow=n,ncol=1)

grad <- matrix(0,nrow=n,ncol=1)

for (j in 1:n){

dx <- (x-x[j])/bw.x
KK <- (1/(sqrt(2*pi)))*exp(-0.5*dx^2)
mhat[j] <- sum(y*KK)/sum(KK)

grad[j] <- ((sum(y*KK*(dx/bw.x))*sum(KK))-(sum(y*KK)*sum(KK*(dx/bw.x))))/((sum(KK))^2)

}

##########

plot(x,mhat,type="l",col="black",xlab="Age",ylab="Log(Wage)",main="LCLS ()",ylim=c(11.8,14.2),lwd=2)

plot(x,grad,type="l",col="black",xlab="Age",ylab="Marginal Effects",main="LCLS () and CI (Wild Bootstrap)",ylim=c(-0.2,0.3),lwd=2)
#上面那張圖的斜率變化在經濟直覺上就是邊際效果(OLS中負協率)
########## LLLS ##########

ones <- rep(1,n)
deltahat <- matrix(0,nrow=n,ncol=2)

for (j in 1:n){

dx <- (x-x[j])/bw.x
KK <- (1/(sqrt(2*pi)))*exp(-0.5*dx^2)

xx <- cbind(ones,(dx*bw.x))
xk <- t(xx*KK)

deltahat[j,] <- solve(xk%*%xx)%*%(xk%*%y)

}

mhat <- deltahat[,1]

grad <- deltahat[,2] 

##########

plot(x,mhat,type="l",col="black",xlab="Age",ylab="Log(Wage)",main="LLLS ()",ylim=c(11.8,14.2),lwd=2)

plot(x,grad,type="l",col="black",xlab="Age",ylab="Marginal Effects",main="LLLS () and CI (Wild Bootstrap)",ylim=c(-0.2,0.3),lwd=2)

#################### step 3 ####################

zero <- array(0,dim=c(length(x),1))

########## LCLS ##########

resid <- y-mhat

nb <- 200

boot.mhat <- matrix(0,nrow=n,ncol=nb)

boot.grad <- matrix(0,nrow=n,ncol=nb)

for (jj in 1:nb){

unif <- runif(length(resid),0,1)

a <- (1-sqrt(5))/2
prob.a <- (sqrt(5)+1)/(2*sqrt(5))
b <- (1+sqrt(5))/2

resid.star <- ifelse(unif<=prob.a,a*resid,b*resid)
y.star <- mhat+resid.star

for (j in 1:n){

dx <- (x-x[j])/bw.x
KK <- (1/(sqrt(2*pi)))*exp(-0.5*dx^2)

boot.mhat[j,jj] <- sum(y.star*KK)/sum(KK)

boot.grad[j,jj] <- ((sum(y.star*KK*(dx/bw.x))*sum(KK))-(sum(y.star*KK)*sum(KK*(dx/bw.x))))/((sum(KK))^2)

}

}

std.error <- sapply(1:n,function(i){sd(boot.mhat[i,])})

mhat.lower <- mhat-1.96*std.error
mhat.upper <- mhat+1.96*std.error

std.error <- sapply(1:n,function(i){sd(boot.grad[i,])})

grad.lower <- grad-1.96*std.error
grad.upper <- grad+1.96*std.error

##########

plot(x,mhat,type="l",col="black",xlab="Age",ylab="Log(Wage)",main="LCLS () and CI (Wild Bootstrap)",ylim=c(11.8,14.2),lwd=2)
lines(x,mhat.lower,lty=5,col="black",lwd=2)
lines(x,mhat.upper,lty=5,col="black",lwd=2)
legend("bottomright",c("Local Constant Estimates","95% Confidence Intervals"),col=c("black","black"),lty=c(1,5),bty="n")

plot(x,grad,type="l",col="black",xlab="Age",ylab="Marginal Effects",main="LCLS () and CI (Wild Bootstrap)",ylim=c(-0.2,0.3),lwd=2)
lines(x,grad.lower,lty=5,col="black",lwd=2)
lines(x,grad.upper,lty=5,col="black",lwd=2)
lines(x,zero,lty=1,col="red",lwd=2)
legend("topright",c("Gradients","95% Confidence Intervals"),col=c("black","black"),lty=c(1,5),bty="n")
#包道0為不顯著，代表年輕的時候表較有顯著的薪資增長經濟意義為年輕時因累積經驗所以加薪快快
########## LLLS ##########

resid <- y-mhat

nb <- 200

boot.mhat <- matrix(0,nrow=n,ncol=nb)

boot.grad <- matrix(0,nrow=n,ncol=nb) 

for (jj in 1:nb){

unif <- runif(length(resid),0,1)

a <- (1-sqrt(5))/2
prob.a <- (sqrt(5)+1)/(2*sqrt(5))
b <- (1+sqrt(5))/2

resid.star <- ifelse(unif<=prob.a,a*resid,b*resid)
y.star <- mhat+resid.star

ones <- rep(1,n)
deltahat <- matrix(0,nrow=n,ncol=2)

for (j in 1:n){

dx <- (x-x[j])/bw.x
KK <- (1/(sqrt(2*pi)))*exp(-0.5*dx^2)

xx <- cbind(ones,(dx*bw.x))
xk <- t(xx*KK)

deltahat[j,] <- solve(xk%*%xx)%*%(xk%*%y.star)

boot.mhat[j,jj] <- deltahat[j,1]

boot.grad[j,jj] <- deltahat[j,2]

}

}

std.error <- sapply(1:n,function(i){sd(boot.mhat[i,])})

mhat.lower <- mhat-1.96*std.error
mhat.upper <- mhat+1.96*std.error

std.error <- sapply(1:n,function(i){sd(boot.grad[i,])})

grad.lower <- grad-1.96*std.error
grad.upper <- grad+1.96*std.error

##########

plot(x,mhat,type="l",col="black",xlab="Age",ylab="Log(Wage)",main="LLLS () and CI (Wild Bootstrap)",ylim=c(11.8,14.2),lwd=2)
lines(x,mhat.lower,lty=5,col="black",lwd=2)
lines(x,mhat.upper,lty=5,col="black",lwd=2)
legend("bottomright",c("Local Linear Estimates","95% Confidence Intervals"),col=c("black","black"),lty=c(1,5),bty="n")

plot(x,grad,type="l",col="black",xlab="Age",ylab="Marginal Effects",main="LLLS () and CI (Wild Bootstrap)",ylim=c(-0.2,0.3),lwd=2)
lines(x,grad.lower,lty=5,col="black",lwd=2)
lines(x,grad.upper,lty=5,col="black",lwd=2)
lines(x,zero,lty=1,col="red",lwd=2)
legend("topright",c("Gradients","95% Confidence Intervals"),col=c("black","black"),lty=c(1,5),bty="n")

#################### answers to questions 2, 3, and 4 ####################

# 2 The dip is present in all of the resulting nonparametric estimates, and Pagan and Ullah (1999) provide an explanation for the dip in the chapter 3.14.2. #

# 3 Based on changes in confidence intervals, the dip appears to be significant. #

# 4 The local linear estimator appears to provide the most appropriate fit to this data, regardless of using the rule-of-thumb or cross-validated bandwidth, #
#   because it has no boundary effects and a smoother fitted curve. #

