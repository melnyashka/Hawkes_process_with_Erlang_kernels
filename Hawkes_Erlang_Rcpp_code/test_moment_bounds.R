
#----------------------------------------------------------------------
# Load Euler-Maruyama and Splitting
#----------------------------------------------------------------------

#Rcpp
library(Rcpp)
library(RcppNumerical)
library(devtools)
find_rtools(T)
sourceCpp(file="Hawkes_EM_Cpp.cpp")
sourceCpp(file="Hawkes_Splitting_Cpp.cpp")
source(file="Hawkes_Matrices.R")

#----------------------------------------------------------------------

c1<--1
c2<-1
ny1<-2 #changed it from 1 to 2!
ny2<-1

N<-20
N1<-N/2
N2<-N/2

p1<-N1/N
p2<-N2/N

eta1<-3
eta2<-2
kappa<-eta1+eta1+2
f2max<-40

#-----------------------------------------
# first moment
#-----------------------------------------

asymptotic_bound_Population1_first_moment <- function(eta1,f2max,c1,ny1,j,N,p2){
  bound<-(c1*f2max)/(ny1^(eta1+2-j))
  return(bound)
}

#-----------------------------------------
# second moment
#-----------------------------------------

asymptotic_bound_Population1_second_moment <- function(eta1,f2max,c1,ny1,j,N,p2){
  cinf<-factorial(2*(eta1+1-j))/( (2^(2*(eta1+1-j)+1))*(ny1^(2*(eta1+1-j)+1)))
  term1<-c1^2*(f2max)^2*(1/(ny1^(2*(eta1+2-j))))
  term2<-(1/N)*(c1^2/p2)*f2max*(1/(factorial(eta1+1-j)^2))*cinf
  term3<-2*((term1*term2)^(1/2))
  bound<-term1+term2+term3
  return(bound)
}

#-----------------------------------------

T<-50

#Splitting: ST
set.seed(1)
h<-10^-3
grid<-seq(0,T,h)
expmatA<-expmatA_SP(h/2,ny1,ny2,eta1,eta2)
startv1<-rep(0,kappa)
sol<-Hawkes_SP_ST_Cpp_(grid,h,startv1,c1,c2,ny1,ny2,p1,p2,N,expmatA,eta1,eta2)

#-----------------------------------------

#first moment
b1<-asymptotic_bound_Population1_first_moment(eta1,f2max,c1,ny1,1,N,p2)
b2<-asymptotic_bound_Population1_first_moment(eta1,f2max,c1,ny1,2,N,p2)
b3<-asymptotic_bound_Population1_first_moment(eta1,f2max,c1,ny1,3,N,p2)
b4<-asymptotic_bound_Population1_first_moment(eta1,f2max,c1,ny1,4,N,p2)

#do a quick check
b1
b2
b3
b4

mean(sol[1,])
mean(sol[2,])
mean(sol[3,])
mean(sol[4,])

#second moment
B1<-asymptotic_bound_Population1_second_moment(eta1,f2max,c1,ny1,1,N,p2)
B2<-asymptotic_bound_Population1_second_moment(eta1,f2max,c1,ny1,2,N,p2)
B3<-asymptotic_bound_Population1_second_moment(eta1,f2max,c1,ny1,3,N,p2)
B4<-asymptotic_bound_Population1_second_moment(eta1,f2max,c1,ny1,4,N,p2)

#do a quick check
B1
B2
B3
B4

mean(sol[1,]^2)
mean(sol[2,]^2)
mean(sol[3,]^2)
mean(sol[4,]^2)

#-------------------

pdf(width=12,height=3.3,"Fig_Moment_Bounds.pdf")
par(mfrow=c(1,3))
par(mar=c(0.3,1, 1,3), oma=c(4,1.5,2,0), mai = c(0.1, 0.5, 0.1, 0.1))

plot(grid,sol[1,],type="l",col="black",ylim=c(-21,4),xlab="",ylab="")
abline(h=b1,col="black",lty=1)
lines(grid,sol[2,],type="l",col="blue",lty=2)
abline(h=b2,col="blue",lty=2)
lines(grid,sol[3,],type="l",col="red",lty=3,lwd=1)
abline(h=b3,col="red",lty=3,lwd=1)
lines(grid,sol[4,],type="l",col="grey",lty=1)
abline(h=b4,col="grey",lty=1)
abline(h=0,col="black",lty=1)
mtext(expression("t"), line = 3, side = 1, outer = F,cex=1.3)
mtext(expression(paste(X^{paste(1,",",j)},"(",t,")")), line = 2, side = 2, outer = F,cex=1.3)
legend("topright", legend=c("first moment bounds"),
       col=c("red","green"), cex=1.7,bg = "white",lty=NULL,bty="n")

plot(grid,sol[1,]^2,type="l",col="black",ylim=c(-1,130),xlab="",ylab="")
abline(h=B1,col="black")
lines(grid,sol[2,]^2,type="l",col="blue",lty=2)
abline(h=B2,col="blue",lty=2)
lines(grid,sol[3,]^2,type="l",col="red",lty=3)
abline(h=B3,col="red",lty=3)
#lines(grid,sol[4,]^2,type="l",col="grey")
#abline(h=B4,col="grey")
mtext(expression("t"), line = 3, side = 1, outer = F,cex=1.3)
mtext(expression(paste("(",X^{paste(1,",",j)},"(",t,"))")^2), line = 2, side = 2, outer = F,cex=1.3)
legend("topright", legend=c("second moment bounds"),
       col=c("red","green"), cex=1.7,bg = "white",lty=NULL,bty="n")

dev.off()

