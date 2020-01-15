
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

#-----------------------------------------

T<-100

c1<--1
c2<-1
ny1<-1
ny2<-1

N<-100
N1<-N/2
N2<-N/2

p1<-N1/N
p2<-N2/N

eta1<-3 #3
eta2<-2 #2
kappa<-eta1+eta2+2

#-----------------------------------------

pdf(width=12,height=3.3,"Fig_phase_EM_SP_Sub1_10_6.pdf")
par(mfrow=c(1,3))
par(mar=c(0.3,1, 1,3), oma=c(4,1.5,2,0), mai = c(0.1, 0.35, 0.1, 0.1))

#Euler-Maruyama: h1
set.seed(1)
h<-10^-2
grid<-seq(0,T,h)
#startv1<-rep(0,kappa)
startv1<-c(0,0,-3.588396,-3.975271,0,1.312606,1.081296)
sol<-Hawkes_EM_Cpp_(grid,h,startv1,c1,c2,ny1,ny2,p1,p2,N,eta1,eta2)
plot(sol[1,],sol[5,],type="l",col="black",xlab="",ylab="",xlim=c(-11,1),ylim=c(0,3),xaxt="n")
mtext(expression(paste(X^"2,1","(",t,")")), line = 2, side = 2, outer = F,cex=1.2)
legend("topright", legend=c("EM",expression(Delta==0.01)),
       col=c("red","green"), cex=1.7,bg = "white",lty=NULL,bty="n")

#Euler-Maruyama: h2
set.seed(1)
h<-5*10^-1
grid<-seq(0,T,h)
#startv1<-rep(0,kappa)
startv1<-c(0,0,-3.588396,-3.975271,0,1.312606,1.081296)
sol<-Hawkes_EM_Cpp_(grid,h,startv1,c1,c2,ny1,ny2,p1,p2,N,eta1,eta2)
plot(sol[1,],sol[5,],type="l",col="black",xlab="",ylab="",xlim=c(-11,1),ylim=c(0,3),xaxt="n",yaxt="n")
legend("topright", legend=c("EM",expression(Delta==0.5)),
       col=c("red","green"), cex=1.7,bg = "white",lty=NULL,bty="n")

#Euler-Maruyama: h3
set.seed(1)
h<-7*10^-1
grid<-seq(0,T,h)
#startv1<-rep(0,kappa)
startv1<-c(0,0,-3.588396,-3.975271,0,1.312606,1.081296)
sol<-Hawkes_EM_Cpp_(grid,h,startv1,c1,c2,ny1,ny2,p1,p2,N,eta1,eta2)
plot(sol[1,],sol[5,],type="l",col="black",xlab="",ylab="",xlim=c(-11,1),ylim=c(0,3),xaxt="n",yaxt="n")
legend("topright", legend=c("EM",expression(Delta==0.7)),
       col=c("red","green"), cex=1.7,bg = "white",lty=NULL,bty="n")

dev.off()

#--------------------------

pdf(width=12,height=3.3,"Fig_phase_EM_SP_Sub2_10_6.pdf")
par(mfrow=c(1,3))
par(mar=c(0.3,1, 1,3), oma=c(4,1.5,2,0), mai = c(0.1, 0.35, 0.1, 0.1))

#Splitting LT: h1
set.seed(1)
h<-10^-2
grid<-seq(0,T,h)
expmatA<-expmatA_SP(h,ny1,ny2,eta1,eta2)
startv1<-c(0,0,-3.588396,-3.975271,0,1.312606,1.081296)
sol<-Hawkes_SP_LT_Cpp_(grid,h,startv1,c1,c2,ny1,ny2,p1,p2,N,expmatA,eta1,eta2)
plot(sol[1,],sol[5,],type="l",col="black",xlab="",ylab="",xlim=c(-11,1),ylim=c(0,3),xaxt="n")
mtext(expression(paste(X^"2,1","(",t,")")), line = 2, side = 2, outer = F,cex=1.2)
legend("topright", legend=c("LT",expression(Delta==0.01)),
       col=c("red","green"), cex=1.7,bg = "white",lty=NULL,bty="n")

#Splitting LT: h2
set.seed(1)
h<-5*10^-1
grid<-seq(0,T,h)
expmatA<-expmatA_SP(h,ny1,ny2,eta1,eta2)
startv1<-c(0,0,-3.588396,-3.975271,0,1.312606,1.081296)
sol<-Hawkes_SP_LT_Cpp_(grid,h,startv1,c1,c2,ny1,ny2,p1,p2,N,expmatA,eta1,eta2)
plot(sol[1,],sol[5,],type="l",col="black",xlab="",ylab="",xlim=c(-11,1),ylim=c(0,3),xaxt="n",yaxt="n")
legend("topright", legend=c("LT",expression(Delta==0.5)),
       col=c("red","green"), cex=1.7,bg = "white",lty=NULL,bty="n")

#Splitting LT: h3
set.seed(1)
h<-7*10^-1
grid<-seq(0,T,h)
expmatA<-expmatA_SP(h,ny1,ny2,eta1,eta2)
startv1<-c(0,0,-3.588396,-3.975271,0,1.312606,1.081296)
sol<-Hawkes_SP_LT_Cpp_(grid,h,startv1,c1,c2,ny1,ny2,p1,p2,N,expmatA,eta1,eta2)
plot(sol[1,],sol[5,],type="l",col="black",xlab="",ylab="",xlim=c(-11,1),ylim=c(0,3),xaxt="n",yaxt="n")
legend("topright", legend=c("LT",expression(Delta==0.7)),
       col=c("red","green"), cex=1.7,bg = "white",lty=NULL,bty="n")

dev.off()

#--------------------------

pdf(width=12,height=3.3,"Fig_phase_EM_SP_Sub3_10_6.pdf")
par(mfrow=c(1,3))
par(mar=c(0.3,1, 1,3), oma=c(4,1.5,2,0), mai = c(0.1, 0.35, 0.1, 0.1))

#Splitting ST: h1
set.seed(1)
h<-10^-2
grid<-seq(0,T,h)
expmatA<-expmatA_SP(h/2,ny1,ny2,eta1,eta2)
startv1<-c(0,0,-3.588396,-3.975271,0,1.312606,1.081296)
sol<-Hawkes_SP_ST_Cpp_(grid,h,startv1,c1,c2,ny1,ny2,p1,p2,N,expmatA,eta1,eta2)
plot(sol[1,],sol[5,],type="l",col="black",xlab="",ylab="",xlim=c(-11,1),ylim=c(0,3))
mtext(expression(paste(X^"1,1","(",t,")")), line = 3, side = 1, outer = F,cex=1.2)
mtext(expression(paste(X^"2,1","(",t,")")), line = 2, side = 2, outer = F,cex=1.2)
legend("topright", legend=c("ST",expression(Delta==0.01)),
       col=c("red","green"), cex=1.7,bg = "white",lty=NULL,bty="n")

#Splitting ST: h2
set.seed(1)
h<-5*10^-1
grid<-seq(0,T,h)
expmatA<-expmatA_SP(h/2,ny1,ny2,eta1,eta2)
startv1<-c(0,0,-3.588396,-3.975271,0,1.312606,1.081296)
sol<-Hawkes_SP_ST_Cpp_(grid,h,startv1,c1,c2,ny1,ny2,p1,p2,N,expmatA,eta1,eta2)
plot(sol[1,],sol[5,],type="l",col="black",xlab="",ylab="",xlim=c(-11,1),ylim=c(0,3),yaxt="n")
mtext(expression(paste(X^"1,1","(",t,")")), line = 3, side = 1, outer = F,cex=1.2)
legend("topright", legend=c("ST",expression(Delta==0.5)),
       col=c("red","green"), cex=1.7,bg = "white",lty=NULL,bty="n")

#Splitting ST: h3
set.seed(1)
h<-7*10^-1
grid<-seq(0,T,h)
expmatA<-expmatA_SP(h/2,ny1,ny2,eta1,eta2)
startv1<-c(0,0,-3.588396,-3.975271,0,1.312606,1.081296)
sol<-Hawkes_SP_ST_Cpp_(grid,h,startv1,c1,c2,ny1,ny2,p1,p2,N,expmatA,eta1,eta2)
plot(sol[1,],sol[5,],type="l",col="black",xlab="",ylab="",xlim=c(-11,1),ylim=c(0,3),yaxt="n")
mtext(expression(paste(X^"1,1","(",t,")")), line = 3, side = 1, outer = F,cex=1.2)
legend("topright", legend=c("ST",expression(Delta==0.7)),
       col=c("red","green"), cex=1.7,bg = "white",lty=NULL,bty="n")

dev.off()

