####################################################
# This R program reproduces Figure 1 in Application
####################################################

#Only create pictures for formulas with small-sample corrections

setwd("...") #need to change the directory here
pdf("Figure1.pdf", width=16, height=12, paper="special")

m=matrix(c(1,2,3,4,5,6,7,7,7), nrow=3, ncol=3, byrow=TRUE)
layout(mat=m, heights = c(5.5, 5.5, 1))
par(mar=c(5.1, 6.1, 4.1, 2.1))


#######################################################################################
# Functions to calculate required number of clusters under different hypothesis tests
#######################################################################################

main.cluster <- function(alpha=0.05, beta=0.2, sigma2y=1, eff=0.25, m_bar, rho=0.01, CV, pi_x=0.5, pi_z=0.5){
  
  nvar2 <- sigma2y*(1+(m_bar-1)*rho)/(pi_x*(1-pi_x)*m_bar) * 
    ( (1+(m_bar-1)*rho)^2/((1+(m_bar-1)*rho)^2-CV^2*m_bar*rho*(1-rho)) + pi_z*(1-rho)*(1+(m_bar-1)*rho)^2/((1-pi_z)*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho))) )
  
  n <- 2
  power <- 0
  while (power < 1-beta){
    n <- n+2
    power <- pt(qt(1-alpha/2, n-2), n-2, ncp=eff/sqrt(nvar2/n), lower.tail = F) + pt(qt(alpha/2, n-2), n-2, ncp=eff/sqrt(nvar2/n))
  }
  return(n)
}


main.ind <- function(alpha=0.05, beta=0.2, sigma2y=1, eff=0.33, m_bar, rho=0.01, CV, pi_x=0.5, pi_z=0.5){
  
  nvar3 <- sigma2y*(1-rho)*(1+(m_bar-1)*rho)^3/((1-pi_z)*pi_x*(1-pi_x)*m_bar*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho)))
  n <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*nvar3/eff^2
  n <- 2*(ceiling(n/2))
  return (n)
  
}


marginal.cluster <- function(alpha=0.05, beta=0.2, sigma2y=1, eff=0.35, m_bar, rho=0.01, CV, pi_x=0.5, pi_z=0.5){
  
  mbar_plus_eta2bar <- m_bar*(1-rho)/( 1+(m_bar-1)*rho )*(1-CV^2*m_bar*rho*(1-rho)/(1+(m_bar-1)*rho)^2)
  omega_x <- sigma2y*(1-rho)/( mbar_plus_eta2bar*pi_x*(1-pi_x) )
  
  n <- 2
  power <- 0
  while (power < 1-beta){
    n <- n+2
    power <- pt(qt(1-alpha/2, n-2), n-2, ncp=eff/sqrt(omega_x/n), lower.tail = F) + pt(qt(alpha/2, n-2), n-2, ncp=eff/sqrt(omega_x/n))
  }
  return(n)
}


marginal.ind <- function(alpha=0.05, beta=0.2, sigma2y=1, eff=0.43, m_bar, rho=0.01, CV, pi_x=0.5, pi_z=0.5){
  mbar_plus_eta1bar <- ( m_bar*(1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar^2*rho^2*(1-rho) )/(1+(m_bar-1)*rho)^3
  omega_z <- sigma2y*(1-rho)/( mbar_plus_eta1bar*pi_z*(1-pi_z) )
  
  n <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*omega_z/eff^2
  n <- 2*(ceiling(n/2))
  return (n)
}


interaction <- function(alpha=0.05, beta=0.2, sigma2y=1, eff=0.2, m_bar, rho=0.01, CV, pi_x=0.5, pi_z=0.5){
  mbar_plus_eta1bar <- ( m_bar*(1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar^2*rho^2*(1-rho) )/(1+(m_bar-1)*rho)^3
  omega_xz <- sigma2y*(1-rho)/( mbar_plus_eta1bar*pi_x*pi_z*(1-pi_x)*(1-pi_z) )
  
  n <- (qnorm(1-alpha/2) + qnorm(1-beta))^2*omega_xz/eff^2
  n <- 2*(ceiling(n/2))
  return (n)
}



#############
# Panel (A1)
#############

Mset <- seq(10, 100, by=2)
v.main.cluster <- Vectorize(main.cluster, c("m_bar"))
t.CV0 <- v.main.cluster(m_bar=Mset, CV=0)
t.CV0.3 <- v.main.cluster(m_bar=Mset, CV=0.3)
t.CV0.6 <- v.main.cluster(m_bar=Mset, CV=0.6)
t.CV0.9 <- v.main.cluster(m_bar=Mset, CV=0.9)

max(c(t.CV0, t.CV0.3, t.CV0.6, t.CV0.9)) #112
min(c(t.CV0, t.CV0.3, t.CV0.6, t.CV0.9)) #18

plot(Mset, t.CV0, ylim=c(6,112), las=1,
     xlab=expression(bar(m)),
     main="(A1) Cluster-level CE",
     ylab="", cex.lab=2.5, cex.axis=2.2, cex.main=2.6, cex=2.2, 
     type="l", lwd=4.4, col="#0F2080",lty=1)
mtext(expression(n[A1]),side=2,las=1,line=3,cex=1.7)
lines(Mset, t.CV0.3,type="l", lwd=4.4, col="#85C0F9",lty=2)
lines(Mset, t.CV0.6,type="l", lwd=4.4, col="#DDCC77",lty=4)
lines(Mset, t.CV0.9,type="l", lwd=4.4, col="#F5793A",lty=5)


#############
# Panel (A2)
#############

v.main.ind <- Vectorize(main.ind, c("m_bar"))
resCV0 <- v.main.ind(m_bar=Mset, CV=0)
resCV0.3 <- v.main.ind(m_bar=Mset, CV=0.3)
resCV0.6 <- v.main.ind(m_bar=Mset, CV=0.6)
resCV0.9 <- v.main.ind(m_bar=Mset, CV=0.9)

max(c(resCV0, resCV0.3, resCV0.6, resCV0.9))#58
min(c(resCV0, resCV0.3, resCV0.6, resCV0.9))#6

plot(Mset, resCV0, ylim=c(6,112), xlab=expression(bar(m)), las=1,
     main="(A2) Individual-level CE",
     ylab="", cex.lab=2.5, cex.axis=2.2, cex.main=2.6, cex=2.2, 
     type="l", lwd=4.4, col="#0F2080",lty=1.2)
mtext(expression(n[A2]),side=2,las=1,line=3,cex=1.7)
lines(Mset, resCV0.3,type="l", lwd=4.4, col="#85C0F9",lty=2)
lines(Mset, resCV0.6,type="l", lwd=4.4, col="#DDCC77",lty=4)
lines(Mset, resCV0.9,type="l", lwd=4.4, col="#F5793A",lty=5)


######### For Blank ################
plot(Mset, resCV0, ylim=c(2,112), xlab=expression(bar(m)), las=1,
     main="(A2)",
     ylab="", cex.lab=2.5, cex.axis=2.2, cex.main=2.6, cex=2.2, 
     type="l", lwd=4.4, col="#0F2080",lty=1.2)
mtext(expression(n[A2]),side=2,las=1,line=3,cex=1.7)
lines(Mset, resCV0.3,type="l", lwd=4.4, col="#85C0F9",lty=2)
lines(Mset, resCV0.6,type="l", lwd=4.4, col="#DDCC77",lty=4)
lines(Mset, resCV0.9,type="l", lwd=4.4, col="#F5793A",lty=5)

#############
# Panel (B1)
#############

v.marginal.cluster <- Vectorize(marginal.cluster, c("m_bar"))
t.CV0 <- v.marginal.cluster(m_bar=Mset, CV=0)
t.CV0.3 <- v.marginal.cluster(m_bar=Mset, CV=0.3)
t.CV0.6 <- v.marginal.cluster(m_bar=Mset, CV=0.6)
t.CV0.9 <- v.marginal.cluster(m_bar=Mset, CV=0.9)

max(c(t.CV0, t.CV0.3, t.CV0.6, t.CV0.9)) #32
min(c(t.CV0, t.CV0.3, t.CV0.6, t.CV0.9)) #8

plot(Mset, t.CV0, ylim=c(2,32), las=1,
     xlab=expression(bar(m)),
     main="(B1) Cluster-level NE",
     ylab="", cex.lab=2.5, cex.axis=2.2, cex.main=2.6, cex=2.2, 
     type="l", lwd=4.4, col="#0F2080",lty=1)
mtext(expression(n[B1]),side=2,las=1,line=3,cex=1.7)
lines(Mset, t.CV0.3,type="l", lwd=4.4, col="#85C0F9",lty=2)
lines(Mset, t.CV0.6,type="l", lwd=4.4, col="#DDCC77",lty=4)
lines(Mset, t.CV0.9,type="l", lwd=4.4, col="#F5793A",lty=5)


#############
# Panel (B2)
#############

v.marginal.ind <- Vectorize(marginal.ind, c("m_bar"))
resCV0 <- v.marginal.ind(m_bar=Mset, CV=0)
resCV0.3 <- v.marginal.ind(m_bar=Mset, CV=0.3)
resCV0.6 <- v.marginal.ind(m_bar=Mset, CV=0.6)
resCV0.9 <- v.marginal.ind(m_bar=Mset, CV=0.9)

max(c(resCV0, resCV0.3, resCV0.6, resCV0.9))#18
min(c(resCV0, resCV0.3, resCV0.6, resCV0.9))#2

plot(Mset, resCV0, ylim=c(2,32), xlab=expression(bar(m)), las=1,
     main="(B2) Individual-level NE",
     ylab="", cex.lab=2.5, cex.axis=2.2, cex.main=2.6, cex=2.2, 
     type="l", lwd=4.4, col="#0F2080",lty=1.2)
mtext(expression(n[B2]),side=2,las=1,line=3,cex=1.7)
lines(Mset, resCV0.3,type="l", lwd=4.4, col="#85C0F9",lty=2)
lines(Mset, resCV0.6,type="l", lwd=4.4, col="#DDCC77",lty=4)
lines(Mset, resCV0.9,type="l", lwd=4.4, col="#F5793A",lty=5)

#############
# Panel (C)
#############

v.int <- Vectorize(interaction, c("m_bar"))
resCV0 <- v.int(m_bar=Mset, CV=0)
resCV0.3 <- v.int(m_bar=Mset, CV=0.3)
resCV0.6 <- v.int(m_bar=Mset, CV=0.6)
resCV0.9 <- v.int(m_bar=Mset, CV=0.9)


plot(Mset, resCV0, ylim=c(32,314), xlab=expression(bar(m)), las=1,
     main="(C) Interaction",
     ylab="", cex.lab=2.5, cex.axis=2.2, cex.main=2.6, cex=2.2, 
     type="l", lwd=4.4, col="#0F2080",lty=1)
mtext(expression(n[C]),side=2,las=1,line=3,cex=1.7)
lines(Mset, resCV0.3,type="l", lwd=4.4, col="#85C0F9",lty=2)
lines(Mset, resCV0.6,type="l", lwd=4.4, col="#DDCC77",lty=4)
lines(Mset, resCV0.9,type="l", lwd=4.4, col="#F5793A",lty=5)
range(c(resCV0, resCV0.3, resCV0.6, resCV0.9))


#Legends
par(mai = c(0, 0, 0, 0))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
plot_colors=c("#0F2080", "#85C0F9", "#DDCC77", "#F5793A")
legend(x="top",inset=0,
       legend = c("CV = 0", "CV = 0.3", "CV = 0.6", "CV = 0.9"), 
       col=plot_colors, lwd=4.4, cex=2.5, lty=c(1,1,1,1),
       horiz = TRUE)

dev.off()



