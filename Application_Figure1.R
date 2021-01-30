####################################################
# This R program reproduces Figure 1 in application
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

#Test for marginal cluster-level treatment effect
marginal.cluster <- function(m_bar, CV, rho=0.01, alpha=0.05, power=0.8, delta_x=0.25, pi_x=0.5, method="z") {
  a <- alpha
  b <- 1-power
  z_a <- qnorm(1-a/2)
  z_b <- qnorm(1-b)
  sigma2_y = 1
  
  eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
  omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
  if (method=="z"){
    n <- (z_a+z_b)^2*omega_x/delta_x^2
    n.final <- ceiling(ceiling(n)/2)*2
  } else if (method=="t"){
    n <- 4
    for (i in 1:2000){
      try.power <- pt(qt(1-a/2, n-2), n-2, ncp=delta_x/sqrt(omega_x/n), lower.tail = F) + pt(qt(a/2, n-2), n-2, ncp=delta_x/sqrt(omega_x/n))
      if (try.power>=power) {break}
      n <- n+2
    }
    n.final <- n
  }
  return(n.final)
}

#Test for marginal individual-level treatment effect
marginal.ind <- function(m_bar, CV, rho=0.01, alpha=0.05, power=0.8, delta_z=0.33, pi_z=0.5) {
  a <- alpha
  b <- 1-power
  z_a <- qnorm(1-a/2)
  z_b <- qnorm(1-b)
  sigma2_y = 1
  
  eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
  omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))
  n <- (z_a+z_b)^2*omega_z/delta_z^2
  n.final <- ceiling(ceiling(n)/2)*2
  return(n.final)
}

#Interaction test
int.test <- function(m_bar, CV, rho=0.01, alpha=0.05, power=0.8, delta_xz=0.3, pi_x=0.5, pi_z=0.5) {
  a <- alpha
  b <- 1-power
  z_a <- qnorm(1-a/2)
  z_b <- qnorm(1-b)
  sigma2_y = 1
  
  eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
  omega_xz <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_x*(1-pi_x)*pi_z*(1-pi_z))
  n <- (z_a+z_b)^2*omega_xz/delta_xz^2
  n.final <- ceiling(ceiling(n)/2)*2
  return(n.final)
}

#Joint test
joint.test <- function(m_bar, CV, rho=0.01, alpha=0.05, power=0.8, delta_x=0.25, delta_z=0.33, pi_x=0.5, pi_z=0.5, method="chi2") {
  a <- alpha
  b <- 1-power
  z_a <- qnorm(1-a/2)
  z_b <- qnorm(1-b)
  sigma2_y = 1
  
  eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
  omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
  eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
  omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))
  
  if (method=="chi2"){
    n <- 2
    for (i in 1:2000){
      theta <- n*(delta_x^2/omega_x + delta_z^2/omega_z)
      try.power <- pchisq(qchisq(1-a, 2), 2, ncp = theta, lower.tail = F)
      if (try.power>=power) {break}
      n <- n+2
    }
    n.final <- n
  } else if (method=="mixed"){
    n <- 4
    for (j in 1:2000){
      set.seed(1234567)
      #Simulate the mixed distribution (CENTRAL) to identify rejection region bound
      f.distn <- rf(10000, 1, n-2)
      chisq.distn <- rchisq(10000, 1)
      mix.distn <- f.distn + chisq.distn
      crt.value <- quantile(mix.distn, 1-a)
      #Simulate the mixed distribution (NONCENTRAL VERSION) for power calculation
      nc.f.distn <- rf(10000, 1, n-2, ncp = n*delta_x^2/omega_x)
      nc.chisq.distn <- rchisq(10000, 1, ncp = n*delta_z^2/omega_z)
      nc.mix.distn <- nc.f.distn + nc.chisq.distn
      try.power <- mean(nc.mix.distn>crt.value)
      if (try.power>=power) {break}
      n <- n+2
    }
    n.final <- n
  }
  return(n.final)
}

#Intersection-union test
IU.test <- function(m_bar, CV, rho=0.01, alpha=0.05, power=0.8, delta_x=0.25, delta_z=0.33, pi_x=0.5, pi_z=0.5, method="z") {
  a <- alpha
  b <- 1-power
  z_a <- qnorm(1-a/2)
  z_b <- qnorm(1-b)
  sigma2_y = 1
  
  eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
  omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
  eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
  omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))
  
  if (method=="z"){
    n <- 2
    for (i in 1:2000){
      wmean.c <- sqrt(n)*delta_x/sqrt(omega_x)
      wmean.i <- sqrt(n)*delta_z/sqrt(omega_z)
      try.power <- pnorm(qnorm(1-a/2), mean = wmean.c, lower.tail = F)*pnorm(qnorm(1-a/2), mean = wmean.i, lower.tail = F) + pnorm(qnorm(1-a/2), mean = wmean.c, lower.tail = F)*pnorm(qnorm(a/2), mean = wmean.i) + pnorm(qnorm(a/2), mean = wmean.c)*pnorm(qnorm(1-a/2), mean = wmean.i, lower.tail = F) + pnorm(qnorm(a/2), mean = wmean.c)*pnorm(qnorm(a/2), mean = wmean.i)
      if (try.power>=power) {break}
      n <- n+2
    }
    n.final <- n
  } else if (method=="mixed"){
    n <- 4
    for (j in 1:2000){
      c.ncp <- sqrt(n)*delta_x/sqrt(omega_x)
      i.mean <- sqrt(n)*delta_z/sqrt(omega_z)
      try.power <- pt(qt(1-a/2, n-2), n-2, c.ncp, lower.tail = F)*pnorm(qnorm(1-a/2), mean = i.mean, lower.tail = F) + pt(qt(1-a/2, n-2), n-2, c.ncp, lower.tail = F)*pnorm(qnorm(a/2), mean = i.mean) + pt(qt(a/2, n-2), n-2, c.ncp)*pnorm(qnorm(1-a/2), mean = i.mean, lower.tail = F) + pt(qt(a/2, n-2), n-2, c.ncp)*pnorm(qnorm(a/2), mean = i.mean)
      if (try.power>=power) {break}
      n <- n+2
    }
    n.final <- n
  }
  return(n.final)
}


#############
# Panel (A1)
#############

Mset <- seq(10, 100, by=2)
v.marcluster <- Vectorize(marginal.cluster, c("m_bar"))
t.CV0 <- v.marcluster(m_bar=Mset, CV=0, method="t")
t.CV0.3 <- v.marcluster(m_bar=Mset, CV=0.3, method="t")
t.CV0.6 <- v.marcluster(m_bar=Mset, CV=0.6, method="t")
t.CV0.9 <- v.marcluster(m_bar=Mset, CV=0.9, method="t")

plot(Mset, t.CV0, ylim=c(4,64), las=1,
     xlab=expression(bar(m)),
     main="(A1)",
     ylab="", cex.lab=2.5, cex.axis=2.2, cex.main=2.6, cex=2.2, 
     type="l", lwd=4.4, col="#0F2080",lty=1)
mtext(expression(n[A1]),side=2,las=1,line=3,cex=1.7)
lines(Mset, t.CV0.3,type="l", lwd=4.4, col="#85C0F9",lty=2)
lines(Mset, t.CV0.6,type="l", lwd=4.4, col="#DDCC77",lty=4)
lines(Mset, t.CV0.9,type="l", lwd=4.4, col="#F5793A",lty=5)


#############
# Panel (A2)
#############

v.marind <- Vectorize(marginal.ind, c("m_bar"))
resCV0 <- v.marind(m_bar=Mset, CV=0)
resCV0.3 <- v.marind(m_bar=Mset, CV=0.3)
resCV0.6 <- v.marind(m_bar=Mset, CV=0.6)
resCV0.9 <- v.marind(m_bar=Mset, CV=0.9)

plot(Mset, resCV0, ylim=c(4,64), xlab=expression(bar(m)), las=1,
     main="(A2)",
     ylab="", cex.lab=2.5, cex.axis=2.2, cex.main=2.6, cex=2.2, 
     type="l", lwd=4.4, col="#0F2080",lty=1.2)
mtext(expression(n[A2]),side=2,las=1,line=3,cex=1.7)
lines(Mset, resCV0.3,type="l", lwd=4.4, col="#85C0F9",lty=2)
lines(Mset, resCV0.6,type="l", lwd=4.4, col="#DDCC77",lty=4)
lines(Mset, resCV0.9,type="l", lwd=4.4, col="#F5793A",lty=5)


#############
# Panel (B)
#############

v.int <- Vectorize(int.test, c("m_bar"))
resCV0 <- v.int(m_bar=Mset, CV=0)
resCV0.3 <- v.int(m_bar=Mset, CV=0.3)
resCV0.6 <- v.int(m_bar=Mset, CV=0.6)
resCV0.9 <- v.int(m_bar=Mset, CV=0.9)

plot(Mset, resCV0, ylim=c(14,140), xlab=expression(bar(m)), las=1,
     main="(B)",
     ylab="", cex.lab=2.5, cex.axis=2.2, cex.main=2.6, cex=2.2, 
     type="l", lwd=4.4, col="#0F2080",lty=1)
mtext(expression(n[B]),side=2,las=1,line=3,cex=1.7)
lines(Mset, resCV0.3,type="l", lwd=4.4, col="#85C0F9",lty=2)
lines(Mset, resCV0.6,type="l", lwd=4.4, col="#DDCC77",lty=4)
lines(Mset, resCV0.9,type="l", lwd=4.4, col="#F5793A",lty=5)
#range(c(resCV0, resCV0.3, resCV0.6, resCV0.9))



#############
# Panel (C)
#############

v.joint <- Vectorize(joint.test, c("m_bar"))
mixed.CV0 <- v.joint(m_bar=Mset, CV=0, method="mixed")
mixed.CV0.3 <- v.joint(m_bar=Mset, CV=0.3, method="mixed")
mixed.CV0.6 <- v.joint(m_bar=Mset, CV=0.6, method="mixed")
mixed.CV0.9 <- v.joint(m_bar=Mset, CV=0.9, method="mixed")

plot(Mset, mixed.CV0, ylim=c(4,64), las=1,
     xlab=expression(bar(m)),
     main="(C)",
     ylab="", cex.lab=2.5, cex.axis=2.2, cex.main=2.6, cex=2.2, 
     type="l", lwd=4.4, col="#0F2080",lty=1)
mtext(expression(n[C]),side=2,las=1,line=3,cex=1.7)
lines(Mset, mixed.CV0.3, type="l", lwd=4.4, col="#85C0F9",lty=2)
lines(Mset, mixed.CV0.6, type="l", lwd=4.4, col="#DDCC77",lty=4)
lines(Mset, mixed.CV0.9, type="l", lwd=4.4, col="#F5793A",lty=5)


#############
# Panel (D)
#############

v.IU <- Vectorize(IU.test, c("m_bar"))
mixed.CV0 <- v.IU(m_bar=Mset, CV=0, method="mixed")
mixed.CV0.3 <- v.IU(m_bar=Mset, CV=0.3, method="mixed")
mixed.CV0.6 <- v.IU(m_bar=Mset, CV=0.6, method="mixed")
mixed.CV0.9 <- v.IU(m_bar=Mset, CV=0.9, method="mixed")

plot(Mset, mixed.CV0, ylim=c(4,64), las=1,
     xlab=expression(bar(m)),
     main="(D)",
     ylab="", cex.lab=2.5, cex.axis=2.2, cex.main=2.6, cex=2.2, 
     type="l", lwd=4.4, col="#0F2080",lty=1)
mtext(expression(n[D]),side=2,las=1,line=3,cex=1.7)
lines(Mset, mixed.CV0.3, type="l", lwd=4.4, col="#85C0F9",lty=2)
lines(Mset, mixed.CV0.6, type="l", lwd=4.4, col="#DDCC77",lty=4)
lines(Mset, mixed.CV0.9, type="l", lwd=4.4, col="#F5793A",lty=5)


##############################
# Create for a blank position
##############################

plot(Mset, mixed.CV0, ylim=c(4,64), las=1,
     xlab=expression(bar(m)),
     main="(D) Intersection-Union Test",
     ylab=expression(n), cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex=1.2, 
     type="l", lwd=3, col="#0F2080",lty=1)
lines(Mset, mixed.CV0.3, type="l", lwd=3, col="#85C0F9",lty=2)
lines(Mset, mixed.CV0.6, type="l", lwd=3, col="#A95AA1",lty=4)
lines(Mset, mixed.CV0.9, type="l", lwd=3, col="#F5793A",lty=5)


#Legends
par(mai = c(0, 0, 0, 0))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
plot_colors=c("#0F2080", "#85C0F9", "#DDCC77", "#F5793A")
legend(x="top",inset=0,
       legend = c("CV = 0", "CV = 0.3", "CV = 0.6", "CV = 0.9"), 
       col=plot_colors, lwd=4.4, cex=2.5, lty=c(1,1,1,1),
       horiz = TRUE)

dev.off()



