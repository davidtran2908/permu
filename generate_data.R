## SIMULATION STUDY - DATA GENERATION FOR TEST POWER

## setup
rm(list=ls())
set.seed(1)
setwd("C:/Users/tran23/Downloads/take_away/permu_project/")

## parameter setting
n <- 4
uisd<-1
tau <- c(0, 0.01, 0.1, 0.3, 0.5, 0.7, 1) 
tau2 <- tau^2
n.rep <- 1000
n.subject <- c("verysmall", "small", "medium")
n<-c(4, 10)

## simulation
df.ma.lst <- list(NULL)
for (n in c(4)) {
  if(n==4){
    beta0 <- c(0.9, 2.3, 0.3, 0.1) 
    beta1 <- c(0.8, 0.7, 0.9, 0.9)
  }
  if(n==10){
    beta0 <- c(0.9, 2.3, 0.3, 0.1, 0.9, 2.3, 0.3, 0.1, 0.9, 2.3) 
    beta1 <- c(0.8, 0.7, 0.9, 0.9, 0.8, 0.7, 0.9, 0.9, 0.8, 0.7)
  }
  for (j in 2:2) {
    if(n==4 & j==1){
      ni <- rbinom(n, 1, 0.8)
      ni <- round(ni*runif(n, min=15, max=30)+(1-ni)*runif(n, min=30, max=100))
    }
    if(n>4 & j==1) next
    if(j==2){
      ni <- round(runif(n, min=30, max=100))
    }
    if(j==3){
      #ni <- round(runif(n, min=30, max=300))
      ni <- round(runif(n, min=100, max=200))
    }
    #sigma <- sqrt(1/ni)
    #sigma <- uisd*rep(1, n)
    sigma <- c(0.94, 0.89, 0.86, 1.39)
    sigma2 <- sigma^2
    p.treat <- runif(n, min=0.5, max=0.7)
    xi <- y0i <- list(NULL)
    for (i in 1:n) {
      xi[[i]] <- rbinom(n=ni[i], size=1, prob=p.treat[i])
      y0i[[i]] <- rnorm(ni[i], mean=4, sd=1)
    }
    for (theta in c(0, 1)) {
      for (k in 1:length(tau)) {
        for (b in 1:n.rep) {
          df.ma <- df.trial <-NULL
          beta2 <- rnorm(n, mean=theta, sd=tau[k]) ## norm

          for (i in 1:n) {
            epsilon.i <- rnorm(ni[i], mean=0, sd=sigma[i]) 
            #epsilon.i <- (sigma[i]/sqrt(3))*rt(ni[i], df=3) ## student t
            #epsilon.i <- (sigma[i]/sqrt((exp(1)-1)*exp(1)))*(rlnorm(ni[i])-exp(1/2)) ## lognormal(0, 1)
            yi <- beta0[i] + beta1[i]*y0i[[i]] + beta2[i]*xi[[i]] + epsilon.i

            df.trial <- data.frame(Study = as.factor(rep(i, ni[i])), Group = xi[[i]], y0 = y0i[[i]], y = yi)
            df.ma <- rbind(df.ma, df.trial)
          } 
          df.ma.lst[[b]] <- df.ma 
          #save.image(paste("sim_dat_epsilonnorm_theta", theta, "_n", n, "_tau", tau[k], "_n.subject", n.subject[j], ".RData", sep=""))
          #save.image(paste("sim_dat_epsilonlognorm_theta", theta, "_n", n, "_tau", tau[k], "_n.subject", n.subject[j], ".RData", sep=""))
          save.image(paste("sim_dat_diffuisd_theta", theta, "_n", n, "_tau", tau[k], "_n.subject", n.subject[j], ".RData", sep=""))
          #save.image(paste("sim_dat_epsilont_theta", theta, "_n", n, "_tau", tau[k], "_n.subject", n.subject[j], ".RData", sep=""))
          print(b)
        }  
      }
    }
  }  
}
