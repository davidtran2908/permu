## SIMULATION STUDY - COMPUTE P-VALUE

## setup
rm(list=ls())
library(nlme)
#library(predictmeans)
#library(lmeresampler)
library(emmeans)
library(tictoc)
library(gtools)
library(parallel)
library(expm)
library(lme4)
library(lmerTest)
set.seed(1)
setwd("/mnt/vast-standard/home/tran23/u14326/permutation/new_system/") #C:/Users/tran23/Downloads/take_away/permu_project/typ_one_err/august_2024
source("/mnt/vast-standard/home/tran23/u14326/permlme2.R") #C:/Users/tran23/Downloads/take_away/permu_project/permlme2.R

## parameter setting
n <- 4
uisd<-1
tau <- c(0, 0.01, 0.1, 0.3, 0.5, 0.7, 1) 
tau2 <- tau^2
n.rep <- 1000
n.subject <- c("verysmall", "small", "medium")
n<-c(4)

for (n in c(4)) {
  for (j in 1:3) {
    for (theta in c(0, 1)) {
      for (k in c(2,3,4,5,7)) {
        ## load data
        #load(paste("sim_dat_epsilonnorm_theta", theta, "_n", n, "_tau", tau[k], "_n.subject", n.subject[j], ".RData", sep=""))
        #load(paste("sim_dat_epsilonlognorm_theta", theta, "_n", n, "_tau", tau[k], "_n.subject", n.subject[j], ".RData", sep=""))
        #load(paste("sim_dat_diffuisd_theta", theta, "_n", n, "_tau", tau[k], "_n.subject", n.subject[j], ".RData", sep=""))
        #load(paste("sim_dat_epsilont_theta", theta, "_n", n, "_tau", tau[k], "_n.subject", n.subject[j], ".RData", sep=""))
        
        ## compute p-value
        mod.n.df <- mod.s.df <- mod.kh.df  <- mod.b.df<- mod.kr.df <-  data.frame(matrix(NA, nrow=n.rep, ncol=6))
        mod.p.df <- data.frame(matrix(NA, nrow=n.rep, ncol=7))
        colnames(mod.n.df) <- colnames(mod.s.df) <- colnames(mod.kh.df) <- colnames(mod.b.df) <- colnames(mod.kr.df) <- c("est", "se", "97.5-perc", "p-value", "low", "hi")
        colnames(mod.p.df) <- c("est", "se", "97.5-perc", "2.5-perc", "p-value", "low", "hi")
        for (b in 1:n.rep) {
          #for (b2 in 1:n.rep) {
          #b<-ind.fail[b2]
          ## model fitting
          mod <- try(lmer(y ~ 0 + Study + y0:Study + Group+(0+Group|Study),
                          na.action = na.omit, data=df.ma.lst[[b]]), silent=TRUE)
          
          if(!(inherits(mod, "try-error"))){
            ## est
            # mod.tau2 <- VarCorr(mod)
            # mod.tau2 <- as.numeric(mod.tau2["Group", "Variance"])
            mod.beta <- summary(mod) ## beta
            mod.beta <- mod.beta$coefficients
            mod.n.df[b,1:2] <- mod.beta["Group", c("Estimate", "Std. Error")]
            mod.n.df[b,3] <- qnorm(0.975)
            mod.n.df[b,4] <- 2*pnorm(abs(mod.beta["Group", "t value"]), lower.tail = FALSE)
            
            ## satterthwaite
            mod.s <- summary(mod, ddf="Satterthwaite")
            mod.s.df[b,1:2] <- mod.s$coefficients["Group", c("Estimate", "Std. Error")]
            mod.s.dof<- mod.s$coefficients["Group", "df"]
            mod.s.df[b,3] <- qt(0.975, df=mod.s.dof)
            mod.s.df[b,4] <- mod.s$coefficients["Group", c("Pr(>|t|)")]
            
            ## kenward-roger
            mod.kr <- summary(mod, ddf="Kenward-Roger")
            mod.kr.df[b,1:2] <- mod.kr$coefficients["Group", c("Estimate", "Std. Error")]
            mod.kr.dof<- mod.kr$coefficients["Group", "df"]
            mod.kr.df[b,3] <- qt(0.975, df=mod.kr.dof)
            mod.kr.df[b,4] <- mod.kr$coefficients["Group", c("Pr(>|t|)")]
            
            ## knapp-hartung
            mod.kh.dof <- n-1
            mod.kh.pv <- 2*pt(abs(mod.beta["Group", "t value"]), df=mod.kh.dof, lower.tail = FALSE)
            mod.kh <- c(mod.beta["Group", c("Estimate", "Std. Error")], qt(0.975, df=mod.kh.dof), mod.kh.pv)
            mod.kh.df[b,-c(5,6)] <- mod.kh
            
            ## permutation
            mod.p <- mod.p2 <- mod.p3 <- rep(NA,4)
            mod0 <- try(lmer(y ~ 0 + Study + y0:Study + (0+Group|Study),
                             na.action = na.omit, data=df.ma.lst[[b]]), silent=TRUE)
            # mod1 <- try(lmer(y ~ 0 + Study + y0:Study + Group+(0+Group|Study),
            #                  na.action = na.omit, data=df.ma.lst[[b]], REML=FALSE), silent=TRUE)
            n.perm <- 10000
            tic()
            if(!(inherits(mod0, "try-error"))){
              suppressWarnings({
                mod.p.pv <- try(permlmer.tplh2(mod0, mod, nperm = n.perm, ncore=16, pvalue=TRUE, perc=TRUE), silent=TRUE) ## thien phuc long hao
              })
              while(inherits(mod.p.pv, "try-error")){
                suppressWarnings({
                  mod.p.pv <- try(permlmer.tplh2(mod0, mod, nperm = n.perm, ncore=16, pvalue=TRUE, perc=TRUE), silent=TRUE) ## thien phuc long hao
                })
              }
            }
            toc()
            if(!(inherits(mod.p.pv, "try-error"))){
              mod.p[3] <- mod.p.pv$`97.5-perc`
              mod.p[4] <- mod.p.pv$`2.5-perc`
              mod.p[5] <- mod.p.pv$`Perm-p`[2]
              
              mod.p[1:2] <- mod.beta["Group", c("Estimate", "Std. Error")]
            }
            
            print(mod.p[3:5])
            
            mod.p.df[b,-c(6,7)] <- mod.p
            
          }
          print(b)
          
          #   }
          #   
          # }
          mod.n.df$low <- mod.n.df$est-mod.n.df$`97.5-perc`*mod.n.df$se
          mod.n.df$hi <- mod.n.df$est+mod.n.df$`97.5-perc`*mod.n.df$se
          
          mod.s.df$low <- mod.s.df$est-mod.s.df$`97.5-perc`*mod.s.df$se
          mod.s.df$hi <- mod.s.df$est+mod.s.df$`97.5-perc`*mod.s.df$se
          
          mod.kr.df$low <- mod.kr.df$est-mod.kr.df$`97.5-perc`*mod.kr.df$se
          mod.kr.df$hi <- mod.kr.df$est+mod.kr.df$`97.5-perc`*mod.kr.df$se
          
          mod.kh.df$low <- mod.kh.df$est-mod.kh.df$`97.5-perc`*mod.kh.df$se
          mod.kh.df$hi <- mod.kh.df$est+mod.kh.df$`97.5-perc`*mod.kh.df$se
          
          mod.p.df$low <- mod.p.df$est-mod.p.df$`97.5-perc`*mod.p.df$se
          mod.p.df$hi <- mod.p.df$est-mod.p.df$`2.5-perc`*mod.p.df$se  
          
          if(b%%10==0){
            #save.image(paste("pvalue_epsilonnorm_theta", theta, "_n", n, "_tau", tau[k], "_n.subject", n.subject[j], "_nperm10000_lmer_TPLH2.RData", sep=""))
            #save.image(paste("pvalue_diffuisd_theta", theta, "_n", n, "_tau", tau[k], "_n.subject", n.subject[j], "_nperm10000_lmer_TPLH2.RData", sep=""))
            #save.image(paste("pvalue_epsilont_theta", theta, "_n", n, "_tau", tau[k], "_n.subject", n.subject[j], "_nperm10000_lmer_TPLH2.RData", sep=""))
            #save.image(paste("pvalue_epsilonlognorm_theta", theta, "_n", n, "_tau", tau[k], "_n.subject", n.subject[j], "_nperm10000_lmer_TPLH2.RData", sep=""))
            
          }
        }
        
      }
    }
  }
}