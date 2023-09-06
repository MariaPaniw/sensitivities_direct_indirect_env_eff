library(boot)
library(popbio)

#### SCENARIO 1 ######

# Intraspecific interactions become more negative with better environmental conditions

# Survival juveniles
S.j <- function(species,env,N1,N2){
  
  if(species==1){
    b1 <- rnorm(1, 0-0.02*env,0.001)
    surv <- inv.logit(3.0 + 0.8*env - b1*N1 - 0.03*N2 - 0.03*env*N1)
  }
    
  
  if(species==2){
    b1 <- rnorm(1, 0-0.01*env,0.001)
    surv <- inv.logit(3.5 + 0.9*env - 0.08*N2 - b1*N1 - 0.02*env*N2)
  }
    
  
  return(surv)
}


# Survival reproductive inidividuals
S.r <- function(species,env,N1,N2){
  
  if(species==1)
    
    surv <- inv.logit(4.5 + 0.2*env - 0.05*N1 - 0.01*N2 - 0.01*env*N1)
  
  if(species==2)
    
    surv <- inv.logit(4.0 + 0.4*env - 0.06*N2 - 0.04*N1 - 0.02*env*N2)
  
  return(surv)
 
}

# Survival non-reproductive inidividuals
S.n <- function(species,env,N1,N2){
  
  if(species==1)
    
    surv <- inv.logit(4.1 + 0.3*env - 0.07*N1 - 0.02*N2)
  
  if(species==2)
    
    surv <- inv.logit(3.8 + 0.3*env - 0.05*N2 - 0.02*N1)
  
  return(surv)

}

# Transition of juveniles to reproductive 
T.j <- function(species,env,N1,N2){
  
  if(species==1){
    
    b1 <- rnorm(1, 0-0.05*env,0.06)
    
    trans <- inv.logit(0.5 + 0.01*env - b1*N1 - 0.01*N2)
  }
    
    
  
  if(species==2){
    
    b1 <- rnorm(1, 0-0.01*env,0.005)
    
    trans <- inv.logit(1 + 0.03*env - 0.1*N2 - b1*N1)
  }
    
    
  
  return(trans)
  
 
}

# Transition of non-reproductive to reproductive 
T.n <- function(species,env,N1,N2){
  
  if(species==1)
    
    trans <- inv.logit(1.2 + 0.01*env - 0.1*N1 - 0.02*N2)
  
  if(species==2)
    
    trans <- inv.logit(1.9 + 0.02*env - 0.1*N2 - 0.02*N1)
  
  return(trans)
 
}

# Number of offsrping (recruitment)
Rec <- function(species,env,N1,N2){
  
  if(species==1)
    
    off <- exp(1.2 + 0.05*env - 0.01*N1 - 0.005*env*N1)
  
  if(species==2)
    
    off <- exp(2 + 0.05*env - 0.01*N1 - 0.02*env*N2)
  
  return(off)

}

### START SIMULATIONS

# Initial condition
years=100 # years (or time steps) of simulations

n = 3 # number of stages
# simulate environment from normal distribution

# set.seed(128)
env=rnorm(years, mean=0, sd=1)

## Initial S1

densS1=c(30,30,25) 

dens.effS1=sum(densS1)

## Initial S2

densS2=c(35,35,20) 

dens.effS2=sum(densS2)


sim.data=array(0,c(2,n,years)) # stage- specific abundances for species and sites

for(i in 1:years){
  
  sim.data[1,,i]= densS1
  sim.data[2,,i]= densS2
  
  mpm.S1= matrix(c(0,0,S.r(1,env[i],dens.effS1,dens.effS2)*Rec(1,env[i],dens.effS1,dens.effS2),
                 S.j(1,env[i],dens.effS1,dens.effS2)*(1-T.j(1,env[i],dens.effS1,dens.effS2)),S.n(1,env[i],dens.effS1,dens.effS2)*(1-T.n(1,env[i],dens.effS1,dens.effS2)),0,
                 S.j(1,env[i],dens.effS1,dens.effS2)*(T.j(1,env[i],dens.effS1,dens.effS2)),S.n(1,env[i],dens.effS1,dens.effS2)*(T.n(1,env[i],dens.effS1,dens.effS2)),S.r(1,env[i],dens.effS1,dens.effS2)),3,3,byrow = T)
  
  
  mpm.S2= matrix(c(0,0,S.r(2,env[i],dens.effS1,dens.effS2)*Rec(2,env[i],dens.effS1,dens.effS2),
                 S.j(2,env[i],dens.effS1,dens.effS2)*(1-T.j(2,env[i],dens.effS1,dens.effS2)),S.n(2,env[i],dens.effS1,dens.effS2)*(1-T.n(2,env[i],dens.effS1,dens.effS2)),0,
                 S.j(2,env[i],dens.effS1,dens.effS2)*(T.j(2,env[i],dens.effS1,dens.effS2)),S.n(2,env[i],dens.effS1,dens.effS2)*(T.n(2,env[i],dens.effS1,dens.effS2)),S.r(2,env[i],dens.effS1,dens.effS2)),3,3,byrow = T) 
  
 
  densS1=mpm.S1%*%densS1
  dens.effS1=sum(densS1) 
  
  densS2=mpm.S2%*%densS2
  dens.effS2=sum(densS2) 

}


library(plyr)
library(ggplot2)

df=adply(sim.data,c(1,2,3))

colnames(df)=c("species","stage","year","density")

levels(df$stage)=c("J","N","R")
levels(df$species)=c("S1","S2")
df$year=as.numeric(df$year)


ggplot(df[df$year>3,],aes(year,density,col=species))+
  geom_line()+
  facet_grid(stage~.,scales = "free")+
  scale_color_manual(name="",values=c("darkgreen","orange"))+
  xlab("Simulation year")+ylab("Abundance")+theme_bw(base_size=20)+
  theme(panel.grid = element_blank())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
