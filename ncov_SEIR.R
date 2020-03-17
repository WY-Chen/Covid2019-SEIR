# install.packages("remotes")
# remotes::install_github("GuangchuangYu/nCov2019")
#get_nCov2019("en")
dat = load_nCov2019("en")
library(ggplot2)
library (deSolve) 
library(plyr)
seir_model = function (current_timepoint, state_values, parameters)
{
  # create state variables (local variables)
  S = state_values [1]        # susceptibles
  E = state_values [2]        # exposed
  I = state_values [3]        # infectious
  R = state_values [4]        # recovered
  N = state_values [5]        # total population, after immigration
  D = state_values [6]        # public perception of risk
  C = state_values [7]        # number of true cumulative cases 
  
  with ( 
    as.list (parameters),     # variable names within parameters can be used 
    {
      # compute derivatives
      b = beta*(1-alpha)*(1-D/N)^kappa
      dS = (-beta * S * F/N - b*S*I/N-mu*S)
      dE = (beta * S * F/N+b*S*I/N) - (sigma+mu )* E
      dI = (sigma * E) - (gamma+mu )* I
      dR = (gamma * I) - mu*R
      dN = -mu*N
      dD = d*gamma*I-lambda*D
      dC = sigma*E
      
      # combine results
      results = c (dS, dE, dI, dR,dN,dD,dC)
      list (results)
    }
  )
}


sim_SEIR<-function(
  region="South Korea", # province in china, or country
  breakpt=NULL,         # change of policy
  alpha=0.4,            # levels of gov reaction, vector input of length(breakpt)+1
  F0=0,                 # number of zootonic cases
  N0=51.47e6,           # population of region 
  I0=1,                 # initial # of cases
  offset=0,             # offset to SEIR model
  window=100,           # "forcast" window: not real forcast, no data is used
  mu=0,                 # rate of emmigration, usually 0. 
  report_rate=1,        # vector: test ability, e.g., (0.2,0.8,1)
  report_breakpt=NULL,  # vector: report rate change date
  explanation=NULL
)
{
  # match region name
  if (region %in% unique(dat$global$country)){
    dat_local = dat$global[which(dat$global[,2]==region),]
  } else {
    dat_local = dat$province[which(dat$province[,2]==region),]
  }
  dat_local<-na.omit(dat_local)
  stopifnot(length(breakpt)+1==length(alpha))
  K = length(alpha)
  
  S0 = 0.9 *N0 # Initial susceptible population
  E0 = I0
  R0 = 0
  D0 = I0
  C0 = I0
  kappa = 1117.3 # Intensity of responds
  mu = 0 # Emigration rate (0,0.0205)
  sigma = 1/3 # 1/Mean latent period
  gamma = 1/5 # Mean infectious period
  d = 0.2 # Proportion of severe cases
  lambda = 1/11.2 # Mean duration of public reaction
  beta0 = 2.8*(sigma+mu)*gamma/sigma # Transmission rate, range (0.59, 1.68)
  
  t=0
  SEIR_output=NULL
  initial_values = c (S = S0, E = E0, I = I0, R = R0,
                      N = N0, D=D0,C=C0)
  for (j in 1:K){
    timepoints = seq (t, c(breakpt,2000)[j]+offset, by=1)
    parameter_list = c(beta = beta0,
                       gamma = gamma,
                       alpha=alpha[j],mu=mu,sigma=sigma,
                       mu-mu,kappa=kappa,d=d,lambda=lambda)
    output = lsoda (initial_values, timepoints, seir_model, parameter_list)
    l = dim(output)[1]
    initial_values = output[l,-1]
    SEIR_output = rbind(SEIR_output,output[-1,])
    t = breakpt[j]+offset
  }
  # offset
  SEIR_output=as.data.frame(SEIR_output)
  SEIR_output= SEIR_output[seq(offset,offset+window),]
  SEIR_output$time=seq(dim(SEIR_output)[1])
  SEIR_output$new = c(0,diff(SEIR_output$C))
  if (length(report_breakpt)){
    report_rate=c(report_rate[1],report_rate)
    report_prd=c(report_breakpt[1],
                 diff(c(report_breakpt,dim(SEIR_output)[1])))
    report_rate = 
      unlist(sapply(
        1:(length(report_rate)-1), 
        function(i){exp(log(seq(report_rate[i],
                                report_rate[i+1],
                                length.out = report_prd[i])))})) 
    # not enough kit in the beginning
  } else {report_rate=rep(1,dim(SEIR_output)[1])}
  SIM_report = data.frame(time = seq(dim(SEIR_output)[1]),
                          new = rep(0,dim(SEIR_output)[1]))
  tot_rpt=0
  for(i in 1:dim(SEIR_output)[1]){
    SIM_report$new[i]= report_rate[i]*(SEIR_output$C[i]-tot_rpt)
    tot_rpt=tot_rpt+report_rate[i]*(SEIR_output$C[i]-tot_rpt)
  }
  
  dat_local$time =seq (0, dim(dat_local)[1]-1, by=1)
  dat_local$new= c(0,diff(dat_local$cum_confirm))
  pl_local=ggplot(mapping = aes(x=time, y=new)) + 
    theme_bw()+
    geom_point(data = dat_local, aes(col = 'True Data')) +
    geom_line(data = SEIR_output, aes(col = 'SEIR'),size=1)+
    geom_line(data = SIM_report, aes(col = 'SEIR_report'))+
    labs(x = "Days since 1st case", y="New cases",
         caption = explanation)+
    ggtitle(paste0("True Data vs. SEIR Simulation in ",region))
  plot(pl_local)
  ggsave(paste0(region,".png"),pl_local)
}

sim_SEIR("Hubei",N0=55e6,alpha=c(0,0.42,0.84),
         breakpt = c(54,64),window=100,
         offset=0,I0=5,F0=10,mu=0.02,
         report_rate = c(0.02,0.3,1),
         report_breakpt = c(54,77),
         explanation = 
           "Gov. Action Strength change on Jan. 23 & Feb 3 (level=0.42 & 0.84), 
         testing capacity change on Jan. 23 & Feb 16 (propotion=1%,30%,100%),
         emigration rate=0.02%")

sim_SEIR("South Korea",N0=51.47e6,alpha=c(0,0.85),
         breakpt = c(41),window=100,
         report_rate = c(0.3,1),
         report_breakpt = c(41),
         offset=0,I0=6,F0=0,
         explanation = "Gov. Action Strength change on Feb. 24 (level=0.85),
         testing capacity change on Feb. 24 by drive-throughs (propotion=30%,100%).")

sim_SEIR("Iran",N0=81.16e6,alpha=c(0),
         breakpt = NULL,window=25,
         report_rate = c(0.2,1),
         report_breakpt = c(9),
         offset=20,I0=20,F0=0,
         explanation = "There might be a lag in reporting, use offset=20days,
         assuming no gov action,
         testing capacity changed on Feb 28 by WHO (proportion=20%,100%)")

sim_SEIR("Italy",N0=60.48e6,alpha=c(0,0.2,0.7),
         breakpt = c(15,33),window=80,
         report_rate = c(0.1,1),
         report_breakpt = c(33),
         offset=20,I0=10,F0=0,
         explanation = 
           "There might be a lag in reporting, use offset=20days,
           Gov. Action Strength change on Feb. 21 & Mar. 10 (levels=0.2,0.7),
         testing capacity change on Mar. 10 (proportion=10%,100%)")

sim_SEIR("Italy",N0=60.48e6,alpha=c(0,0.24,0.6),
         breakpt = c(15,33),window=80,
         report_rate = c(0.1,1),
         report_breakpt = c(33),
         offset=14,I0=25,F0=0,
         explanation = 
           "Offset=2 weeks. Gov. Action Strength change on Feb. 21 & Mar. 10 (level=0.6),
         testing capacity change on Mar. 10")

sim_SEIR("United States",N0=327.2e6,alpha=c(0,0.4),
         breakpt = 50,window=80,
         report_rate = c(0.3,0.5),
         report_breakpt = c(50),
         offset=0,I0=1,F0=0,
         explanation = 
           "Gov. Action Strength change on Mar. 13 (level=0.4),
         testing capacity change on Mar. 13 (proportion=30%,50%)")



sim_SEIR("China",N0=14e9,alpha=c(0,0.42,0.84),
         breakpt = c(54,64),window=100,
         offset=0,I0=5,F0=10,mu=0.02,
         report_rate = c(0.02,0.3,1),
         report_breakpt = c(54,77),
         explanation = 
           "Gov. Action Strength change on Jan. 23 & Feb 3 (level=0.42 & 0.84), 
         testing capacity change on Jan. 23 & Feb 16 (propotion=1%,30%,100%)")


sim_SEIR("Hong Kong",N0=7.4e6,alpha=c(0.5),
         breakpt = NULL,window=100,
         offset=0,I0=1,F0=0,mu=0,
         report_rate = c(0.8),
         report_breakpt = NULL,
         explanation = 
           "Gov. Action Strength  (level=0.5), 
         testing capacity (propotion=80%)")


sim_SEIR("France",N0=66.99e6,alpha=0,
         breakpt = NULL,window=60,
         offset=0,I0=2,F0=0,mu=0,
         report_rate = c(0.8),
         report_breakpt = NULL,
         explanation = 
           "Gov. Action Strength  (level=0), 
         testing capacity (propotion=80%)")

sim_SEIR("Spain",N0=46.66e6,alpha=0,
         breakpt = NULL,window=50,
         offset=14,I0=2,F0=0,mu=0,
         report_rate = c(0.6),
         report_breakpt = NULL,
         explanation = 
           "Gov. Action Strength  (level=0), 
         testing capacity (propotion=80%)")


sim_SEIR("United Kingdom",N0=66.44e6,alpha=0,
         breakpt = NULL,window=50,
         offset=0,I0=3,F0=0,mu=0,
         report_rate = c(0.8,1),
         report_breakpt = 50,
         explanation = 
           "Gov. Action Strength  (level=0), 
         testing capacity (propotion=80%)")
