

install.packages(TMB)
install.packages(rstan)





df <- read.csv("data/examplesr.csv")


compile("src/Ricker_simple.cpp")
dyn.load(dynlib("src/Ricker_simple"))





## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simPar$species),sep = "_"))



parameters_simple<- list(
    alpha=srm$coefficients[1],
    logbeta = log(ifelse(-srm$coefficients[2]<0,1e-08,-srm$coefficients[2])),
    logsigobs=log(.4)
    )

  obj_simple <- MakeADFun(SRdata,parameters_simple,DLL="Ricker_simple")
  newtonOption(obj_simple, smartsearch=FALSE)

  opt_simple <- nlminb(obj_simple$par,obj_simple$fn,obj_simple$gr)
  rep_simple <- obj_simple$report()



 datm = list(N=nrow(df),
                R_S =df$logRS,
                S=df$S)

fit1 <- stan(
  file = "src/ricker_linear.stan",  # Stan program
  data = df,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 2              # number of cores (could use one per chain)

  )


