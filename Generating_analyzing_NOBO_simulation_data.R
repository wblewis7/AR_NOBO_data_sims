# Generating simulation data --------------------------------------------------

  library(truncnorm)
  
  ## Simulation Parameters -----------------------------------------------------
  nsims <- 100 # Number of simulations
  AvgDens <- seq(0.1,2.5,by=0.4)
  len.AvgDens <- length(AvgDens)
  Ngrids <- c(2,4,6,8,10)
  ProportionArea <- round(Ngrids*600*600*pi/10000/2023.428,2) # 11%, 22%, 34%, 45% and 56% of area
  len.surveys <- length(Ngrids)

  Nvisits <- c(2,3,4,5) # Number of repeat visits for point counts
  NdaysARU <- c(7,14,21,28) # Number of days of ARU deployment
  Propsurvey <- c(0,0.5,1)
  # Number of visits, number of ARU days, and proportion of ARUs/counts is not factorial. e.g., can't have
  #   different visits without point counts or different ARU deployments without ARUs
  scenarios <- as.data.frame(matrix(NA,nrow=0,ncol=4))
  colnames(scenarios) <- c("PC","ARU","Visits","ARUday")
  for(i in 1:length(Propsurvey)){
    for(j in 1:length(Propsurvey)){
      if((Propsurvey[i] + Propsurvey[j])>=1){
        if(Propsurvey[i]>0 & Propsurvey[j]>0){
        scenarios.temp <- data.frame(PC=Propsurvey[i],
                                       ARU=Propsurvey[j],
                                       Visits=rep(Nvisits,times=length(NdaysARU)),
                                       ARUday=rep(NdaysARU,each=length(Nvisits)))
        } else if(Propsurvey[i]>0 & Propsurvey[j]==0){
         scenarios.temp <- data.frame(PC=Propsurvey[i],
                                       ARU=Propsurvey[j],
                                       Visits=Nvisits,
                                       ARUday=NA)
        } else{
          scenarios.temp <- data.frame(PC=Propsurvey[i],
                                      ARU=Propsurvey[j],
                                      Visits=NA,
                                      ARUday=NdaysARU)
        }
        scenarios <- rbind(scenarios, scenarios.temp)
      }
    }
  }
  N.scenarios <- nrow(scenarios)

  
  ## State Parameters ----------------------------------------------------------
  
  # Need to convert from bird density to covey density. Average of 12 birds/covey
  avg.covey.size <- 12
  
  # Point counts truncated at 600m radius, preliminary analysis of playback of bobwhite
  #   calls at known distance from ARU suggests that detection radius is 520m.
  PointAreaARU <- pi * 520 * 520
  PointAreaPC <- pi * 600 * 600
  PointAreaDiff <- PointAreaPC - PointAreaARU
  
  # Point count detection process. Half-normal detection function
  PC.sigma <- 193
  db <- seq(0,600,by=100) 
  nB <- length(db) - 1                   
  xg <- NA 
  for(j in 1:(length(db)-1)){
    xg[j] <- (db[j] + db[j+1])/2
  }
  delta.sig <- NA    
  for(j in 1:length(db)-1){
    delta.sig[j] <- db[j+1] - db[j]
  }
  pix <- (2*xg*delta.sig)/(max(db)*max(db))         # bin area proportions
  pd_bins <- rep(NA, times=nB)
  for(b in 1:nB){
    pd_bins[b] <- (PC.sigma^2*(1-exp(-db[b+1]^2/(2*PC.sigma^2)))-PC.sigma^2*(1-exp(-db[b]^2/(2*PC.sigma^2))))*2*3.1416/(PointAreaPC*pix[b])
  }
  
  # ARU hurdle model covariates. From analysis of 2022 Arkansas data.
  alpha.0 <- -3.23 # Probability of detecting a false positive when N=0
  alpha.1 <- 1.3 # How probability of detecting at least one call changes with abundance
  
  # Putting covariates on delta for background noise. Simulating values each day. 
  cov.mean <- 0
  cov.sd <- 1
  
  omega <- 10.5 # Mean number of false positives per recording
  
  # Modifying methods of Doser et al. (2021) for generating the number of ARU detections
  #   per day. Makes more biologically-realistic data.
  delta.b0 <- -1.03 # Intercept (logit scale) for probability of at least one detection on a daily ARU recording at an occupied site
  delta.b1 <- 0.87 # Effect (logit scale) of abundance on probability of at least one detection on a daily ARU recording at an occupied site
  delta.mean <- 21.7 # Mean delta (average calling rate/detection parameter) when at least one detection is recorded
  delta.sd <- 15.1 # Standard deviation of delta when at least one detection is recorded
  delta.min <- 1.4 # 2.5% quantile of delta when at least one detection is recorded
  delta.max <- 55.7 # 97.5% quantile of delta when at least one detection is recorded
  
  # Effect of standardized background noise on delta
  gamma.2 <- -0.23
  
  # Only a proportion of days for ARU deployment have suitable weather (not raining and wind speed <16 kmph)
  p.suitweather <- 0.81
  
  
  ## Simulations ---------------------------------------------------------------

  for(z in 1:len.AvgDens){
    
    # Calculating covery density from bird density
    dens.c <- AvgDens[z]/avg.covey.size
    
    for(i in 1:len.surveys){
      for(j in 1:N.scenarios){
        
        PC.N.det.data <- PC.det.distbin.data <- PC.N.distbin.act <- ARU.noise.data <- ARU.y.data <- ARU.v.data <- N.calls.ARU.data <- omega.testing.data <- delta.act <- GridwARU.act <- A.times.data <- n.A.times.data <- GridwPC.act <- vector(mode="list", length=nsims)
        covey.N.ARU <- covey.N.PC <- matrix(NA, nrow=nsims, ncol=Ngrids[i])
        Ndaysgood.dat <- rep(NA, times=nsims)
        
        for(x in 1:nsims){
          
          # Simulating abundance for ARUs and count areas. Ensuring that
          #    abundance on counts is equal or greater than ARUs
          covey.N.ARU[x,] <- rpois(Ngrids[i], PointAreaARU/10000 * dens.c)
          covey.N.PC[x,] <- rpois(Ngrids[i], PointAreaDiff/10000 * dens.c) + covey.N.ARU[x,]
          
          
          ### Point counts (if doing) ------------------------------------------
          if(scenarios$PC[j]>0){
            
            nGridwPC <-scenarios$PC[j] * Ngrids[i]
            # Selecting which Grid locations get point counts
            GridwPC <- sample(1:Ngrids[i], nGridwPC, replace=F)
            GridwPC <- GridwPC[order(GridwPC)]
            
            # Allowing distance band of coveys to vary by visit
            PC.N.distbin <- PC.det.distbin <- array(NA, dim=c(nGridwPC,nB,scenarios$Visits[j]))
            PC.N.det <- matrix(NA, nrow=nGridwPC, ncol=scenarios$Visits[j])
            for(t in 1:nGridwPC){
              for(m in 1:scenarios$Visits[j]){
                PC.N.distbin[t,1:nB,m] <- rmultinom(1,covey.N.PC[x,GridwPC[t]],pix)
                for(r in 1:nB){
                  PC.det.distbin[t,r,m] <- rbinom(1,PC.N.distbin[t,r,m],pd_bins[r])
                }
                PC.N.det[t,m] <- sum(PC.det.distbin[t,1:nB,m])
              }
            }
            
            GridwPC.act[[x]] <- GridwPC
            PC.N.det.data[[x]] <- PC.N.det
            PC.det.distbin.data[[x]] <- PC.det.distbin
            PC.N.distbin.act[[x]] <- PC.N.distbin
          }
          
          
          ### ARUs (if doing) --------------------------------------------------
          if(scenarios$ARU[j]>0){
            
            nGridwARU <-scenarios$ARU[j] * Ngrids[i]
            # If ratio is 50/50 with PC and ARUs, ensuring that ARU locations are those not surveyed by PC
            if(scenarios$PC[j]==0.5 & scenarios$ARU[j]==0.5){
              GridwARU <- (1:Ngrids[i])[!(1:Ngrids[i] %in% GridwPC)]
            } else{
              GridwARU <- sample(1:Ngrids[i], nGridwARU, replace=F)
            }
            GridwARU <- GridwARU[order(GridwARU)]
            
            # Simulating days of suitable weather during the deployment period
            gooddays <- matrix(rbinom(scenarios$ARUday[j],1,p.suitweather),nrow=1)
            gooddays[gooddays==0] <- NA
            Ndaysgood.dat[x] <- sum(gooddays,na.rm=T)
            
            # Simulating standardized background noise for each day of recording
            ARU.noise <- matrix(rnorm(nGridwARU * scenarios$ARUday[j], cov.mean, cov.sd), nrow=nGridwARU)
            for(b in 1:nrow(ARU.noise)){
              ARU.noise[b,] <- ARU.noise[b,] * gooddays
            }
            
            ARU.y <- delta.temp <- delta.z <- delta2 <- omega.dat <- delta.val <- ARU.v <- A.times <- matrix(NA, nrow=nGridwARU, ncol=scenarios$ARUday[j])
            n.A.times <- rep(NA, times=nGridwARU)
            for(t in 1:nGridwARU){
              
              # Probability of detecting at least one call on a recording
              pa.temp <- plogis(alpha.0 + alpha.1*covey.N.ARU[x,GridwARU[t]])
              
              for(q in 1:scenarios$ARUday[j]){
                
                if(!is.na(gooddays[q])){
                  ARU.y[t,q] <- rbinom(1, 1, pa.temp)
                  
                  # Delta will be 0 if no calls detected or abundance is 0
                  # Adjusting mean based on effect of wind speed
                  delta.temp[t,q] <- ifelse(covey.N.ARU[x,GridwARU[t]]==0,0,ARU.y[t,q]*rtruncnorm(1, a=delta.min, b=delta.max, mean=delta.mean+gamma.2*ARU.noise[t,q], sd=delta.sd))
                  
                  # For those sites with N > 1, sometimes will only have false positives (number of detections based on omega)
                  #     and sometimes have actual calls (number of detections based on delta and omega)
                  delta.z[t,q] <- ifelse(ARU.y[t,q]==1 & covey.N.ARU[x,GridwARU[t]],rbinom(1,1,plogis(delta.b0 + delta.b1*covey.N.ARU[x,GridwARU[t]])),0)
                  delta2[t,q] <- delta.temp[t,q] * delta.z[t,q] * ARU.y[t,q]
                  
                  # Number of false positives per recording
                  omega.dat[t,q] <- rpois(1, omega) * ARU.y[t,q]
                  
                  # Actual delta values
                  delta.val[t,q] <- ifelse(ARU.y[t,q]==1,delta2[t,q],NA)
                  
                  # ARU call detection data
                  ARU.v[t,q] <- rpois(1, delta2[t,q]*covey.N.ARU[x,GridwARU[t]]) + omega.dat[t,q]
                }
              }
              
              n.A.times[t] <- length(which(ARU.v[t,]>0))
              if(n.A.times[t]>0){
                A.times[t,1:n.A.times[t]] <- which(ARU.v[t,]>0)
              }
            }
            
            N.calls.ARU <- rowSums(ARU.v, na.rm=T)
            
            
            ### False positive validation -------------------------------------
            
            # Randomly picking 5 (or max, if <5) recordings with detections to assess false positives/recording (omega)
            omega.test <- which(omega.dat>0)
            if(length(omega.test)>0){
              if(length(omega.test)==1){
                N.omega.test.samp <- omega.test
              } else{
                N.omega.test.samp <- sample(omega.test, min(5,length(omega.test)), replace=F)
              }
              omega.col <- floor(N.omega.test.samp/nrow(omega.dat)-0.0000001)+1
              omega.row <- N.omega.test.samp - nrow(omega.dat)*(omega.col-1)
              omega.testing <- matrix(NA, nrow=length(N.omega.test.samp), ncol=4)
              for(p in 1:nrow(omega.testing)){
                omega.testing[p,1:2] <- c(omega.row[p],omega.col[p])
                omega.testing[p,3] <- ARU.v[omega.row[p],omega.col[p]]
                omega.testing[p,4] <- omega.dat[omega.row[p],omega.col[p]]
              }
              omega.testing <- as.data.frame(omega.testing)
            } else {
              omega.testing <- as.data.frame(matrix(NA, nrow=0, ncol=4))
            }
            colnames(omega.testing) <- c("Grid","Day","Ncalls","FP")
            
            
            GridwARU.act[[x]] <- GridwARU
            ARU.noise.data[[x]] <- ARU.noise
            ARU.y.data[[x]] <- ARU.y
            ARU.v.data[[x]] <- ARU.v
            A.times.data[[x]] <- A.times
            n.A.times.data[[x]] <- n.A.times
            N.calls.ARU.data[[x]] <- N.calls.ARU
            omega.testing.data[[x]] <- omega.testing
            delta.act[[x]] <- delta.val
          }
        }
        
        ### Saving out ---------------------------------------------------------
        NOBO.sim.PC.ARU.dat <- list(Avg.Density=AvgDens[z],
                                    N.Grid=Ngrids[i],
                                    ProportionArea.PC=ProportionArea[i],
                                    N.PC=Ngrids[i]*scenarios$PC[j],
                                    Nvisits.PC=scenarios$Visits[j],
                                    N.ARUs=Ngrids[i]*scenarios$ARU[j],
                                    ARU.days=scenarios$ARUday[j],
                                    Nsims=nsims,
                                    PointAreaARU.ha=PointAreaARU/10000,
                                    PointAreaDiff.ha=PointAreaDiff/10000,
                                    PointAreaPC.ha=PointAreaPC/10000,
                                    db=db,
                                    nB=nB,
                                    pix=pix,
                                    PC.sigma.act=PC.sigma,
                                    avg.covey.size=avg.covey.size,
                                    N.PC.act=covey.N.PC,
                                    Grid.w.PC.act=GridwPC.act,
                                    PC.N.distbin.act=PC.N.distbin.act,
                                    PC.det.distbin.data=PC.det.distbin.data,
                                    PC.N.det.data=PC.N.det.data,
                                    N.ARU.act=covey.N.ARU,
                                    Grid.w.ARU.act=GridwARU.act,                                    
                                    omega.act=omega,
                                    ARU.alpha.0.act=alpha.0,
                                    ARU.alpha.1.act=alpha.1,
                                    ARU.delta.mean.act=delta.mean,
                                    ARU.delta.sd.act=delta.sd,
                                    ARU.delta.min.act=delta.min,
                                    ARU.delta.max.act=delta.max,
                                    ARU.gamma.2.act=gamma.2,
                                    ARU.delta.b0.act=delta.b0,
                                    ARU.delta.b1.act=delta.b1,
                                    delta.act=delta.act,
                                    ARU.backgroundnoise.data=ARU.noise.data,
                                    ARU.y.data=ARU.y.data,
                                    ARU.v.data=ARU.v.data,
                                    A.times.data=A.times.data,
                                    n.A.times.data=n.A.times.data,
                                    omega.testing.data=omega.testing.data,
                                    N.calls.ARU.data=N.calls.ARU.data)
        
        filename <- paste0("NOBO_simdata_",AvgDens[z],"Dens_",Ngrids[i],"Grids_",scenarios$PC[j]*Ngrids[i],"PC_",scenarios$Visits[j],"visits_",scenarios$ARU[j]*Ngrids[i],"ARUs_",scenarios$ARUday[j],"ARUdays.gzip")
        save(NOBO.sim.PC.ARU.dat, file=filename)
        
      }
    }
  }
  
  
  
  
  
  
# Analyzing simulated data

  library(dclone)
  
  for(z in 1:len.AvgDens){
    for(i in 1:len.surveys){
      for(j in 1:N.scenarios){
        
        load(paste0("NOBO_simdata_",AvgDens[z],"Dens_",Ngrids[i],"Grids_",scenarios$PC[j]*Ngrids[i],"PC_",scenarios$Visits[j],"visits_",scenarios$ARU[j]*Ngrids[i],"ARUs_",scenarios$ARUday[j],"ARUdays.gzip"))
        
        N_PC.samps <- N_ARU.samps <- b.0.samps <- gamma.1.samps <- sigma.samps <- omega.samps <- GR <- vector(mode="list",length=NOBO.sim.PC.ARU.dat$Nsims)
        
        # Looping through simulations
        for(x in 1:NOBO.sim.PC.ARU.dat$Nsims){
          
          print(paste0("Doing simulation ",x,": ",NOBO.sim.PC.ARU.dat$Avg.Density,"Dens_",NOBO.sim.PC.ARU.dat$N.PC,"PC_",NOBO.sim.PC.ARU.dat$Nvisits,"visits_",NOBO.sim.PC.ARU.dat$N.ARUs,"ARUs_",NOBO.sim.PC.ARU.dat$ARU.days,"ARUdays",": ",Sys.time()))
          
          ### If have point counts and ARUs -----------------------------------
          if(NOBO.sim.PC.ARU.dat$N.ARUs>0 & NOBO.sim.PC.ARU.dat$N.PC>0){
            sink("NOBO_AVC_sim.jags")
            cat("
            model{
    
              b.0 ~ dnorm(0, 0.001)
              beta ~ dunif(0,10)
              log(sigma) <- beta
              a.phi ~ dunif(0,100)
              mu.alpha ~ dunif(0,1)
              alpha.0 <- logit(mu.alpha)
              alpha.1 ~ dunif(0,1000)
              gamma.0 ~ dnorm(0,0.001)
              gamma.1 ~ dunif(-1000,0)
              omega.log ~ dunif(0,10)
              log(omega) <- omega.log
              for(x in 1:nvalidate){
                FP[x] ~ dpois(omega)
              }
              for(b in 1:nB){
                pd_bins[b] <- (sigma^2*(1-exp(-db[b+1]^2/(2*sigma^2)))-sigma^2*(1-exp(-db[b]^2/(2*sigma^2))))*2*3.1416/(point.area.PC*10000*pix[b])
                pd_bins_adj[b] <- pd_bins[b]*pix[b]
              }
              pd <- sum(pd_bins_adj[1:nB])
            
              ### State Model
              for(i in 1:n.Grids){
                N_diff[i] ~ dpois(lambda_diff[i])
                lambda_diff[i] <- exp(b.0 + log(point.area.DIFF))
                N_ARU[i] ~ dpois(lambda_ARU[i])
                lambda_ARU[i] <- exp(b.0 + log(point.area.ARU))
                N_PC[i] <- N_ARU[i] + N_diff[i]
              }
            
              ### Observation Process

              #### Point Counts
              for(n in 1:n.PC){
                for(j in 1:nvisits){
                  y.PC[n,j] ~ dbin(pd, N_PC[GridwPC[n]])
                  ydb.PC[n,1:nB,j] ~ dmulti(pd_bins_adj[1:nB],y.PC[n,j])
                }
              }

              ##### ARUs
              for(g in 1:n.ARUs){
                logit(p.a[g]) <- alpha.0 + alpha.1*N_ARU[GridwARU[g]]
                for(k in 1:ARU.days){
                  phi[g, k] ~ dgamma(a.phi, a.phi)
                  y.ARU[g, k] ~ dbin(p.a[g], 1)
                }
                for(p in 1:n.A.times[g]){
                  log(delta[g, p]) <- gamma.0 + gamma.1*ARU.noise[g, A.times[g, p]]
                  v.ARU[g, A.times[g, p]] ~ dpois((delta[g, p] * N_ARU[GridwARU[g]] + omega) * phi[g, A.times[g,p]] * y.ARU[g, A.times[g, p]])T(1, )
                }
              }
            }
          ",fill=TRUE)
            sink()
            
            data <- list(y.PC=NOBO.sim.PC.ARU.dat$PC.N.det.data[[x]],
                         ydb.PC=NOBO.sim.PC.ARU.dat$PC.det.distbin.data[[x]],
                         nB=NOBO.sim.PC.ARU.dat$nB,
                         db=NOBO.sim.PC.ARU.dat$db,
                         pix=NOBO.sim.PC.ARU.dat$pix,
                         point.area.PC=NOBO.sim.PC.ARU.dat$PointAreaPC.ha,
                         point.area.ARU=NOBO.sim.PC.ARU.dat$PointAreaARU.ha,
                         point.area.DIFF=NOBO.sim.PC.ARU.dat$PointAreaDiff.ha,
                         n.Grids=NOBO.sim.PC.ARU.dat$N.Grid,
                         n.PC=NOBO.sim.PC.ARU.dat$N.PC,
                         nvisits=NOBO.sim.PC.ARU.dat$Nvisits.PC,
                         nvalidate=nrow(NOBO.sim.PC.ARU.dat$omega.testing.data[[x]]),
                         FP=NOBO.sim.PC.ARU.dat$omega.testing.data[[x]]$FP,
                         n.ARUs=NOBO.sim.PC.ARU.dat$N.ARUs,
                         ARU.days=NOBO.sim.PC.ARU.dat$ARU.days,
                         n.A.times=NOBO.sim.PC.ARU.dat$n.A.times.data[[x]],
                         A.times=NOBO.sim.PC.ARU.dat$A.times.data[[x]],
                         v.ARU=NOBO.sim.PC.ARU.dat$ARU.v[[x]],
                         y.ARU=NOBO.sim.PC.ARU.dat$ARU.y[[x]],
                         GridwPC=NOBO.sim.PC.ARU.dat$Grid.w.PC.act[[x]],
                         GridwARU=NOBO.sim.PC.ARU.dat$Grid.w.ARU.act[[x]],
                         ARU.noise=NOBO.sim.PC.ARU.dat$ARU.backgroundnoise.data[[x]])
            
            inits.fun <- function() list(b.0=runif(1,-2,2),
                                         beta=runif(1,4.9,5.5),
                                         N_ARU=rep((max(NOBO.sim.PC.ARU.dat$PC.N.det.data[[x]],na.rm=T)+1),times=NOBO.sim.PC.ARU.dat$N.Grid),
                                         N_diff=rep((max(NOBO.sim.PC.ARU.dat$PC.N.det.data[[x]],na.rm=T)+1),times=NOBO.sim.PC.ARU.dat$N.Grid),
                                         mu.alpha=runif(1,0.1,0.9),
                                         alpha.1=runif(1,0,1),
                                         gamma.0=runif(1,-2,2),
                                         gamma.1=runif(1,-10,0),
                                         omega.log=runif(1,0,3),
                                         phi=matrix(runif(NOBO.sim.PC.ARU.dat$N.ARUs*NOBO.sim.PC.ARU.dat$ARU.days,0.01,2),nrow=NOBO.sim.PC.ARU.dat$N.ARUs))
            params <- c("N_PC","N_ARU","b.0","gamma.1","sigma","omega")
            
            cl3 <- makePSOCKcluster(3)
            NOBO.model.cs <- jags.parfit(cl3, data=data, params=params,
                                         model="NOBO_AVC_sim.jags", inits=inits.fun(),
                                         n.chains=3, n.iter=120000, n.burnin=30000,n.adapt=30000,
                                         thin=25)
            
            stopCluster(cl3)
            
            NOBO.mcmc.samps <- rbind(NOBO.model.cs[[1]],NOBO.model.cs[[2]],NOBO.model.cs[[3]])
            N_PC.samps[[x]] <- NOBO.mcmc.samps[,grep("N_PC",colnames(NOBO.mcmc.samps))]
            N_ARU.samps[[x]] <- NOBO.mcmc.samps[,grep("N_ARU",colnames(NOBO.mcmc.samps))]
            b.0.samps[[x]] <- NOBO.mcmc.samps[,"b.0"]
            gamma.1.samps[[x]] <- NOBO.mcmc.samps[,"gamma.1"]
            sigma.samps[[x]] <- NOBO.mcmc.samps[,"sigma"]
            omega.samps[[x]] <- NOBO.mcmc.samps[,"omega"]
            GR[[x]] <- gelman.diag(NOBO.model.cs, multivariate=F)$psrf[,1]
            
          } else if(NOBO.sim.PC.ARU.dat$N.ARUs>0 & NOBO.sim.PC.ARU.dat$N.PC==0){
            
            ### ARUs but not counts --------------------------------------------
            
            sink("NOBO_AV_sim.jags")
            cat("
            model{
    
              b.0 ~ dnorm(0, 0.001)
              a.phi ~ dunif(0,100)
              mu.alpha ~ dunif(0,1)
              alpha.0 <- logit(mu.alpha)
              alpha.1 ~ dunif(0,1000)
              gamma.0 ~ dnorm(0,0.001)
              gamma.1 ~ dunif(-1000,0)
              omega.log ~ dunif(0,10)
              log(omega) <- omega.log
              for(x in 1:nvalidate){
                FP[x] ~ dpois(omega)
              }
            
              ### State Model
              for(i in 1:n.Grids){
                N_diff[i] ~ dpois(lambda_diff[i])
                lambda_diff[i] <- exp(b.0 + log(point.area.DIFF))
                N_ARU[i] ~ dpois(lambda_ARU[i])
                lambda_ARU[i] <- exp(b.0 + log(point.area.ARU))
                N_PC[i] <- N_ARU[i] + N_diff[i]
              }
            
              ### Observation Process

              ##### ARUs
              for(g in 1:n.ARUs){
                logit(p.a[g]) <- alpha.0 + alpha.1*N_ARU[GridwARU[g]]
                for(k in 1:ARU.days){
                  phi[g, k] ~ dgamma(a.phi, a.phi)
                  y.ARU[g, k] ~ dbin(p.a[g], 1)
                }
                for(p in 1:n.A.times[g]){
                  log(delta[g, p]) <- gamma.0 + gamma.1*ARU.noise[g, A.times[g, p]]
                  v.ARU[g, A.times[g, p]] ~ dpois((delta[g, p] * N_ARU[GridwARU[g]] + omega) * phi[g, A.times[g,p]] * y.ARU[g, A.times[g, p]])T(1, )
                }
              }
            }
          ",fill=TRUE)
            sink()
            
            data <- list(point.area.ARU=NOBO.sim.PC.ARU.dat$PointAreaARU.ha,
                         point.area.DIFF=NOBO.sim.PC.ARU.dat$PointAreaDiff.ha,
                         n.Grids=NOBO.sim.PC.ARU.dat$N.Grid,
                         nvalidate=nrow(NOBO.sim.PC.ARU.dat$omega.testing.data[[x]]),
                         FP=NOBO.sim.PC.ARU.dat$omega.testing.data[[x]]$FP,
                         n.ARUs=NOBO.sim.PC.ARU.dat$N.ARUs,
                         ARU.days=NOBO.sim.PC.ARU.dat$ARU.days,
                         n.A.times=NOBO.sim.PC.ARU.dat$n.A.times.data[[x]],
                         A.times=NOBO.sim.PC.ARU.dat$A.times.data[[x]],
                         v.ARU=NOBO.sim.PC.ARU.dat$ARU.v[[x]],
                         y.ARU=NOBO.sim.PC.ARU.dat$ARU.y[[x]],
                         GridwARU=NOBO.sim.PC.ARU.dat$Grid.w.ARU.act[[x]],
                         ARU.noise=NOBO.sim.PC.ARU.dat$ARU.backgroundnoise.data[[x]])
            
            inits.fun <- function() list(b.0=runif(1,-2,2),
                                         N_ARU=sample(1:20,NOBO.sim.PC.ARU.dat$N.Grid,replace=T),
                                         N_diff=sample(1:20,NOBO.sim.PC.ARU.dat$N.Grid,replace=T),
                                         mu.alpha=runif(1,0.1,0.9),
                                         alpha.1=runif(1,0,1),
                                         gamma.0=runif(1,-1,1),
                                         gamma.1=runif(1,-1,0),
                                         omega.log=runif(1,0,2.5),
                                         phi=matrix(runif(NOBO.sim.PC.ARU.dat$N.ARUs*NOBO.sim.PC.ARU.dat$ARU.days,0.01,2),nrow=NOBO.sim.PC.ARU.dat$N.ARUs))
            params <- c("N_PC","N_ARU","b.0","gamma.1","omega")
            
            cl3 <- makePSOCKcluster(3)
            NOBO.model.cs <- jags.parfit(cl3, data=data, params=params,
                                         model="NOBO_AV_sim.jags", inits=inits.fun(),
                                         n.chains=3, n.iter=120000, n.burnin=30000,n.adapt=30000,
                                         thin=25)
            
            stopCluster(cl3)
            
            NOBO.mcmc.samps <- rbind(NOBO.model.cs[[1]],NOBO.model.cs[[2]],NOBO.model.cs[[3]])
            N_PC.samps[[x]] <- NOBO.mcmc.samps[,grep("N_PC",colnames(NOBO.mcmc.samps))]
            N_ARU.samps[[x]] <- NOBO.mcmc.samps[,grep("N_ARU",colnames(NOBO.mcmc.samps))]
            b.0.samps[[x]] <- NOBO.mcmc.samps[,"b.0"]
            gamma.1.samps[[x]] <- NOBO.mcmc.samps[,"gamma.1"]
            omega.samps[[x]] <- NOBO.mcmc.samps[,"omega"]
            GR[[x]] <- gelman.diag(NOBO.model.cs, multivariate=F)$psrf[,1]
            
          } else {
            
            
            ### Point counts but no ARUs ---------------------------------------
            
            sink("NOBO_C_sim.jags")
            cat("
            model{
    
              b.0 ~ dnorm(0, 0.001)
              beta ~ dunif(0,10)
              log(sigma) <- beta
              for(b in 1:nB){
                pd_bins[b] <- (sigma^2*(1-exp(-db[b+1]^2/(2*sigma^2)))-sigma^2*(1-exp(-db[b]^2/(2*sigma^2))))*2*3.1416/(point.area.PC*10000*pix[b])
                pd_bins_adj[b] <- pd_bins[b]*pix[b]
              }
              pd <- sum(pd_bins_adj[1:nB])
            
              ### State Model
              for(i in 1:n.Grids){
                N_diff[i] ~ dpois(lambda_diff[i])
                lambda_diff[i] <- exp(b.0 + log(point.area.DIFF))
                N_ARU[i] ~ dpois(lambda_ARU[i])
                lambda_ARU[i] <- exp(b.0 + log(point.area.ARU))
                N_PC[i] <- N_ARU[i] + N_diff[i]
              }
            
              ### Observation Process

              #### Point Counts
              for(n in 1:n.PC){
                for(j in 1:nvisits){
                  y.PC[n,j] ~ dbin(pd, N_PC[GridwPC[n]])
                  ydb.PC[n,1:nB,j] ~ dmulti(pd_bins_adj[1:nB],y.PC[n,j])
                }
              }
            }
          ",fill=TRUE)
            sink()
            
            data <- list(y.PC=NOBO.sim.PC.ARU.dat$PC.N.det.data[[x]],
                         ydb.PC=NOBO.sim.PC.ARU.dat$PC.det.distbin.data[[x]],
                         nB=NOBO.sim.PC.ARU.dat$nB,
                         db=NOBO.sim.PC.ARU.dat$db,
                         pix=NOBO.sim.PC.ARU.dat$pix,
                         point.area.PC=NOBO.sim.PC.ARU.dat$PointAreaPC.ha,
                         point.area.ARU=NOBO.sim.PC.ARU.dat$PointAreaARU.ha,
                         point.area.DIFF=NOBO.sim.PC.ARU.dat$PointAreaDiff.ha,
                         n.Grids=NOBO.sim.PC.ARU.dat$N.Grid,
                         n.PC=NOBO.sim.PC.ARU.dat$N.PC,
                         nvisits=NOBO.sim.PC.ARU.dat$Nvisits.PC,
                         GridwPC=NOBO.sim.PC.ARU.dat$Grid.w.PC.act[[x]])
            
            inits.fun <- function() list(b.0=runif(1,-2,2),
                                         beta=runif(1,4.9,5.5),
                                         N_ARU=rep((max(NOBO.sim.PC.ARU.dat$PC.N.det.data[[x]],na.rm=T)+1),times=NOBO.sim.PC.ARU.dat$N.Grid),
                                         N_diff=rep((max(NOBO.sim.PC.ARU.dat$PC.N.det.data[[x]],na.rm=T)+1),times=NOBO.sim.PC.ARU.dat$N.Grid))
            params <- c("N_PC","N_ARU","b.0","sigma")
            
            cl3 <- makePSOCKcluster(3)
            NOBO.model.cs <- jags.parfit(cl3, data=data, params=params,
                                         model="NOBO_C_sim.jags", inits=inits.fun(),
                                         n.chains=3, n.iter=120000, n.burnin=30000,n.adapt=30000,
                                         thin=25)
            
            stopCluster(cl3)
            
            NOBO.mcmc.samps <- rbind(NOBO.model.cs[[1]],NOBO.model.cs[[2]],NOBO.model.cs[[3]])
            N_PC.samps[[x]] <- NOBO.mcmc.samps[,grep("N_PC",colnames(NOBO.mcmc.samps))]
            N_ARU.samps[[x]] <- NOBO.mcmc.samps[,grep("N_ARU",colnames(NOBO.mcmc.samps))]
            b.0.samps[[x]] <- NOBO.mcmc.samps[,"b.0"]
            sigma.samps[[x]] <- NOBO.mcmc.samps[,"sigma"]
            GR[[x]] <- gelman.diag(NOBO.model.cs, multivariate=F)$psrf[,1]
            
          }
          
          rm(list=c("data","NOBO.model.cs","inits.fun","params"))
          
        }
        
        NOBO.sim.data.out <- list(Avg.Density=NOBO.sim.PC.ARU.dat$Avg.Density,
                                  n.Grid=NOBO.sim.PC.ARU.dat$N.Grid,
                                  Prop.Area=NOBO.sim.PC.ARU.dat$ProportionArea.PC,
                                  n.PC=NOBO.sim.PC.ARU.dat$N.PC,
                                  n.visits.PC=NOBO.sim.PC.ARU.dat$Nvisits.PC,
                                  n.ARU=NOBO.sim.PC.ARU.dat$N.ARUs,
                                  n.days.ARU=NOBO.sim.PC.ARU.dat$ARU.days,
                                  n.sims=NOBO.sim.PC.ARU.dat$Nsims,
                                  ARUarea=NOBO.sim.PC.ARU.dat$PointAreaARU.ha,
                                  PCarea=NOBO.sim.PC.ARU.dat$PointAreaPC.ha,
                                  avg.covey.size=NOBO.sim.PC.ARU.dat$avg.covey.size,
                                  N_PC.act=NOBO.sim.PC.ARU.dat$N.PC.act,
                                  N_ARU.act=NOBO.sim.PC.ARU.dat$N.ARU.act,
                                  sigma.act=NOBO.sim.PC.ARU.dat$PC.sigma.act,
                                  Grid.w.PC.act=NOBO.sim.PC.ARU.dat$Grid.w.PC.act,
                                  Grid.w.ARU.act=NOBO.sim.PC.ARU.dat$Grid.w.ARU.act,
                                  omega.act=NOBO.sim.PC.ARU.dat$omega.act,
                                  gamma.1.act=NOBO.sim.PC.ARU.dat$ARU.gamma.2.act,
                                  N_PC.obs=N_PC.samps,
                                  N_ARU.obs=N_ARU.samps,
                                  b.0.obs=b.0.samps,
                                  sigma.obs=sigma.samps,
                                  gamma.1.obs=gamma.1.samps,
                                  omega.obs=omega.samps,
                                  GR=GR)
        
        filename <- paste0("NOBO_simresults_",NOBO.sim.PC.ARU.dat$Avg.Density,"Dens_",NOBO.sim.PC.ARU.dat$N.Grid,"Grids_",NOBO.sim.PC.ARU.dat$N.PC,"PC_",NOBO.sim.PC.ARU.dat$Nvisits.PC,"visits_",NOBO.sim.PC.ARU.dat$N.ARUs,"ARUs_",NOBO.sim.PC.ARU.dat$ARU.days,"ARUdays.gzip")
        save(NOBO.sim.data.out, file=filename)
        
        rm(list=c("NOBO.sim.PC.ARU.dat","N_PC.samps","N_ARU.samps","b.0.samps",
                  "gamma.1.samps","sigma.samps","omega.samps","GR"))
        
      }
    }
  }
  
  