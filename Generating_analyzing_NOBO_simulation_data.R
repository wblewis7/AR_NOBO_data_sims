# Code for generating and analyzing simulation data using point counts and autonomous
#   recording units (ARUs) to estimate abundance of northern bobwhite (Colinus virginianus).
# From: Lewis, W. B., C. Johnson, and J. A. Martin. Assessing the efficacy and 
#   cost-effectiveness of integrating autonomous recording units and point-count
#   surveys for population monitoring of northern bobwhite (Colinus virginianus)
# Contact: William Lewis, wblewis7@gmail.com, University of Georgia


# For simulations, we imagined a hypothetical 2016-ha study area from which to estimate 
#   bobwhite density. We divided this area into 14 1200 x 1200 m grids at which
#   surveys could occur. 
# Simulations vary in average bobwhite density (0.1, 0.5, 0.9, 1.3, 1.7, 2.1, and 2.5 
#   birds/ha), three monitoring scenarios (model AV:only ARUs with manual validation, 
#   model C: only point counts, or model AVC: both), the proportion of grids surveyed
#   (22, 34, 45, 56% of survey area, corresponding to 4, 6, 8, or 10 locations), the
#   number of repeat point count surveys (2-4), and ARU deployment periods of either
#   14 or 28 days. We simulated data 100 times for each of the 308 scenarios.
# Simulation parameters are based on analysis of point count and ARU data on bobwhite 
#   abundance from Arkansas, USA, 2022.
# Analysis of the Arkansas data suggested that the number of detections across recordings
#   on the same ARU was highly variable, and this variation was poorly explained by estimates
#   of the parameter delta, which governs the calling rate and detection probability of 
#   bobwhite calls on ARUS. We therefore used a slightly different procedure between the 
#   data generation process and the statistical model for calculating delta. We believe 
#   that this method of data generation creates more biologically realistic data, thus 
#   better assessing the performance of the statistical model on ecological datasets;
#   however, we acknowledge that the simulations are unlikely to be able to recover
#   parameters associated with the ARU detection process. We therefore only assessed 
#   model performance of density estimates and the associated parameter β_0.

# We analyze simulated abundance data in a similar framework as in Analyzing_AR_NOBO_data.R. Data
#   are analyzed via Bayesian hierarchical models based on the models of Doser et al, 2021
#   (Integrating automated acoustic vocalization data and point count surveys for estimation
#   of bird abundance), with modifications by Nolan et al. 2024 (Effects of management practices
#   on Northern Bobwhite Colinus virginianus density in privately owned working forests across 
#   the Southeastern United States). Point count data are incorporated via distance-sampling
#   models modified for repeat surveys, while ARU data are incorporated via a zero-truncated
#   hurdle model to account for false positives and false negatives. The statistical model
#   varies depending on if points are surveyed only by point counts (model C), only ARUs with
#   manual validation of some detections (model AV), or both (model AVC).


require(truncnorm)
require(extraDistr)
require(dclone)
require(rjags)






# Generating simulation data --------------------------------------------------

  
  ## Simulation Parameters -----------------------------------------------------
  nsims <- 100 # Number of simulations
  AvgDens <- seq(0.1,2.5,by=0.4)
  len.AvgDens <- length(AvgDens)
  Ngrids <- c(4,6,8,10)
  ProportionArea <- round(Ngrids*600*600*pi/10000/2023.428,2) # 22%, 34%, 45% and 56% of area
  len.surveys <- length(Ngrids)
  Nvisits <- 2:4 # Number of repeat visits for point counts
  NdaysARU <- c(14,28) # Number of days of ARU deployment
  Propsurvey <- c(0,1)
  
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
  
  # Point counts truncated at 600m radius, preliminary analysis of playback of NOBO
  #   calls at known distance from ARU suggests that detection radius is 520m. This
  #   was determined by estimating sigma at ARU from playback and sigma at point
  #   counts through hierarchical distance-sampling model, and determining radius
  #   for ARUs at which detection probability equals point count detection probability
  #   at 600m.
  areaARU <- pi * 520 * 520
  areaPC <- pi * 600 * 600
  areaDiff <- areaPC - areaARU
  
  # Point count detection process
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
    pd_bins[b] <- (PC.sigma^2*(1-exp(-db[b+1]^2/(2*PC.sigma^2)))-PC.sigma^2*(1-exp(-db[b]^2/(2*PC.sigma^2))))*2*3.1416/(areaPC*pix[b])
  }
  
  # ARU hurdle model covariates
  alpha.0 <- -3.23 # Probability of detecting a false positive when N=0
  alpha.1 <- 1.3 # How probability of detecting at least one call changes with abundance
  
  # Putting covariates on delta for background noise. Simulating values each day.
  noise.mean <- 0
  noise.sd <- 1
  
  # Effect of standardized background noise on delta
  gamma.1 <- -0.18
  
  # Governs number of false positives on recordings
  omega <- 6.77 
  
  # Predicting number of calls on recordings. For unoccupied sites, all detections
  #   are false positives and so v.ARU will be based on omega. For occupied sites,
  #   predicting the probability that the recordings on a day will contain true 
  #   positives (delta.b0, delta.b1). If there are no true positives, v.ARU is 
  #   based on omega. If there are true positives, v.ARU is based on omega and delta 
  #   calculation below.
  delta.b0 <- -0.42 # Intercept (logit scale) for probability of at least one detection on a daily ARU recording at an occupied site
  delta.b1 <- 0.40 # Effect (logit scale) of abundance on probability of at least one detection on a daily ARU recording at an occupied site
  
  delta.mean <- 13.1 # Mean delta when at least one detection is recorded
  delta.sd <- 9.7 # Standard deviation of delta when at least one detection is recorded
  delta.min <- 0.7 # 5% quantile of delta when at least one detection is recorded
  delta.max <- 32.8 # 95% quantile of delta when at least one detection is recorded
  
  # Only a proportion of days for ARU deployment have suitable weather (not raining and wind speed <16kmph)
  p.suitweather <- 0.81 # From AR data
  
  
  ## Simulations ---------------------------------------------------------------

  for(z in 1:len.AvgDens){
    
    # Density is in birds/ha, but surveys assess covey density
    # Changing to coveys/ha by dividing by average covey size
    dens.c <- AvgDens[z]/avg.covey.size
    
    for(i in 1:len.surveys){
      for(j in 1:N.scenarios){
        
        PC.N.det.data <- PC.det.distbin.data <- PC.N.distbin.act <- ARU.noise.data <- ARU.y.data <- ARU.v.data <- which.v.ARU.data <- delta.act <- GridwARU.act <- A.times.data <- n.A.times.data <- GridwPC.act <- k.val.data <- n.val.data <- ARUID.val.data <- ARUday.val.data <- ARUID.sub.val.data <- ARU.Atimes.ID.val.data <- vector(mode="list", length=nsims)
        covey.N.ARU <- covey.N.PC <- matrix(NA, nrow=nsims, ncol=Ngrids[i])
        n.ARU.dates <- n.which.v.ARU.data <- PercentValid.data <- n.ARU.val.data <- rep(NA, times=nsims)
        
        for(x in 1:nsims){
          
          # Simulating abundance for ARUs and count areas. Ensuring that
          #    abundance on counts is equal or greater than ARUs
          covey.N.ARU[x,] <- rpois(Ngrids[i], areaARU/10000 * dens.c)
          covey.N.PC[x,] <- rpois(Ngrids[i], areaDiff/10000 * dens.c) + covey.N.ARU[x,]
          
          
          # Point counts (if doing)
          #######################################################################
          if(scenarios$PC[j]>0){
            
            nGridwPC <-scenarios$PC[j] * Ngrids[i]
            # Selecting which Grid locations get point counts
            GridwPC <- sample(1:Ngrids[i], nGridwPC, replace=F)
            GridwPC <- GridwPC[order(GridwPC)]
            
            # Allowing distance bin of coveys to vary by visit
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
          
          
          # ARUs (if doing)
          ########################################################################
          if(scenarios$ARU[j]>0){
            
            nGridwARU <-scenarios$ARU[j] * Ngrids[i]
            GridwARU <- sample(1:Ngrids[i], nGridwARU, replace=F)
            GridwARU <- GridwARU[order(GridwARU)]
            
            # Going to simulate days of suitable weather during the deployment period
            # Assume that this will be similar across days at the same site
            n.ARU.dates[x] <- rbinom(1,scenarios$ARUday[j],p.suitweather)
            
            # Simulating standardized background noise for each day of recording
            ARU.noise <- matrix(rnorm(nGridwARU * n.ARU.dates[x], noise.mean, noise.sd), nrow=nGridwARU)
            
            ARU.y <- delta <- TP <- ARU.v <- A.times <- matrix(NA, nrow=nGridwARU, ncol=n.ARU.dates[x])
            n.A.times <- rep(NA, times=nGridwARU)
            which.v.ARU <- c()
            
            for(t in 1:nGridwARU){
              
              # Probability of detecting at least one call on a recording, incorporates false positives, higher probability for higher abundance, and negative effect of background noise
              pa.temp <- plogis(alpha.0 + alpha.1*covey.N.ARU[x,GridwARU[t]])
              delta.dat <- c()
              
              for(q in 1:n.ARU.dates[x]){
                
                ARU.y[t,q] <- rbinom(1, 1, pa.temp)
                
                if(ARU.y[t,q] == 1){
                  
                  # For occupied sites, sometimes will only have false positives (number of calls based on omega)
                  #     and sometimes have actual calls (number of calls based on delta and omega)
                  delta.z <- ifelse(covey.N.ARU[x,GridwARU[t]] > 0,rbinom(1,1,plogis(delta.b0 + delta.b1*covey.N.ARU[x,GridwARU[t]])),0)
                  
                  # Delta (calls detected/bird), will be 0 if no calls detected or abundance is 0
                  # Adjusting mean based on effect of wind speed
                  delta.temp <- rtruncnorm(1, a=delta.min, b=delta.max, mean=delta.mean+gamma.1*ARU.noise[t,q], sd=delta.sd) * delta.z
                  delta.dat <- c(delta.dat, delta.temp)
                  
                  # ARU call detection data. Must be > 0
                  ARU.v[t,q] <- rtpois(1, delta.temp*covey.N.ARU[x,GridwARU[t]] + omega, a=0, b=Inf)
                  # True positives have to be 0 if delta.z = 0, otherwise must be at least 1
                  if(delta.z==0){
                    TP[t,q] <- 0
                  } else {
                    dtemp <- max(delta.temp*covey.N.ARU[x,GridwARU[t]] - 1, 0)
                    TP[t,q] <- 1 + rbinom(1, prob=dtemp/(dtemp + omega), ARU.v[t,q]-1)
                  }
                }
              }
              
              n.A.times[t] <- length(which(ARU.v[t,]>0))
              if(n.A.times[t]>0){
                A.times[t,1:n.A.times[t]] <- which(ARU.v[t,]>0)
                which.v.ARU <- c(which.v.ARU, t)
                delta[t,1:n.A.times[t]] <- delta.dat
              }
            }
            
            n.which.v.ARU <- length(which.v.ARU)
            delta <- delta[!is.na(delta[,1]),colSums(!apply(delta, 2, is.na)) > 0]
            
            
            # Validating false positives. Selecting up to 10 recordings, preferentially ones with true positives, then
            #   randomly selects 25% of recordings for validation. Preferentially selecting
            #   recordings with true positives because the Hypergeometric formulation 
            #   for incorporating manual validation estimates a true positive rate,
            #   which is only informative if true positives are present.
            
            k.val <- n.val <- ARUID.val <- ARUday.val <- ARUID.sub.val <- ARU.Atimes.ID.val <- c()
            n.ARU.val <- 0
            PercentValid <- NA
            
            if(sum(ARU.v, na.rm=T) > 0){
              
              k.val.tmp <- ARU.v.tmp <- ARUID.val.tmp <- ARUday.val.tmp <- ARUID.sub.val.tmp <- ARU.Atimes.ID.val.tmp <- c()
              vals.tmp <- vector(mode="list", length=0)
              
              for(t in 1:n.which.v.ARU){
                for(q in 1:n.A.times[which.v.ARU[t]]){
                  ARU.v.tmp <- c(ARU.v.tmp, ARU.v[which.v.ARU[t], A.times[which.v.ARU[t], q]])
                  k.val.tmp <- c(k.val.tmp, TP[which.v.ARU[t], A.times[which.v.ARU[t], q]])
                  vals.tmp <- append(vals.tmp, list(c(rep(1, times=TP[which.v.ARU[t], A.times[which.v.ARU[t], q]]), rep(0, times=ARU.v[which.v.ARU[t], A.times[which.v.ARU[t], q]] - TP[which.v.ARU[t], A.times[which.v.ARU[t], q]]))))
                  ARUID.val.tmp <- c(ARUID.val.tmp, which.v.ARU[t])
                  ARUday.val.tmp <- c(ARUday.val.tmp, A.times[which.v.ARU[t], q])
                }
              } 
              
              ARUID.sub.val.tmp <- ARU.Atimes.ID.val.tmp <- rep(NA, times=length(ARUID.val.tmp))
              for(r in 1:length(ARUID.sub.val.tmp)){
                ARUID.sub.val.tmp[r] <- which(which.v.ARU == ARUID.val.tmp[r])
                ARU.Atimes.ID.val.tmp[r] <- which(A.times[ARUID.val.tmp[r],1:n.A.times[ARUID.val.tmp[r]]] == ARUday.val.tmp[r])
              }
              
              n.to.val <- min(length(ARU.Atimes.ID.val.tmp), 10)
              
              # Going to preferentially select sites with true positives, barring those, will pick 
              #   from those with only false positives
              records.use <- which(k.val.tmp > 0)
              if(length(records.use) > n.to.val){
                records.use <- records.use[sample(x=1:length(records.use),size=n.to.val)]
              }
              if(length(records.use) < n.to.val){
                if(length(which(k.val.tmp == 0)) == 1){
                  records.use <- c(records.use, which(k.val.tmp == 0))
                } else{
                  records.use <- c(records.use, sample(x=which(k.val.tmp == 0), size=n.to.val - length(records.use)))
                }
              }
              k.val <- n.val <-rep(NA, times=length(records.use))
              # Going to select 25% of records from each, or maximum (if < 10)
              for(b in 1:length(records.use)){
                Nselect <- ceiling(length(vals.tmp[[records.use[b]]])*0.25)
                if(Nselect < 10){
                  Nselect <- ifelse(length(vals.tmp[[records.use[b]]]) > 10, 10, length(vals.tmp[[records.use[b]]]))
                }
                det.use <- sample(x=1:length(vals.tmp[[records.use[b]]]), size=Nselect)
                k.val[b] <- sum(vals.tmp[[records.use[b]]][det.use])
                n.val[b] <- Nselect
              }
              ARUID.val <- ARUID.val.tmp[records.use]
              ARUday.val <- ARUday.val.tmp[records.use]
              ARUID.sub.val <- ARUID.sub.val.tmp[records.use]
              ARU.Atimes.ID.val <- ARU.Atimes.ID.val.tmp[records.use]
              n.ARU.val <- length(ARUID.val)
              PercentValid <- sum(n.val)/sum(ARU.v, na.rm=T)
            } else {
              k.val <- n.val <- ARUID.val <- ARUday.val <- ARUID.sub.val <- ARU.Atimes.ID.val <- 0
            }
            
            
            if(n.which.v.ARU==0){
              ARU.v <- matrix(0, nrow=nrow(ARU.v), ncol=ncol(ARU.v))
              A.times <- matrix(0, nrow=nrow(A.times), ncol=ncol(A.times))
            }
            if(is.null(which.v.ARU)){
              which.v.ARU <- 0
            }
            
            GridwARU.act[[x]] <- GridwARU
            ARU.noise.data[[x]] <- ARU.noise
            ARU.y.data[[x]] <- ARU.y
            ARU.v.data[[x]] <- ARU.v
            A.times.data[[x]] <- A.times
            n.A.times.data[[x]] <- n.A.times
            n.which.v.ARU.data[x] <- n.which.v.ARU
            which.v.ARU.data[[x]] <- which.v.ARU
            delta.act[[x]] <- delta
            n.ARU.val.data[x] <- n.ARU.val
            PercentValid.data[x] <- PercentValid
            k.val.data[[x]] <- k.val
            n.val.data[[x]] <- n.val
            ARUID.val.data[[x]] <- ARUID.val
            ARUID.sub.val.data[[x]] <- ARUID.sub.val
            ARU.Atimes.ID.val.data[[x]] <- ARU.Atimes.ID.val
            ARUday.val.data[[x]] <- ARUday.val
          }  
        }
        
        sim.params <- list(nsims = nsims,
                           AvgDens = AvgDens[z],
                           CoveyDens = AvgDens[z]/avg.covey.size,
                           Nsites = Ngrids[i],
                           ProportionArea = ProportionArea[i],
                           Nvisits = scenarios$Visits[j],
                           NPC = scenarios$PC[j]*Ngrids[i],
                           PercentPC = scenarios$PC[j],
                           PercentARU = scenarios$ARU[j],
                           NARU = scenarios$ARU[j]*Ngrids[i],
                           NARUdays = scenarios$ARUday[j],
                           avg.covey.size = avg.covey.size,
                           areaARU = areaARU,
                           areaDiff = areaDiff,
                           areaPC = areaPC,
                           PC.sigma = PC.sigma,
                           pix = pix,
                           nB = nB,
                           db = db,
                           alpha.0 = alpha.0,
                           alpha.1 = alpha.1,
                           noise.mean = noise.mean,
                           noise.sd = noise.sd,
                           omega = omega,
                           delta.b0 = delta.b0,
                           delta.b1 = delta.b1,
                           delta.mean = delta.mean,
                           delta.sd = delta.sd,
                           delta.min = delta.min,
                           delta.max = delta.max,
                           gamma.1 = gamma.1,
                           p.suitweather = p.suitweather) 
        
        PC.data <- list(covey.N.PC = covey.N.PC,
                        GridwPC = GridwPC.act,
                        y.PC = PC.N.det.data,
                        ydb.PC = PC.det.distbin.data,
                        ydb.PC.act = PC.N.distbin.act)
        
        ARU.data <- list(covey.N.ARU = covey.N.ARU,
                         GridwARU.act = GridwARU.act,
                         n.ARU.dates = n.ARU.dates,
                         ARU.noise = ARU.noise.data,
                         y.ARU = ARU.y.data,
                         v.ARU = ARU.v.data,
                         A.times = A.times.data,
                         n.A.times = n.A.times.data,
                         n.which.v.ARU = n.which.v.ARU.data,
                         which.v.ARU = which.v.ARU.data,
                         delta.act = delta.act,
                         n.ARU.val = n.ARU.val.data,
                         PercentValid = PercentValid.data,
                         k.val = k.val.data,
                         n.val= n.val.data,
                         ARUID.val = ARUID.val.data,
                         ARUID.sub.val = ARUID.sub.val.data,
                         ARUday.val = ARUday.val.data,
                         ARU.Atimes.ID.val = ARU.Atimes.ID.val.data)
        
        # Saving out
        NOBO.sim.PC.ARU.dat <- list(sim.params = sim.params,
                                    PC.data = PC.data,
                                    ARU.data = ARU.data)
        
        filename <- paste0("NOBO_simdata_",AvgDens[z],"Dens_",Ngrids[i],"Grids_",scenarios$PC[j]*Ngrids[i],"PC_",scenarios$Visits[j],"visits_",scenarios$ARU[j]*Ngrids[i],"ARUs_",scenarios$ARUday[j],"ARUdays.gzip")
        save(NOBO.sim.PC.ARU.dat, file=filename)
      }
    }
  }
  
  
  
  
  
  
  
# Analyzing simulated data -----------------------------------------------------
  
  for(z in 1:len.AvgDens){
    for(i in 1:len.surveys){
      for(j in 1:N.scenarios){
        
        load(paste0("NOBO_simdata_",AvgDens[z],"Dens_",Ngrids[i],"Grids_",scenarios$PC[j]*Ngrids[i],"PC_",scenarios$Visits[j],"visits_",scenarios$ARU[j]*Ngrids[i],"ARUs_",scenarios$ARUday[j],"ARUdays.gzip"))
        
        N_PC.sum <- N_ARU.sum <-  b.0.samps <- GR <- vector(mode="list",length=NOBO.sim.PC.ARU.dat$sim.params$nsims)
        sigma.sum <- alpha.0.sum <- alpha.1.sum <- gamma.0.sum <- gamma.1.sum <- omega.sum <- matrix(NA, nrow=NOBO.sim.PC.ARU.dat$sim.params$nsims, 5)
        colnames(sigma.sum) <- colnames(alpha.0.sum) <- colnames(alpha.1.sum) <- colnames(gamma.0.sum) <- colnames(gamma.1.sum) <- colnames(omega.sum) <- c("Q2.5","Q50","Q97.5","Mean","SD")
        
        # Looping through simulations
        for(x in 1:NOBO.sim.PC.ARU.dat$sim.params$nsims){
          
          print(paste0("Doing simulation ",x,": ",NOBO.sim.PC.ARU.dat$sim.params$AvgDens,"Dens_",NOBO.sim.PC.ARU.dat$sim.params$ProportionArea*100,"%SURV_",ifelse(is.na(NOBO.sim.PC.ARU.dat$sim.params$Nvisits),0,NOBO.sim.PC.ARU.dat$sim.params$Nvisits),"PC.SURV_",ifelse(is.na(NOBO.sim.PC.ARU.dat$sim.params$NARUdays),0,NOBO.sim.PC.ARU.dat$sim.params$NARUdays),"ARU.SURV",": ",Sys.time()))
          
          
          if(NOBO.sim.PC.ARU.dat$sim.params$NARU > 0 & NOBO.sim.PC.ARU.dat$sim.params$NPC > 0){
            
            ## If have point counts and ARUs (model AVC) ----------------------------
            
            sink("NOBO_AVC_sim.jags")
              cat("
                model{
                  
                  # Priors
                  b.0 ~ dnorm(0, 0.001)
                  gamma.0 ~ dnorm(0, 0.001)
                  beta ~ dunif(0,10)
                  log(sigma) <- beta
                  mu.alpha ~ dunif(0,1)
                  alpha.0 <- logit(mu.alpha)
                  alpha.1 ~ dunif(0,1000)
                  gamma.1 ~ dunif(-1000,0)
                  omega ~ dunif(0.0001, 1000)
            
            
                  # Distance detection function for point counts
                  for(b in 1:nB){
                    pd_bins[b] <- (sigma^2*(1-exp(-db[b+1]^2/(2*sigma^2)))-sigma^2*(1-exp(-db[b]^2/(2*sigma^2))))*2*3.1416/(areaPC*pix[b])
                    pd_bins_adj[b] <- pd_bins[b]*pix[b]
                  }
                  p <- sum(pd_bins_adj[1:nB])
            
            
                  # Manual validation of a subset of ARU detections
                  for(v in 1:n.ARU.val){
                    tp[v] <- (delta[ARUID.sub.val[v], ARU.Atimes.ID.val[v]] * N_ARU[ARUID.val[v]]) / (delta[ARUID.sub.val[v], ARU.Atimes.ID.val[v]] * N_ARU[ARUID.val[v]] + omega)
                    K[v] ~ dbin(tp[v], v.ARU[ARUID.val[v], ARUday.val[v]])
                    k.val[v] ~ dhyper(K[v], v.ARU[ARUID.val[v], ARUday.val[v]] - K[v], n.val[v], 1)
                  }
            
            
                  ### State Model ----------------------------------------------
                  for(i in 1:npoints){
                    N_diff[i] ~ dpois(lambda_diff[i])
                    lambda_diff[i] <- exp(b.0 + log(areaDiff/10000))
                    N_ARU[i] ~ dpois(lambda_ARU[i])
                    lambda_ARU[i] <- exp(b.0 + log(areaARU/10000))
                    N_PC[i] <- N_ARU[i] + N_diff[i]
                    
                    
                    #### Point Count Observation Process -----------------------
                    for(j in 1:nvisits){
                      y.PC[i,j] ~ dbin(p, N_PC[i])
                      ydb.PC[i,1:nB,j] ~ dmulti(pd_bins_adj[1:nB],y.PC[i,j])
                    }
                    
                    
                    #### ARU Observtion Process: Hurdle Process ----------------
                    logit(p.a[i]) <- alpha.0 + alpha.1*N_ARU[i]
                    for(k in 1:n.ARU.days){
                      phi[i, k] ~ dgamma(1,1)
                      y.ARU[i, k] ~ dbin(p.a[i], 1)
                    }
                  }
                  
                  #### ARU Observation Process: Zero-truncated Poisson ---------
                  for(s in 1:n.which.v.ARU){
                    for(p in 1:n.A.times[which.v.ARU[s]]){
                      log(delta[s, p]) <- gamma.0 + gamma.1*ARU.noise[which.v.ARU[s], A.times[which.v.ARU[s], p]]
                      v.ARU[which.v.ARU[s], A.times[which.v.ARU[s], p]] ~ dpois((delta[s, p] * N_ARU[which.v.ARU[s]] + omega) * phi[which.v.ARU[s], A.times[which.v.ARU[s],p]] * y.ARU[which.v.ARU[s], A.times[which.v.ARU[s], p]])T(1, )
                    }
                  }
            
              }
            ",fill=TRUE)
          sink()
            
          data.mod.AVC <- list(y.PC = NOBO.sim.PC.ARU.dat$PC.data$y.PC[[x]],
                               ydb.PC = NOBO.sim.PC.ARU.dat$PC.data$ydb.PC[[x]],
                               nB = NOBO.sim.PC.ARU.dat$sim.params$nB,
                               db = NOBO.sim.PC.ARU.dat$sim.params$db,
                               pix = NOBO.sim.PC.ARU.dat$sim.params$pix,
                               areaPC = NOBO.sim.PC.ARU.dat$sim.params$areaPC,
                               areaARU = NOBO.sim.PC.ARU.dat$sim.params$areaARU,
                               areaDiff = NOBO.sim.PC.ARU.dat$sim.params$areaDiff,
                               n.ARU.val = NOBO.sim.PC.ARU.dat$ARU.data$n.ARU.val[[x]],
                               ARUID.sub.val = NOBO.sim.PC.ARU.dat$ARU.data$ARUID.sub.val[[x]],
                               ARU.Atimes.ID.val = NOBO.sim.PC.ARU.dat$ARU.data$ARU.Atimes.ID.val[[x]],
                               ARUID.val = NOBO.sim.PC.ARU.dat$ARU.data$ARUID.val[[x]],
                               ARUday.val = NOBO.sim.PC.ARU.dat$ARU.data$ARUday.val[[x]],
                               k.val = NOBO.sim.PC.ARU.dat$ARU.data$k.val[[x]],
                               n.val = NOBO.sim.PC.ARU.dat$ARU.data$n.val[[x]],
                               npoints = NOBO.sim.PC.ARU.dat$sim.params$Nsites,
                               nvisits = NOBO.sim.PC.ARU.dat$sim.params$Nvisits,
                               n.ARU.days = NOBO.sim.PC.ARU.dat$ARU.data$n.ARU.dates[[x]],
                               y.ARU = NOBO.sim.PC.ARU.dat$ARU.data$y.ARU[[x]],
                               n.which.v.ARU = NOBO.sim.PC.ARU.dat$ARU.data$n.which.v.ARU[x],
                               n.A.times = NOBO.sim.PC.ARU.dat$ARU.data$n.A.times[[x]],
                               which.v.ARU = NOBO.sim.PC.ARU.dat$ARU.data$which.v.ARU[[x]],
                               ARU.noise = NOBO.sim.PC.ARU.dat$ARU.data$ARU.noise[[x]],
                               A.times = NOBO.sim.PC.ARU.dat$ARU.data$A.times[[x]],
                               v.ARU = NOBO.sim.PC.ARU.dat$ARU.data$v.ARU[[x]])
            
            K.inits <- rep(NA, times=data.mod.AVC$n.ARU.val)
            for(m in 1:length(K.inits)){
              K.inits[m] <- max(data.mod.AVC$v.ARU[data.mod.AVC$ARUID.val[m], data.mod.AVC$ARUday.val[m]] - (data.mod.AVC$n.val[m] - data.mod.AVC$k.val[m]), data.mod.AVC$k.val[m])
            }
            inits.fun.mod.AVC <- function() list(b.0=runif(1,-2,2),
                                                 beta=runif(1,4.5,5.5),
                                                 N_ARU=rep((max(data.mod.AVC$y.PC,na.rm=T)+1),times=data.mod.AVC$npoints),
                                                 N_diff=rep((max(data.mod.AVC$y.PC,na.rm=T)+1),times=data.mod.AVC$npoints),
                                                 mu.alpha=runif(1,0.1,0.9),
                                                 alpha.1=runif(1,0,1),
                                                 gamma.0=runif(1,-1,1),
                                                 gamma.1=runif(1,-1,0),
                                                 omega=runif(1,0.1,3),
                                                 K=K.inits,
                                                 phi=matrix(runif(max(data.mod.AVC$n.ARU.days)*NOBO.sim.PC.ARU.dat$sim.params$NARU,0.9,1.1),nrow=NOBO.sim.PC.ARU.dat$sim.params$NARU))
            
            params.mod.AVC <- c("N_PC","N_ARU","b.0","alpha.0","alpha.1","gamma.0","gamma.1","sigma","omega")
            
            cl <- makePSOCKcluster(3)
            NOBO.model.cs.AVC <- jags.parfit(cl, data=data.mod.AVC, params=params.mod.AVC,
                                             model='NOBO_AVC_sim.jags', inits=inits.fun.mod.AVC,
                                             n.chains=3, n.iter=200000, n.burnin=60000,n.adapt=60000,
                                             thin=50)
            
            stopCluster(cl)
            

            # Only saving out samples for b.0. For other parameters, only doing the summaries to save memory.
            b.0.samps[[x]] <- do.call(c,NOBO.model.cs.AVC[,"b.0"])
            
            alpha.0.sum[x,] <- c(quantile(do.call(c,NOBO.model.cs.AVC[,"alpha.0"]), probs=c(0.025,0.5,0.975)),
                                      mean(do.call(c,NOBO.model.cs.AVC[,"alpha.0"])),
                                      sd(do.call(c,NOBO.model.cs.AVC[,"alpha.0"])))
            alpha.1.sum[x,] <- c(quantile(do.call(c,NOBO.model.cs.AVC[,"alpha.1"]), probs=c(0.025,0.5,0.975)),
                                      mean(do.call(c,NOBO.model.cs.AVC[,"alpha.1"])),
                                      sd(do.call(c,NOBO.model.cs.AVC[,"alpha.1"])))
            gamma.0.sum[x,] <- c(quantile(do.call(c,NOBO.model.cs.AVC[,"gamma.0"]), probs=c(0.025,0.5,0.975)),
                                      mean(do.call(c,NOBO.model.cs.AVC[,"gamma.0"])),
                                      sd(do.call(c,NOBO.model.cs.AVC[,"gamma.0"])))
            gamma.1.sum[x,] <- c(quantile(do.call(c,NOBO.model.cs.AVC[,"gamma.1"]), probs=c(0.025,0.5,0.975)),
                                      mean(do.call(c,NOBO.model.cs.AVC[,"gamma.1"])),
                                      sd(do.call(c,NOBO.model.cs.AVC[,"gamma.1"])))
            sigma.sum[x,] <- c(quantile(do.call(c,NOBO.model.cs.AVC[,"sigma"]), probs=c(0.025,0.5,0.975)),
                                    mean(do.call(c,NOBO.model.cs.AVC[,"sigma"])),
                                    sd(do.call(c,NOBO.model.cs.AVC[,"sigma"])))
            omega.sum[x,] <- c(quantile(do.call(c,NOBO.model.cs.AVC[,"omega"]), probs=c(0.025,0.5,0.975)),
                                    mean(do.call(c,NOBO.model.cs.AVC[,"omega"])),
                                    sd(do.call(c,NOBO.model.cs.AVC[,"omega"])))
            
            
            N_PC.sum[[x]] <- rbind(apply(do.call(rbind, NOBO.model.cs.AVC[,paste0("N_PC[",1:data.mod.AVC$npoints,"]")]), 2, quantile, probs=c(0.025,0.5,0.975)),
                                        apply(do.call(rbind, NOBO.model.cs.AVC[,paste0("N_PC[",1:data.mod.AVC$npoints,"]")]), 2, mean),
                                        apply(do.call(rbind, NOBO.model.cs.AVC[,paste0("N_PC[",1:data.mod.AVC$npoints,"]")]), 2, sd))
            N_ARU.sum[[x]] <- rbind(apply(do.call(rbind, NOBO.model.cs.AVC[,paste0("N_ARU[",1:data.mod.AVC$npoints,"]")]), 2, quantile, probs=c(0.025,0.5,0.975)),
                                         apply(do.call(rbind, NOBO.model.cs.AVC[,paste0("N_ARU[",1:data.mod.AVC$npoints,"]")]), 2, mean),
                                         apply(do.call(rbind, NOBO.model.cs.AVC[,paste0("N_ARU[",1:data.mod.AVC$npoints,"]")]), 2, sd))
            rownames(N_PC.sum[[x]]) <- rownames(N_ARU.sum[[x]]) <- c("Q2.5","Q50","Q97.5","Mean","SD")
            
            GR[[x]] <- gelman.diag(NOBO.model.cs.AVC, multivariate=F)$psrf
            
            rm(list=c("NOBO.model.cs.AVC", "params.mod.AVC", "inits.fun.mod.AVC", "K.inits", "data.mod.AVC"))
            
            
          } else if(NOBO.sim.PC.ARU.dat$sim.params$NPC > 0){
            
            ## If have only point counts (model C) -----------------------------
            
            sink("NOBO_C_sim.jags")
            cat("
                model{
                  
                  # Priors
                  b.0 ~ dnorm(0, 0.001)
                  gamma.0 ~ dnorm(0, 0.001)
                  beta ~ dunif(0,10)
                  log(sigma) <- beta
                  mu.alpha ~ dunif(0,1)
                  alpha.0 <- logit(mu.alpha)
                  alpha.1 ~ dunif(0,1000)
                  gamma.1 ~ dunif(-1000,0)
                  omega ~ dunif(0.0001, 1000)
            
            
                  # Distance detection function for point counts
                  for(b in 1:nB){
                    pd_bins[b] <- (sigma^2*(1-exp(-db[b+1]^2/(2*sigma^2)))-sigma^2*(1-exp(-db[b]^2/(2*sigma^2))))*2*3.1416/(areaPC*pix[b])
                    pd_bins_adj[b] <- pd_bins[b]*pix[b]
                  }
                  p <- sum(pd_bins_adj[1:nB])
            
            
                  ### State Model ----------------------------------------------
                  for(i in 1:npoints){
                    N_diff[i] ~ dpois(lambda_diff[i])
                    lambda_diff[i] <- exp(b.0 + log(areaDiff/10000))
                    N_ARU[i] ~ dpois(lambda_ARU[i])
                    lambda_ARU[i] <- exp(b.0 + log(areaARU/10000))
                    N_PC[i] <- N_ARU[i] + N_diff[i]
                    
                    
                    #### Point Count Observation Process -----------------------
                    for(j in 1:nvisits){
                      y.PC[i,j] ~ dbin(p, N_PC[i])
                      ydb.PC[i,1:nB,j] ~ dmulti(pd_bins_adj[1:nB],y.PC[i,j])
                    }
                  }
              }
            ",fill=TRUE)
            sink()
            
            data.mod.C <- list(y.PC = NOBO.sim.PC.ARU.dat$PC.data$y.PC[[x]],
                                 ydb.PC = NOBO.sim.PC.ARU.dat$PC.data$ydb.PC[[x]],
                                 nB = NOBO.sim.PC.ARU.dat$sim.params$nB,
                                 db = NOBO.sim.PC.ARU.dat$sim.params$db,
                                 pix = NOBO.sim.PC.ARU.dat$sim.params$pix,
                                 areaPC = NOBO.sim.PC.ARU.dat$sim.params$areaPC,
                                 areaARU = NOBO.sim.PC.ARU.dat$sim.params$areaARU,
                                 areaDiff = NOBO.sim.PC.ARU.dat$sim.params$areaDiff,
                                 npoints = NOBO.sim.PC.ARU.dat$sim.params$Nsites,
                                 nvisits = NOBO.sim.PC.ARU.dat$sim.params$Nvisits)

            inits.fun.mod.C <- function() list(b.0=runif(1,-2,2),
                                                 beta=runif(1,4.5,5.5),
                                                 N_ARU=rep((max(data.mod.C$y.PC,na.rm=T)+1),times=data.mod.C$npoints),
                                                 N_diff=rep((max(data.mod.C$y.PC,na.rm=T)+1),times=data.mod.C$npoints),
                                                 mu.alpha=runif(1,0.1,0.9),
                                                 alpha.1=runif(1,0,1),
                                                 gamma.0=runif(1,-1,1),
                                                 gamma.1=runif(1,-1,0),
                                                 omega=runif(1,0.1,3))
            
            params.mod.C <- c("N_PC","N_ARU","b.0","alpha.0","alpha.1","gamma.0","gamma.1","sigma","omega")
            
            cl <- makePSOCKcluster(3)
            NOBO.model.cs.C <- jags.parfit(cl, data=data.mod.C, params=params.mod.C,
                                             model='NOBO_C_sim.jags', inits=inits.fun.mod.C,
                                             n.chains=3, n.iter=200000, n.burnin=60000,n.adapt=60000,
                                             thin=50)
            
            stopCluster(cl)
            
            
            # Only saving out samples for b.0. For other parameters, only doing the summaries to save memory.
            # Posterior samples from the ARU observation process (all but b.0 and sigma) are not informed
            #   by data and are just based on priors.
            b.0.samps[[x]] <- do.call(c,NOBO.model.cs.C[,"b.0"])
            
            alpha.0.sum[x,] <- c(quantile(do.call(c,NOBO.model.cs.C[,"alpha.0"]), probs=c(0.025,0.5,0.975)),
                                 mean(do.call(c,NOBO.model.cs.C[,"alpha.0"])),
                                 sd(do.call(c,NOBO.model.cs.C[,"alpha.0"])))
            alpha.1.sum[x,] <- c(quantile(do.call(c,NOBO.model.cs.C[,"alpha.1"]), probs=c(0.025,0.5,0.975)),
                                 mean(do.call(c,NOBO.model.cs.C[,"alpha.1"])),
                                 sd(do.call(c,NOBO.model.cs.C[,"alpha.1"])))
            gamma.0.sum[x,] <- c(quantile(do.call(c,NOBO.model.cs.C[,"gamma.0"]), probs=c(0.025,0.5,0.975)),
                                 mean(do.call(c,NOBO.model.cs.C[,"gamma.0"])),
                                 sd(do.call(c,NOBO.model.cs.C[,"gamma.0"])))
            gamma.1.sum[x,] <- c(quantile(do.call(c,NOBO.model.cs.C[,"gamma.1"]), probs=c(0.025,0.5,0.975)),
                                 mean(do.call(c,NOBO.model.cs.C[,"gamma.1"])),
                                 sd(do.call(c,NOBO.model.cs.C[,"gamma.1"])))
            sigma.sum[x,] <- c(quantile(do.call(c,NOBO.model.cs.C[,"sigma"]), probs=c(0.025,0.5,0.975)),
                               mean(do.call(c,NOBO.model.cs.C[,"sigma"])),
                               sd(do.call(c,NOBO.model.cs.C[,"sigma"])))
            omega.sum[x,] <- c(quantile(do.call(c,NOBO.model.cs.C[,"omega"]), probs=c(0.025,0.5,0.975)),
                               mean(do.call(c,NOBO.model.cs.C[,"omega"])),
                               sd(do.call(c,NOBO.model.cs.C[,"omega"])))
            
            
            N_PC.sum[[x]] <- rbind(apply(do.call(rbind, NOBO.model.cs.C[,paste0("N_PC[",1:data.mod.C$npoints,"]")]), 2, quantile, probs=c(0.025,0.5,0.975)),
                                   apply(do.call(rbind, NOBO.model.cs.C[,paste0("N_PC[",1:data.mod.C$npoints,"]")]), 2, mean),
                                   apply(do.call(rbind, NOBO.model.cs.C[,paste0("N_PC[",1:data.mod.C$npoints,"]")]), 2, sd))
            N_ARU.sum[[x]] <- rbind(apply(do.call(rbind, NOBO.model.cs.C[,paste0("N_ARU[",1:data.mod.C$npoints,"]")]), 2, quantile, probs=c(0.025,0.5,0.975)),
                                    apply(do.call(rbind, NOBO.model.cs.C[,paste0("N_ARU[",1:data.mod.C$npoints,"]")]), 2, mean),
                                    apply(do.call(rbind, NOBO.model.cs.C[,paste0("N_ARU[",1:data.mod.C$npoints,"]")]), 2, sd))
            rownames(N_PC.sum[[x]]) <- rownames(N_ARU.sum[[x]]) <- c("Q2.5","Q50","Q97.5","Mean","SD")
            
            GR[[x]] <- gelman.diag(NOBO.model.cs.C, multivariate=F)$psrf
            
            rm(list=c("NOBO.model.cs.C", "params.mod.C", "inits.fun.mod.C", "data.mod.C"))
            
            
          } else{
            
            ## If have only ARUs (model AV) -----------------------------
            
            sink("NOBO_AV_sim.jags")
            cat("
                model{
                  
                  # Priors
                  b.0 ~ dnorm(0, 0.001)
                  gamma.0 ~ dnorm(0, 0.001)
                  beta ~ dunif(0,10)
                  log(sigma) <- beta
                  mu.alpha ~ dunif(0,1)
                  alpha.0 <- logit(mu.alpha)
                  alpha.1 ~ dunif(0,1000)
                  gamma.1 ~ dunif(-1000,0)
                  omega ~ dunif(0.0001, 1000)
            

                  # Manual validation of a subset of ARU detections
                  for(v in 1:n.ARU.val){
                    tp[v] <- (delta[ARUID.sub.val[v], ARU.Atimes.ID.val[v]] * N_ARU[ARUID.val[v]]) / (delta[ARUID.sub.val[v], ARU.Atimes.ID.val[v]] * N_ARU[ARUID.val[v]] + omega)
                    K[v] ~ dbin(tp[v], v.ARU[ARUID.val[v], ARUday.val[v]])
                    k.val[v] ~ dhyper(K[v], v.ARU[ARUID.val[v], ARUday.val[v]] - K[v], n.val[v], 1)
                  }
            
            
                  ### State Model ----------------------------------------------
                  for(i in 1:npoints){
                    N_diff[i] ~ dpois(lambda_diff[i])
                    lambda_diff[i] <- exp(b.0 + log(areaDiff/10000))
                    N_ARU[i] ~ dpois(lambda_ARU[i])
                    lambda_ARU[i] <- exp(b.0 + log(areaARU/10000))
                    N_PC[i] <- N_ARU[i] + N_diff[i]
                    
                    
                    #### ARU Observtion Process: Hurdle Process ----------------
                    logit(p.a[i]) <- alpha.0 + alpha.1*N_ARU[i]
                    for(k in 1:n.ARU.days){
                      phi[i, k] ~ dgamma(1,1)
                      y.ARU[i, k] ~ dbin(p.a[i], 1)
                    }
                  }
                  
                  #### ARU Observation Process: Zero-truncated Poisson ---------
                  for(s in 1:n.which.v.ARU){
                    for(p in 1:n.A.times[which.v.ARU[s]]){
                      log(delta[s, p]) <- gamma.0 + gamma.1*ARU.noise[which.v.ARU[s], A.times[which.v.ARU[s], p]]
                      v.ARU[which.v.ARU[s], A.times[which.v.ARU[s], p]] ~ dpois((delta[s, p] * N_ARU[which.v.ARU[s]] + omega) * phi[which.v.ARU[s], A.times[which.v.ARU[s],p]] * y.ARU[which.v.ARU[s], A.times[which.v.ARU[s], p]])T(1, )
                    }
                  }
            
              }
            ",fill=TRUE)
            sink()
            
            data.mod.AV <- list(areaARU = NOBO.sim.PC.ARU.dat$sim.params$areaARU,
                                 areaDiff = NOBO.sim.PC.ARU.dat$sim.params$areaDiff,
                                 n.ARU.val = NOBO.sim.PC.ARU.dat$ARU.data$n.ARU.val[[x]],
                                 ARUID.sub.val = NOBO.sim.PC.ARU.dat$ARU.data$ARUID.sub.val[[x]],
                                 ARU.Atimes.ID.val = NOBO.sim.PC.ARU.dat$ARU.data$ARU.Atimes.ID.val[[x]],
                                 ARUID.val = NOBO.sim.PC.ARU.dat$ARU.data$ARUID.val[[x]],
                                 ARUday.val = NOBO.sim.PC.ARU.dat$ARU.data$ARUday.val[[x]],
                                 k.val = NOBO.sim.PC.ARU.dat$ARU.data$k.val[[x]],
                                 n.val = NOBO.sim.PC.ARU.dat$ARU.data$n.val[[x]],
                                 npoints = NOBO.sim.PC.ARU.dat$sim.params$Nsites,
                                 n.ARU.days = NOBO.sim.PC.ARU.dat$ARU.data$n.ARU.dates[[x]],
                                 y.ARU = NOBO.sim.PC.ARU.dat$ARU.data$y.ARU[[x]],
                                 n.which.v.ARU = NOBO.sim.PC.ARU.dat$ARU.data$n.which.v.ARU[x],
                                 n.A.times = NOBO.sim.PC.ARU.dat$ARU.data$n.A.times[[x]],
                                 which.v.ARU = NOBO.sim.PC.ARU.dat$ARU.data$which.v.ARU[[x]],
                                 ARU.noise = NOBO.sim.PC.ARU.dat$ARU.data$ARU.noise[[x]],
                                 A.times = NOBO.sim.PC.ARU.dat$ARU.data$A.times[[x]],
                                 v.ARU = NOBO.sim.PC.ARU.dat$ARU.data$v.ARU[[x]])
            
            K.inits <- rep(NA, times=data.mod.AV$n.ARU.val)
            for(m in 1:length(K.inits)){
              K.inits[m] <- max(data.mod.AV$v.ARU[data.mod.AV$ARUID.val[m], data.mod.AV$ARUday.val[m]] - (data.mod.AV$n.val[m] - data.mod.AV$k.val[m]), data.mod.AV$k.val[m])
            }
            inits.fun.mod.AV <- function() list(b.0=runif(1,-2,2),
                                                 beta=runif(1,4.5,5.5),
                                                 N_ARU=sample(1:20,NOBO.sim.PC.ARU.dat$sim.params$Nsites,replace=T),
                                                 N_diff=sample(1:20,NOBO.sim.PC.ARU.dat$sim.params$Nsites,replace=T),
                                                 mu.alpha=runif(1,0.1,0.9),
                                                 alpha.1=runif(1,0,1),
                                                 gamma.0=runif(1,-1,1),
                                                 gamma.1=runif(1,-1,0),
                                                 omega=runif(1,0.1,3),
                                                 K=K.inits,
                                                 phi=matrix(runif(max(data.mod.AV$n.ARU.days)*NOBO.sim.PC.ARU.dat$sim.params$NARU,0.9,1.1),nrow=NOBO.sim.PC.ARU.dat$sim.params$NARU))
            
            params.mod.AV <- c("N_PC","N_ARU","b.0","alpha.0","alpha.1","gamma.0","gamma.1","sigma","omega")
            
            cl <- makePSOCKcluster(3)
            NOBO.model.cs.AV <- jags.parfit(cl, data=data.mod.AV, params=params.mod.AV,
                                             model='NOBO_AV_sim.jags', inits=inits.fun.mod.AV,
                                             n.chains=3, n.iter=200000, n.burnin=60000,n.adapt=60000,
                                             thin=50)
            
            stopCluster(cl)
            
            
            # Only saving out samples for b.0. For other parameters, only doing the summaries to save memory.
            # Posterior samples from point count observation process (sigma) are not informed
            #   by data and are just based on priors. N_PC can be estimated
            #   from b0 but is not directly related to data.
            b.0.samps[[x]] <- do.call(c,NOBO.model.cs.AV[,"b.0"])
            
            alpha.0.sum[x,] <- c(quantile(do.call(c,NOBO.model.cs.AV[,"alpha.0"]), probs=c(0.025,0.5,0.975)),
                                 mean(do.call(c,NOBO.model.cs.AV[,"alpha.0"])),
                                 sd(do.call(c,NOBO.model.cs.AV[,"alpha.0"])))
            alpha.1.sum[x,] <- c(quantile(do.call(c,NOBO.model.cs.AV[,"alpha.1"]), probs=c(0.025,0.5,0.975)),
                                 mean(do.call(c,NOBO.model.cs.AV[,"alpha.1"])),
                                 sd(do.call(c,NOBO.model.cs.AV[,"alpha.1"])))
            gamma.0.sum[x,] <- c(quantile(do.call(c,NOBO.model.cs.AV[,"gamma.0"]), probs=c(0.025,0.5,0.975)),
                                 mean(do.call(c,NOBO.model.cs.AV[,"gamma.0"])),
                                 sd(do.call(c,NOBO.model.cs.AV[,"gamma.0"])))
            gamma.1.sum[x,] <- c(quantile(do.call(c,NOBO.model.cs.AV[,"gamma.1"]), probs=c(0.025,0.5,0.975)),
                                 mean(do.call(c,NOBO.model.cs.AV[,"gamma.1"])),
                                 sd(do.call(c,NOBO.model.cs.AV[,"gamma.1"])))
            sigma.sum[x,] <- c(quantile(do.call(c,NOBO.model.cs.AV[,"sigma"]), probs=c(0.025,0.5,0.975)),
                               mean(do.call(c,NOBO.model.cs.AV[,"sigma"])),
                               sd(do.call(c,NOBO.model.cs.AV[,"sigma"])))
            omega.sum[x,] <- c(quantile(do.call(c,NOBO.model.cs.AV[,"omega"]), probs=c(0.025,0.5,0.975)),
                               mean(do.call(c,NOBO.model.cs.AV[,"omega"])),
                               sd(do.call(c,NOBO.model.cs.AV[,"omega"])))
            
            
            N_PC.sum[[x]] <- rbind(apply(do.call(rbind, NOBO.model.cs.AV[,paste0("N_PC[",1:data.mod.AV$npoints,"]")]), 2, quantile, probs=c(0.025,0.5,0.975)),
                                   apply(do.call(rbind, NOBO.model.cs.AV[,paste0("N_PC[",1:data.mod.AV$npoints,"]")]), 2, mean),
                                   apply(do.call(rbind, NOBO.model.cs.AV[,paste0("N_PC[",1:data.mod.AV$npoints,"]")]), 2, sd))
            N_ARU.sum[[x]] <- rbind(apply(do.call(rbind, NOBO.model.cs.AV[,paste0("N_ARU[",1:data.mod.AV$npoints,"]")]), 2, quantile, probs=c(0.025,0.5,0.975)),
                                    apply(do.call(rbind, NOBO.model.cs.AV[,paste0("N_ARU[",1:data.mod.AV$npoints,"]")]), 2, mean),
                                    apply(do.call(rbind, NOBO.model.cs.AV[,paste0("N_ARU[",1:data.mod.AV$npoints,"]")]), 2, sd))
            rownames(N_PC.sum[[x]]) <- rownames(N_ARU.sum[[x]]) <- c("Q2.5","Q50","Q97.5","Mean","SD")
            
            GR[[x]] <- gelman.diag(NOBO.model.cs.AV, multivariate=F)$psrf
            
            rm(list=c("NOBO.model.cs.AV", "params.mod.AV", "inits.fun.mod.AV", "K.inits", "data.mod.AV"))
            
          }
        }
        
        NOBO.sim.sum <- list(N_PC.sum. = N_PC.sum,
                             N_ARU.sum = N_PC.sum,
                             b.0.samps = b.0.samps,
                             sigma.sum = sigma.sum,
                             alpha.0.sum = alpha.0.sum,
                             alpha.1.sum = alpha.1.sum,
                             gamma.0.sum = gamma.0.sum,
                             gamma.1.sum = gamma.1.sum,
                             omega.sum = omega.sum,
                             GR = GR)
        
        NOBO.sim.out <- list(sim.params = NOBO.sim.PC.ARU.dat$sim.params,
                             NOBO.sim.sum = NOBO.sim.sum)
        
        filename <- paste0("NOBO_simresults_",AvgDens[z],"Dens_",Ngrids[i],"Grids_",scenarios$PC[j]*Ngrids[i],"PC_",scenarios$Visits[j],"visits_",scenarios$ARU[j]*Ngrids[i],"ARUs_",scenarios$ARUday[j],"ARUdays.gzip")
        save(NOBO.sim.out, file=filename)
        
      }
    }
  }
