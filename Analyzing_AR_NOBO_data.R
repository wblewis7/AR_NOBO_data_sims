require(rjags)

load("AR_NOBO_data.gzip")

# Model
sink("NOBO_AR_integrated.jags")
cat("
    model{
    
    b.0 ~ dnorm(0, 0.001)
    for(a in 1:nareaeffects){
      b.areas[a] ~ dnorm(0, 0.001)
    }
    beta ~ dunif(0,10)
    log(sigma) <- beta
    
    a.phi ~ dunif(0,100)
    
    # alpha.0 and alpha.1 are probability of detections on an unoccupied site, and additional probaility of detections
    #    with birds at site. Constraining alpha.1 to be positive, since can't have lower chance of detecting a call
    #    when birds are present
    mu.alpha ~ dunif(0,1)
    alpha.0 <- logit(mu.alpha)
    alpha.1 ~ dunif(0,1000)

    
    # Each area gets different intercept for delta. Intercept is first area
    gamma.0 ~ dnorm(0, 0.001)
    for(a in 1:nareaeffects){
      gamma.1[a] ~ dnorm(0, 0.001)
    }
    # Effect of background noise on delta. Expect negative effect
    gamma.2 ~ dunif(-1000,0)
    
   
    # Omega is the mean number of false positive acoustic detections. Modeling with poisson distribution
    omega.log ~ dunif(0,10)
    log(omega) <- omega.log
    for(x in 1:nvalidate){
      FP[x] ~ dpois(omega)
    }
      
    for(b in 1:nB){
        pd_bins[b] <- (sigma^2*(1-exp(-db[b+1]^2/(2*sigma^2)))-sigma^2*(1-exp(-db[b]^2/(2*sigma^2))))*2*3.1416/(point.area.PC*pix[b])
        pd_bins_adj[b] <- pd_bins[b]*pix[b]
      }
    pd <- sum(pd_bins_adj[1:nB])
    
    
    
    
    
    
    for(i in 1:npoints){


      ############################################################################
      ## State Model
      ############################################################################
    
      #Need to model differences between ARU and point count detection area, point count
      #   area is larger than ARU, so abundance on point counts must be larger than or equal to
      #   abundance on ARUs. Including offset of area
    
      N_diff[i] ~ dpois(lambda_diff[i]) # Number of birds in area between ARU area and point count area
      lambda_diff[i] <- exp(b.0 + sum(b.areas[1:nareaeffects]*area[i,1:nareaeffects]) + log(point.area.DIFF/10000))
      
      N_ARU[i] ~ dpois(lambda_ARU[i]) # Number of birds in area surveyed by both ARU and point counts
      lambda_ARU[i] <- exp(b.0 + sum(b.areas[1:nareaeffects]*area[i,1:nareaeffects]) + log(point.area.ARU/10000))
      
      N_PC[i] <- N_ARU[i] + N_diff[i]
      
      
      
      
      ##########################################################################
      ## Observation Process
      ##########################################################################
      
      
      ##### Point Counts ######################################################
      
      for(j in 1:nvisits[i]){
      
          # Number detected
          
          y.PC[i,j] ~ dbin(pd, N_PC[i])
      
          # Bin of detection
          
          ydb.PC[i,1:nB,j] ~ dmulti(pd_bins_adj[1:nB],y.PC[i,j])
      }
    }
    
    
    
    ###### ARUs ###############################################################
    
    for(g in 1:n.ARU){
    
      # Probability of detecting at least one detection on a recording
      logit(p.a[g]) <- alpha.0 + alpha.1*N_ARU[ARU.pointID[g]]
      
      for(r in 1:maxdates){
        # Estimating overdispersion parameter for each visit based on a hyperparameter
        phi[g, r] ~ dgamma(a.phi, a.phi)
      }
      
      for(k in 1:n.ARU.dates[g]){
      
        # Presence/absence of detections at ARUs
        
        y.ARU[g, k] ~ dbin(p.a[g], 1)
        
      }
      
      for(p in 1:n.A.times[g]){
      
        # Each recording gets its own delta, based on area and background noise
        log(delta[g, p]) <- gamma.0 + sum(gamma.1[1:nareaeffects]*area[ARU.pointID[g],1:nareaeffects]) + gamma.2*ARU.noise[g, A.times[g, p]]
      
        
      
        # The number of detections, conditional on at least one being detected, are a function of 
        #     delta (calling rate/detectability), abundance, omega (false positives), 
        #     and presence/absence of detections
        # Adding in overdispersion parameter (phi) to aid convergence
        # Truncating to minimum of 1 to prevent observations with y.ARU=1 from having v.ARU=0
        
        v.ARU[g, A.times[g, p]] ~ dpois((delta[g, p] * N_ARU[ARU.pointID[g]] + omega) * phi[g, A.times[g,p]] * y.ARU[g, A.times[g, p]])T(1, )
      }
    } 
    
    }
    ",fill=TRUE)
sink()


data <- AR_NOBO_data


inits.fun <- function() list(b.0=runif(1,-2,2),
                             b.areas=runif(data$nareaeffects,-2,2),
                             beta=runif(1,4.5,5.5),
                             N_ARU=rep((max(data$y.PC,na.rm=T)+1),times=data$npoint),
                             N_diff=rep((max(data$ydb.PC,na.rm=T)+1),times=data$npoint),
                             mu.alpha=runif(1,0.1,0.9),
                             alpha.1=runif(1,0,1),
                             gamma.0=runif(1,-2,2),
                             gamma.1=runif(data$nareaeffects,-2,2),
                             gamma.2=runif(1,-10,0),
                             omega.log=runif(1,0,3),
                             phi=matrix(runif(max(data$n.ARU.dates)*data$n.ARU,0.01,2),nrow=data$n.ARU))


params <- c("b.0","b.areas","N_ARU","N_PC","sigma","pd","delta","gamma.0",
            "gamma.1","gamma.2","omega","p.a","alpha.0","alpha.1","a.phi","phi")

NOBO.AR.model <- jags.model(data=data, inits = inits.fun, 
                                   file='NOBO_AR_integrated_model.jags',
                                   n.chain = 3, n.adapt = 10000)
NOBO.AR.model.cs <- coda.samples(NOBO.AR.model, variable.names=params, n.iter=100000, n.burn=10000, thin=10)
