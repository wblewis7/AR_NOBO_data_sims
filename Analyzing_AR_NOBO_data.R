# Sample code for integrating point counts and autonomous recording units (ARUs)
#   to estimate abundance of northern bobwhite (Colinus virginianus) in Arkansas,
#   USA, 2022.
# From: Lewis, W. B., C. Johnson, and J. A. Martin. Assessing the efficacy and 
#   cost-effectiveness of integrating autonomous recording units and point-count
#   surveys for population monitoring of northern bobwhite (Colinus virginianus)
# Contact: William Lewis, wblewis7@gmail.com, University of Georgia


# Coveys (groups of bobwhites) were surveyed during the fall of 2022 using traditional
#   distance-sampling point counts at 29 survey locations across 8 study sites in 
#   Arkansas, USA. Point count surveys were repeated 1 - 3 times at each location, and 
#   distance to observer values were grouped into 6 distance bins. ARUs were deployed
#   at 26 of these sites, set to record once per day from ~ October 1 - November 10.
#   Bobwhite covey calls were detected on ARU recordings using a Convolutional Neural
#   Network automated classifier (Nolan et al. 2023. The development of a convolutional
#   neural network for the automatic detection of Northern Bobwhite Colinus virginianus
#   covey calls) and a score threshold of 0.95.

# Point count and ARU data are integrated in a Bayesian hierarchical framework. Both
#   datasets jointly estimate the state (abundance) process but have separate observation
#   processes. We use a modified version of the statistical model of Doser et al. 2021
#   (Integrating automated acoustic vocalization data and point count surveys for estimation
#   of bird abundance), with modifications by Nolan et al. 2024 (Effects of management practices
#   on Northern Bobwhite Colinus virginianus density in privately owned working forests across 
#   the Southeastern United States). Point counts are incorporated via a hierarchical distance-sampling
#   model (Royle et al. 2004 Modeling abundance effects in distance sampling) with a Half-Normal
#   detection function (Yeiser et al. 2021. Addressing temporal variability in bird calling with design and estimation:
#   A Northern Bobwhite example), modified to account for repeat surveys. ARU data are incorporated via
#   a zero-truncated hurdle model and account for false positives and false negatives. We model abundance as varying by 
#   study site, while the ARU detection process varies by site and background noise on recordings.
#   A subset of detections are manually validated and incorporated with a Hypergeometric formulation
#   to correct for false positives.

# We show example code for the statistical model integrating point counts and ARU data, with a subset
#   of ARU detections manually verified (model AVC in Doser et al. 2021).


require(dclone)
require(rjags)



load("AR_NOBO_data.gzip")

data <- AR_NOBO_data



sink("NOBO_AR_AVC_model.jags")
cat("
    model{
    
    # b.0 are the site-specific intercepts for abundance modeling.
    # gamma.0 are the site-specific intercepts for modeling the parameter governing
    #   the call generation/call detection process on ARU recordings (delta).
    for(a in 1:nsites){
      b.0[a] ~ dnorm(0, 0.001)
      gamma.0[a] ~ dnorm(0, 0.001)
    }
    
    # sigma governs the scale of the Half-Normal detection function for point counts.
    beta ~ dunif(0, 10)
    log(sigma) <- beta
    
    # omega refers to the number of false positive detections on ARU recordings.
    omega ~ dunif(0.0001, 1000)
    
    # alpha.0 is the probability of detecting at least one call on ARU recordings at an
    #   unocuppied site (i.e., false positives).
    mu.alpha ~ dunif(0,1)
    alpha.0 <- logit(mu.alpha)
    
    # alpha.1 represents how the probability of detecting at least one call on an
    #   ARU recording increases with abundance. We constrained alpha.1 to be positive.
    #   as in Doser et al. (2021)
    alpha.1 ~ dunif(0,1000)

    
    # gamma.1 is the effect of background noise for modeling the parameter governing
    #   the call generation/call detection process on ARU recordings (delta).
    #   We constrained gamma.1 parameter to be negative.
    gamma.1 ~ dunif(-1000,0)

    
    # Average of site effects
    b.siteavg <- sum(b.0[1:nsites])/(nsites)
    
    # Total abundance
    N_ARU_Tot <- sum(N_ARU[1:npoints])
    N_PC_Tot <- sum(N_PC[1:npoints])
    
    # Calculating the probability of detecting a bird in each distance bin from 
    #   point counts. This is a function of the Half-Normal detection function
    #   and the relative area of each distance bin.
    for(b in 1:nB){
        pd_bins[b] <- (sigma^2*(1-exp(-db[b+1]^2/(2*sigma^2)))-sigma^2*(1-exp(-db[b]^2/(2*sigma^2))))*2*3.1416/(areaPC*pix[b])
        pd_bins_adj[b] <- pd_bins[b]*pix[b]
    }
    # p is the overall detection probability on point counts.
    p <- sum(pd_bins_adj[1:nB])
    
    
    # Manual validation of a subset of ARU detections. This is used to estimate
    #   the true positive rate at an ARU (tp), which in term estimates omega.
    for(v in 1:n.ARU.val){
      tp[v] <- (delta[ARUID.sub.val[v], ARU.Atimes.ID.val[v]] * N_ARU[ARU.pointID[ARUID.val[v]]]) / (delta[ARUID.sub.val[v], ARU.Atimes.ID.val[v]] * N_ARU[ARU.pointID[ARUID.val[v]]] + omega)
      K[v] ~ dbin(tp[v], v.ARU[ARUID.val[v], ARUday.val[v]])
      k.val[v] ~ dhyper(K[v], v.ARU[ARUID.val[v], ARUday.val[v]] - K[v], n.val[v], 1)
    }
    
    
    for(i in 1:npoints){


      ############################################################################
      ## State Model
      ############################################################################
      
      # ARUs and point counts have different detection radii, with point counts surverying
      #   the entire ARU area plus an additional area. We model abundance within the ARU survey
      #   area (N_ARU) and the expanded area surveyed by point counts but not ARUs (N_diff).
      #   We nclude survey area as an offset. This method ensures that abundance on point counts
      #   (N_PC) is greater than or equal to N_ARU.
    
      N_diff[i] ~ dpois(lambda_diff[i])
      lambda_diff[i] <- exp(b.0[site[i]] + log(areaDiff/10000))
      
      N_ARU[i] ~ dpois(lambda_ARU[i])
      lambda_ARU[i] <- exp(b.0[site[i]] + log(areaARU/10000))
      
      N_PC[i] <- N_ARU[i] + N_diff[i]
      
      
      
      
      ##########################################################################
      ## Observation Process
      ##########################################################################
      
      
      ##### Point Counts ######################################################
      
      for(j in 1:nvisits[i]){
      
          # Number detected
          
          y.PC[i,j] ~ dbin(p, N_PC[i])
          
          # Model fit of y.PC
          y.PC.pred[i,j] ~ dbin(p, N_PC[i])
          y.PC.resid[i,j] <- pow(pow(y.PC[i, j], 0.5) - pow(p * N_PC[i], 0.5), 2)
          y.PC.resid.pred[i,j] <- pow(pow(y.PC.pred[i, j], 0.5) - pow(p * N_PC[i], 0.5), 2)
          
          # Bin of detection
          
          ydb.PC[i,1:nB,j] ~ dmulti(pd_bins_adj[1:nB],y.PC[i,j])
          
          # Model fit of ydb.PC
          ydb.PC.pred[i,1:nB,j] ~ dmulti(pd_bins_adj[1:nB],y.PC[i,j])
          for(n in 1:nB){
            ydb.PC.resid.nB[i,n,j] <- pow(pow(ydb.PC[i,n,j], 0.5) - pow(pd_bins_adj[n] * y.PC[i,j], 0.5), 2)
            ydb.PC.resid.nB.pred[i,n,j] <- pow(pow(ydb.PC.pred[i,n,j], 0.5) - pow(pd_bins_adj[n] * y.PC[i,j], 0.5), 2)
          }
      }
      
      # Model fit for y.PC and ydb.PC
      tmp.y.PC[i] <- sum(y.PC.resid[i,1:nvisits[i]])
      tmp.y.PC.pred[i] <- sum(y.PC.resid.pred[i,1:nvisits[i]])
      tmp.ydb.PC[i] <- sum(ydb.PC.resid.nB[i,1:nB,1:nvisits[i]])
      tmp.ydb.PC.pred[i] <- sum(ydb.PC.resid.nB.pred[i,1:nB,1:nvisits[i]])
    }
    
    # Bayesian P-values for point counts
    fit_y.PC <- sum(tmp.y.PC[1:npoints])
    fit_y.PC_pred <- sum(tmp.y.PC.pred[1:npoints])
    fit_ydb.PC <- sum(tmp.ydb.PC[1:npoints])
    fit_ydb.PC_pred <- sum(tmp.ydb.PC.pred[1:npoints])
    bp_y.PC <- step(fit_y.PC_pred - fit_y.PC)
    bp_ydb.PC <- step(fit_ydb.PC_pred - fit_ydb.PC)
    
    
    
    ###### ARUs ###############################################################
    
    for(g in 1:n.ARU){
    
      # Probability of at least one detection on a recording
      logit(p.a[g]) <- alpha.0 + alpha.1*N_ARU[ARU.pointID[g]]
      
      for(r in 1:maxdates){
        # Estimating overdispersion parameter for each recording
        phi[g, r] ~ dgamma(1,1)
      }
      
      for(k in 1:n.ARU.dates[g]){
      
        # Presence/absence of detections on ARU recordings
        
        y.ARU[g, k] ~ dbin(p.a[g], 1)
        
        # Model fit of y.ARU
        y.ARU.pred[g,k] ~ dbin(p.a[g], 1)
        y.ARU.resid[g,k] <- pow(pow(y.ARU[g, k], 0.5) - pow(p.a[g], 0.5), 2)
        y.ARU.resid.pred[g,k] <- pow(pow(y.ARU.pred[g, k], 0.5) - pow(p.a[g], 0.5), 2)
      }
      
      # Model fit of y.ARU
      tmp.y.ARU[g] <- sum(y.ARU.resid[g,1:n.ARU.dates[g]])
      tmp.y.ARU.pred[g] <- sum(y.ARU.resid.pred[g,1:n.ARU.dates[g]])
    }  
    
    for(q in 1:n.which.v.ARU){
      for(p in 1:n.A.times[which.v.ARU[q]]){
      
        # Each recording gets its own delta, which governs the call generation/detection rate.
        # This is modeled via effects of study site and background noise
        log(delta[q, p]) <- gamma.0[site[ARU.pointID[which.v.ARU[q]]]] + gamma.1*ARU.noise[which.v.ARU[q], A.times[which.v.ARU[q], p]]
      
        
        # The number of detections, conditional on at least one being detected, are a function of 
        #     delta (calling rate/detectability), abundance, omega (false positives), 
        #     and presence/absence of detections
        # Adding in overdispersion parameter (phi) to aid convergence
        # Truncating to minimum of 1 to prevent observations with y.ARU=1 from having v.ARU=0
        v.ARU[which.v.ARU[q], A.times[which.v.ARU[q], p]] ~ dpois((delta[q, p] * N_ARU[ARU.pointID[which.v.ARU[q]]] + omega) * phi[which.v.ARU[q], A.times[which.v.ARU[q],p]] * y.ARU[which.v.ARU[q], A.times[which.v.ARU[q], p]])T(1, )
        
        # Model fit for v.ARU
        v.ARU.pred[q, p] ~ dpois((delta[q, p] * N_ARU[ARU.pointID[which.v.ARU[q]]] + omega) * phi[which.v.ARU[q], A.times[which.v.ARU[q],p]] * y.ARU[which.v.ARU[q], A.times[which.v.ARU[q], p]])T(1, )
        mu.v.ARU[q, p] <- ((delta[q, p] * N_ARU[ARU.pointID[which.v.ARU[q]]] + omega) * phi[which.v.ARU[q], A.times[which.v.ARU[q],p]]) / (1 - exp(-1 * ((delta[q, p] * N_ARU[ARU.pointID[which.v.ARU[q]]] + omega) * phi[which.v.ARU[q], A.times[which.v.ARU[q],p]])))
        v.ARU.resid[q,p] <- pow(pow(v.ARU[which.v.ARU[q], A.times[which.v.ARU[q], p]], 0.5) - pow(mu.v.ARU[q, p], 0.5), 2)
        v.ARU.resid.pred[q,p] <- pow(pow(v.ARU.pred[q, p], 0.5) - pow(mu.v.ARU[q, p], 0.5), 2)
      }
      
      # Model fit for v.ARU
      tmp.v.ARU[q] <- sum(v.ARU.resid[q,1:n.A.times[which.v.ARU[q]]])
      tmp.v.ARU.pred[q] <- sum(v.ARU.resid.pred[q,1:n.A.times[which.v.ARU[q]]])
    }

    
    # Bayesian P-values for ARUs
    fit_y.ARU <- sum(tmp.y.ARU[1:n.ARU])
    fit_y.ARU_pred <- sum(tmp.y.ARU.pred[1:n.ARU])
    fit_v.ARU <- sum(tmp.v.ARU[1:n.which.v.ARU])
    fit_v.ARU_pred <- sum(tmp.v.ARU.pred[1:n.which.v.ARU])
    bp_y.ARU <- step(fit_y.ARU_pred - fit_y.ARU)
    bp_v.ARU <- step(fit_v.ARU_pred - fit_v.ARU)
    }
    ",fill=TRUE)
sink()


# Initial values
K.inits <- rep(NA, times=data$n.ARU.val)
for(i in 1:length(K.inits)){
  K.inits[i] <- max(data$n.val[i] - 1, data$k.val[i])
}
inits.fun <- function() list(b.0=runif(AR_NOBO_data$nsites,-2,2),
                              beta=runif(1,4.5,5.5),
                              N_ARU=rep((max(AR_NOBO_data$y.PC,na.rm=T)+1),times=AR_NOBO_data$npoints),
                              N_diff=rep((max(AR_NOBO_data$y.PC,na.rm=T)+1),times=AR_NOBO_data$npoints),
                              mu.alpha=runif(1,0.1,0.9),
                              alpha.1=runif(1,0,1),
                              gamma.0=runif(AR_NOBO_data$nsites,-1,1),
                              gamma.1=runif(1,-1,0),
                              omega=runif(1,0.1,3),
                              K=K.inits,
                              phi=matrix(runif(max(AR_NOBO_data$n.ARU.dates)*AR_NOBO_data$n.ARU,0.9,1.1),nrow=AR_NOBO_data$n.ARU))


params <- c("b.0","N_ARU","N_PC","sigma","p","b.sitesavg","delta","gamma.0",
             "gamma.1","omega","p.a","alpha.0","alpha.1","phi","bp_y.PC",
             "bp_ydb.PC", "bp_y.ARU", "bp_v.ARU", "N_ARU_Tot", "N_PC_Tot")



cl <- makePSOCKcluster(3)
NOBO.model.AR.AVC <- jags.parfit(cl, data=data, params=params,
                                        model='NOBO_AR_AVC_model.jags', inits=inits.fun(),
                                        n.chains=3, n.iter=200000, n.burnin=60000,n.adapt=60000,
                                        thin=50)

stopCluster(cl)
