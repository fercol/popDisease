# ============================== CODE METADATA =============================== #
# AUTHOR: Fernando Colchero <colchero@imada.sdu.dk>
# DATE CREATED: 2021-03-24
# DESCRIPTION: 
# COMMENTS: 
# ================================ CODE START ================================ #
# =============================== #
# ==== POPULATION PROJECTION: ====
# =============================== #
# Function to run SIRS projection model:
ProjPopSIRS <- function(tFin = 10, demog, infIni = 1) {
  # Time sequence:
  tseq <- seq(0, tFin, demog$dx)
  nt <- length(tseq)
  
  # States names:
  stateName <- c("sus", "inf", "imm", "dea")
  stateFullName <- c(sus = "Susceptible", inf = "Infected", imm = "Immune", 
                     dea = "Dead")
  nStates <- length(stateName)
  
  # Matrices for each state (suscept., infected, recovered, dead, immune):
  states <- list()
  for (st in 1:nStates) {
    states[[stateName[st]]] <- array(0, dim = c(demog$nx, nt, 2),
                                     dimnames = list(NULL, NULL, c("F", 'M')))
  }
  
  # Immunity duration:
  immDur <- demog$disease$immDur
  
  # Additional immune array to account for immune response decline in time:
  immState <- array(0, dim = c(demog$nx, immDur, 2),
                    dimnames = list(NULL, NULL, c("F", 'M')))
  
  # # Fill up initial states:
  states$sus[, 1, ] <- demog$susStart
  
  # Assign initial infected individual:
  infStart <- AssignInfected(tt = 1, nInf = infIni, demog = demog,
                             states = states)
  states$sus[, 1, ] <- states$sus[, 1, ] - infStart
  states$inf[, 1, ] <- infStart
  
  # Index of states without death:
  idNoDI <- which(stateName != "dea" & stateName != "imm")
  
  # ------------------------- #
  # Run projection model ----
  # ------------------------- #
  for (tt in 2:nt) {
    # Baseline survival:
    for (sx in c("F", "M")) {
      # Susceptibles:
      sust0 <- states$sus[-demog$nx, tt - 1, sx]
      sust1 <- rbinom(n = demog$nx - 1, size = sust0, 
                      prob = demog$px[[sx]][-demog$nx])
      deat1 <- sust0 - sust1
      
      # Infected:
      infsurv <- 1 - demog$disease$qx
      inft0 <- states$inf[-demog$nx, tt - 1, sx]
      inft1 <- rbinom(n = demog$nx - 1, size = inft0, 
                      prob = demog$px[[sx]][-demog$nx] * infsurv[-1])
      deat1 <- deat1 + inft0 - inft1
      
      # Immune:
      for (tim in 1:immDur) {
        immSurvt <- rbinom(n = demog$nx - 1, 
                           size = immState[-demog$nx, tim, sx], 
                           prob = demog$px[[sx]][-demog$nx])
        
        # Assign deaths:
        deat1 <- deat1 + immState[-demog$nx, tim, sx] - immSurvt
        
        # Assign to immune table or to succeptibles:
        if (tim == 1) {
          sust1 <- sust1 + immSurvt
        } else {
          immState[-1, tim - 1, sx] <- immSurvt
        }
      }
      
      # New immune individuals:
      immune <- rbinom(n = demog$nx - 1, size = inft1, 
                       prob = demog$disease$immProb[-1])
      immState[-1, immDur, sx] <- immune
      inft1 <- inft1 - immune
      sust1 <- sust1 + inft1
      
      # Assign to states:
      states$sus[-1, tt, sx] <- sust1
      states$dea[-1, tt, sx] <- deat1
      if (immDur > 1) {
        states$imm[, tt, sx] <- c(immState[,, sx] %*% rep(1, immDur))
      } else {
        states$imm[, tt, sx] <- immState[,, sx]
      }
    }
    
    # New infections:
    nSucept <- sum(states$sus[, tt, ])
    nInf0 <- sum(states$inf[, tt - 1, ])
    nInft <- min(rpois(1, nInf0 * demog$disease$R0), nSucept)
    infTabt <- AssignInfected(tt = tt, nInf = nInft, demog = demog,
                              states = states)
    states$sus[, tt, ] <- states$sus[, tt, ] - infTabt
    states$inf[, tt, ] <- infTabt
    
    # Number of new-born babies:
    nMales <- sum(states$sus[demog$idAdult, tt, "M"] + 
                    states$imm[demog$idAdult, tt, "M"])
    if (nMales > 1) {
      nFems <- states$sus[, tt, "F"] + states$imm[, tt, "F"]
      nBabies <- sum(rbinom(demog$nx,size = nFems, prob = demog$fx))
      nBabiesF <- rbinom(n = 1, size = nBabies, prob = demog$sexRatio["F"])
      nBabiesM <- nBabies - nBabiesF
      states$sus[1, tt, "F"] <- nBabiesF
      states$sus[1, tt, "M"] <- nBabiesM
    } else {
      nBabiesF <- 0
      nBabiesM <- 0
    }
  }
  
  # Extract population sizes:
  Mt <- rep(0, nt)
  for (st in 1:nStates) {
    Mt <- Mt + apply(states[[st]], 2, sum)
  }
  
  # Output:
  outSIRS <- list(states = states, tseq = tseq, nt = nt, nStates = nStates,
                  stateName = stateName, stateFullName = stateFullName, 
                  Mt = Mt)
  
  # Attribute:
  attr(outSIRS, "type") <- "singleSIRS"
  
  return(outSIRS)  
}

# =============================== #
# ==== DEMOGRAPHIC FUNCTIONS: ====
# =============================== #
# Age-specific fecundity:
fecfun <- function(x, b) exp(b[1] + b[2] * x - b[3] * x^2 + b[4] / (x + 1))

# Age-specific hazards rate (mortality function):
mortfun <- function(x, b) exp(b[1] - b[2] * x) + b[3] + exp(b[4] + b[5] * x)

# Cumulative survival:
survfun <- function(x, b) {
  exp(exp(b[1]) / b[2] * (exp(-b[2] * x) - 1) - b[3] * x + exp(b[4]) / b[5] * 
        (1 - exp(b[5] * x)))
}

# ============================ #
# ==== RESULTS MANAGEMENT: ====
# ============================ #
# Function to run multiple SIRS projection models in parallel:
ParalProjSIRS <- function(nsim, ncpus, tFin = 10, demog, infIni = 1) {
  
  # internal parallel function:
  IntParalFun <- function(sim, ...) ProjPopSIRS(...)
  
  # List of objects to be uploaded to the cores:
  cpuvars <- c("AssignInfected", "ProjPopSIRS", "demog")
  
  # Start computation time:
  Start <- Sys.time()
  
  # Run multiple simulations in parallel:
  sfInit(parallel = TRUE, cpus = ncpus)
  
  # Export variables to cpus:
  sfExport(list = cpuvars)
  
  # Run common prep functions on all cpus:
  # sfSource("prepfunct.R")
  
  # Run parallel executions:
  outParalSIRS <- sfClusterApplyLB(x = 1:nsim, fun = IntParalFun, tFin = tFin,
                                   demog = demog, infIni = infIni)
  
  # Stop paralell executions:
  sfStop()
  
  # End computation time:
  End <- Sys.time()
  
  # Output list:
  outParallelSIRS <- list(out = outParalSIRS, nsim = nsim, tFin = tFin, 
                          demog = demog, infIni = infIni, 
                          compTime = End - Start)
  
  # Assign attribute:
  attr(outParallelSIRS, "type") <- "parallelSIRS"
  
  # Return output:
  return(outParallelSIRS)
}

# Function to extract outputs from Parallel projection SIRS:
ExtractParalSIRS <- function(outParal) {
  # Number of simulations:
  nsim <- outParal$nsim
  
  # Number of time steps:
  tseq <- outParal$out[[1]]$tseq
  ntseq <- length(tseq)
  
  # Number of states:
  nStates <- outParal$out[[1]]$nStates
  
  # Calculate mean state values per time interval, 
  # total population sizes and number of extinct populations:
  statesMean <- outParal$out[[1]]$states
  nExt <- rep(0, ntseq)
  Mtmat <- matrix(0, nsim, ntseq)
  for (ns in 2:nsim) {
    statens <- outParal$out[[ns]]$states[[1]] * 0
    for (ii in 1:nStates) {
      statesMean[[ii]] <- statesMean[[ii]] +
        outParal$out[[ns]]$states[[ii]]
      statens <- statens + outParal$out[[ns]]$states[[ii]]
    }
    sumstates <- outParal$out[[ns]]$Mt
    id0 <- which(sumstates == 0)[1]
    if (length(id0) > 0) {
      nExt[id0] <- nExt[id0] + 1
    }
    Mtmat[ns, ] <- sumstates
  }
  
  # Calculate mean and 95% quantiles of population size:
  Mtq <- rbind(apply(Mtmat, 2, mean), 
               apply(Mtmat, 2, quantile, c(0.025, 0.33, 0.5, 0.67, 0.975)))
  
  # Calculate sd of state values per time interval:
  statesVar <- outParal$out[[1]]$states
  for (ns in 1:nsim) {
    for (ii in 1:nStates) {
      if (ns == 1) {
        statesVar[[ii]] <- statesVar[[ii]] * 0
      }
      idNo0 <- which(outParal$out[[ns]]$states[[ii]] != 0)
      statesVar[[ii]][idNo0] <- (statesVar[[ii]] +
                                   (outParal$out[[ns]]$states[[ii]] - 
                                      statesMean[[ii]])^2)[idNo0]
    }
  }
  
  # Calculate mean and variance:
  for (ii in 1:nStates) {
    statesMean[[ii]] <- statesMean[[ii]] / nsim
    statesVar[[ii]] <- statesVar[[ii]] / (nsim - 1)
  }
  
  # Summary output list:
  summOut <- list(statesMean = statesMean, statesVar = statesVar, 
                  settings = c(nsim = nsim, tFin = outParal$tFin, 
                               nt = outParal$out[[1]]$nt,
                               infIni = outParal$infIni),
                  demog = outParal$demog, nStates = nStates, 
                  stateName = outParal$out[[1]]$stateFullName,
                  tseq = outParal$out[[1]]$tseq, Mt = Mtq, Mtmat = Mtmat, 
                  nExt = nExt)
  
  # Assign attribute:
  attr(summOut, "type") <- "summSIRS"
  
  return(summOut)
}

# Plotting function:
plotSIRS <- function(x, type = "l", sim = 1, ...) {
  # Object to reset plotting parameters:
  op <- par(no.readonly = TRUE)
  
  # dotted arguments:
  argsl <- list(...)
  
  # Extract data based on attribute:
  if (attr(x, "type") == "singleSIRS") {
    projStates <- x$states
    nStates <- x$nStates
    tseq <- x$tseq
    nt <- x$nt
    stateFullName <- x$stateFullName
  } else if (attr(x, "type") == "parallelSIRS") {
    projStates <- x$out[[sim]]$states
    nStates <- x$out[[sim]]$nStates
    tseq <- x$out[[sim]]$tseq
    nt <- x$out[[sim]]$nt
    stateFullName <- x$out[[sim]]$stateFullName
  } else if (attr(x, "type") == "summSIRS") {
    projStates <- x$statesMean
    nStates <- x$nStates
    tseq <- x$tseq
    nt <- x$settings["nt"]
    stateFullName <- x$stateName
  }
  
  # Plot results:
  if (type == "l") {
    par(mfrow = c(4, 2), mar = c(4, 4, 1, 1))
    for (ii in 1:nStates) {
      for (sx in c("F", "M")) {
        plot(tseq, apply(projStates[[ii]][,, sx], 2, sum), type = 'l', 
             lwd = 2, xlab = "Time (years)", ylab = "Number of individuals")
        mtext(sprintf("%s - %s", stateFullName[ii], sx), side = 3, line = -2)
      }
    }
    
  } else {
    # Extract totals:
    tots <- lapply(c("F", "M"), function(sx) {
      t(cbind(sus = apply(projStates$sus[,, sx], 2, sum), 
              inf = apply(projStates$inf[,, sx], 2, sum), 
              imm = apply(projStates$imm[,, sx], 2, sum), 
              dea = apply(projStates$dea[,, sx], 2, sum)))
    })
    names(tots) <- c("F", "M")
    
    # ylim argument:
    if ("ylim" %in% names(argsl)) {
      ylim <- argsl$ylim
    } else {
      maxy <- max(sapply(c("F", "M"), function(sx) {
        max(apply(tots[[sx]], 2, sum))
      }))
      ylim <- c(0, maxy)
      ylim[2] <- ceiling(maxy / 10) * 10
    }
    
    if ("col" %in% names(argsl)) {
      col <- argsl$col
    } else {
      col <- brewer.pal(nrow(tots$F), "Set1")
      names(col) <- rownames(tots$F)
    }
    labXaxis <- 0:max(tseq)
    
    par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
    for (sx in c("F", "M")) {
      if (sx == "M") {
        legText <- stateFullName
      } else {
        legText <- FALSE
      }
      bp <- barplot(tots[[sx]], col = col, border = col, 
                    legend.text = legText, ylim = ylim, axes = FALSE)
      Mbp <- max(bp)
      mbp <- min(bp)
      a <- (Mbp - mbp) / diff(range(tseq))
      b <- Mbp - a * max(tseq)
      mtext(c(F = "Females", M = "Males")[sx], side = 3, line = 0)
      atx <- labXaxis * a + b
      axis(at = atx, labels = labXaxis, side = 1)
      Axis(ylim, side = 2, las = 2)
      mtext("Number of Inds.", side = 2, line = 3)
    }
    mtext("Time (years)", side = 1, line = 2)
  }
  par(op)
}

# ============================ #
# ==== INTERNAL FUNCTIONS: ====
# ============================ #
# Function to find available ages:
AssignInfected <- function(tt, nInf, demog, states) {
  infTab <- states$inf[, tt, ] * 0
  rownames(infTab) <- as.character(demog$x)
  if (nInf > 0) {
    ageVec <- c()
    sexVec <- c()
    totSus <- sum(states$sus[-1, tt, ])
    totImm <- sum(states$imm[-1, tt, ])
    totAvail <- totSus + totImm
    nAvailTot <- rbinom(1, size = nInf, prob = totSus / totAvail)
    
    for (sx in c("F", "M")) {
      idx <- which(states$sus[, tt, sx] > 0 & demog$x > 0)
      xAvail <- demog$x[idx]
      nAvail <- states$sus[idx, tt, sx]
      ageVec <- c(ageVec, rep(xAvail, nAvail))
      sexVec <- c(sexVec, rep(sx, sum(nAvail)))
    }
    idInf <- sample(1:length(ageVec), size = nAvailTot, 
                    replace = FALSE)
    
    tabInf <- table(ageVec[idInf], sexVec[idInf])
    ageInf <- rownames(tabInf)
    
    # Assign new infections:
    for (sx in unique(sexVec[idInf])) {
      infTab[ageInf, sx] <- tabInf[, sx]
    }
  }
  return(infTab)
}

# ================================= CODE END ================================= #
