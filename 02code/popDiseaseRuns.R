# ============================== CODE METADATA =============================== #
# AUTHOR: Fernando Colchero <colchero@imada.sdu.dk>
# DATE CREATED: 2021-03-24
# DESCRIPTION: 
# COMMENTS: 
# ================================ CODE START ================================ #
# ================ #
# ==== SETUP: ==== 
# ================ #
# Libraries:
library(snowfall)
library(colorRamps)
library(RColorBrewer)

# Working directory:
# setwd("Path to main directory")

# Logical for saving results:
saveResults <- FALSE

# Sourced files:
source("02code/popDiseaseFunctions.R")

# ================================= #
# ==== DEMOGRAPHIC DATA OBJECT: ====
# ================================= #
# ---------------------------- #
# Demographic parameters: ----
# ---------------------------- #
# Create list of demographic parameters, values and rates:
demog <- list()

# Read-in parameter values:
demog$pars <- list()

# Mortality parameters (as from study population):
demog$pars$mort <- cbind(F = c(a0 = -0.229, a1 = 1.882, c = 0.013, 
                               b0 = -10.275, b1 = 0.223), 
                         M = c(a0 = 0.004, a1 = 2.382, c = 0.021, 
                               b0 = -7.52, b1 = 0.189))

# Fecundity:
demog$pars$fec <- c(r0 = -1.216, r1 = 0.022, r2 = 0.001, r3 = -0.972)

# -------------------------------- #
# Calculate monthly survival: ----
# -------------------------------- #
# Total number of individuals in pop:
N <- c(F = 60, M = 50)

# Increments in months:
dx <- 1 / 12

# Maximum longevity:
omega <- 60

# Daily age vector in year units: 
demog$x <- seq(0, 60, dx)
demog$nx <- length(demog$x)
demog$dx <- dx

# Index of adult ages:
demog$idAdult <- which(demog$x >= 8)
demog$nxf <- length(demog$idAdult)

# Calculate monthly age-specific survival per sex:
demog$px <- list()
for (sx in c("F", "M")) {
  theta <- demog$pars$mort[, sx]
  demog$px[[sx]] <- survfun(demog$x + demog$dx, theta) / 
    survfun(demog$x, theta)
}

# Calculate monthly reproduction output
demog$fx <- fecfun(demog$x, demog$pars$fec) * demog$dx
demog$fx[demog$x < 8] <- 0

# Sex ratios:
propM <- 0.523
propF <- 1 - propM
demog$sexRatio <- c(F = propF, M = propM)

# initial year:
yearIni <- 2021

# Fill up initial susceptible states:
demog$susStart <- matrix(0, demog$nx, 2, dimnames = list(NULL, c("F", "M")))
xAgeSt <- 0:omega
for (sx in c("F", "M")) {
  theta <- demog$pars$mort[, sx]
  Sx <- survfun(xAgeSt, theta)
  nAgeVec <- rbinom(n = length(xAgeSt), size = N[sx], prob = Sx / sum(Sx))
  for (xx in 1:length(xAgeSt)) {
    idx <- which(abs(demog$x - xAgeSt[xx]) == 
                   min(abs(demog$x - xAgeSt[xx])))[1]
    nAge <- rbinom(n = 1, size = N[sx], prob = Sx[xx] / sum(Sx))
    demog$susStart[idx, sx] <- nAgeVec[xx]
  }
}

# ------------------------------ #
# Setup disease parameters: ----
# ------------------------------ #
demog$disease <- list()

# Contagion rate (R0):
demog$disease$R0 <- 2

# incubation time:
demog$disease$incTime <- 14

# Mortality probability:
demog$disease$qx <- 0.3 / (1 + exp(-0.2 * (demog$x - omega * 0.4)))

# Immunity rate:
demog$disease$immProb <- 0.8 / (1 + exp(-0.1 * (demog$x - omega * 0.4)))

#Immunity duration (in months):
demog$disease$immDur <- 10

# ========================== #
# ==== SINGLE MODEL RUN: ====
# ========================== #
# ---------------------------- #
# Setup projection model: ----
# ---------------------------- #
# Start and end times and sequence of time steps:
t0 <- 0
tFin <- 10

# Initial number of infected individuals:
infIni <- 1

# Run function:
disProj <- ProjPopSIRS(tFin = tFin, demog = demog, infIni = infIni)

# Plot results:
plotSIRS(disProj)
plotSIRS(disProj, type = "ab")

# ============================= #
# ==== MULTIPLE MODEL RUNS: ====
# ============================= #
# Vary R0:
R0vec <- c(0.5, 1, 2, 3)

# Vary maximum immunity:
maxImProb <- c(0.2, 0.4, 0.6, 0.8)

# Vary immunity duration:
immDurVec <- c(1, 3, 6, 12)

# Vary maximum mortality:
maxQxVec <- c(0.3, 0.6)

# Grid of new values:
covidVars <- expand.grid(R0 = R0vec, MImProb = maxImProb, immDur = immDurVec,
                         maxQx = maxQxVec)

# Output:
outMulti <- list()

for (vv in 1:nrow(covidVars)) {
  # Create new demog object:
  demogn <- demog
  
  # Replace reproductive value:
  demogn$disease$R0 <- covidVars$R0[vv]
  
  # Replace immunity probability
  demogn$disease$immProb <- covidVars$MImProb[vv] / 
    (1 + exp(-0.1 * (demog$x - omega * 0.4)))
  
  # Replace immunity duration:
  demogn$disease$immDur <- covidVars$immDur[vv]
  
  # Mortality probability:
  demogn$disease$qx <- covidVars$maxQx[vv] / 
    (1 + exp(-0.2 * (demog$x - omega * 0.4)))
  
  # Run multiple models in parallel:
  outParal <- ParalProjSIRS(nsim = 2000, ncpus = 4, tFin = 10, demog = demogn, 
                            infIni = infIni)
  
  # Extract average values:
  meanOutParal <- ExtractParalSIRS(outParal = outParal)
  
  # Plot barplots:
  cat(sprintf("%s - R0 = %s; ImmPr = %s; ImmDur = %s\n", vv, covidVars$R0[vv], 
              covidVars$MImProb[vv], covidVars$immDur[vv]))
  plotSIRS(meanOutParal, type = "bp")
  
  # Store results:
  outMulti[[vv]] <- meanOutParal
}


# ================================= CODE END ================================= #
