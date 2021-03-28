# ============================== CODE METADATA =============================== #
# AUTHOR: Fernando Colchero <colchero@imada.sdu.dk>
# DATE CREATED: 2021-03-24
# DESCRIPTION: Code to run the population dynamics-epidemiological model 
#              proposed by Colchero, Eckardt, Stoinski (2021).
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

# Sourced files:
source("02code/popDiseaseFunctions.R")

# ================================= #
# ==== DEMOGRAPHIC DATA OBJECT: ====
# ================================= #
# ---------------------------- #
# Demographic parameters: ----
# ---------------------------- #
# Mortality parameters:
theta <- cbind(F = c(a0 = -0.229, a1 = 1.882, c = 0.013, 
                     b0 = -10.275, b1 = 0.223), 
               M = c(a0 = 0.004, a1 = 2.382, c = 0.021, 
                     b0 = -7.52, b1 = 0.189))

# Reproduction parameters:
beta <- c(r0 = -1.216, r1 = 0.022, r2 = 0.001, r3 = -0.972)

# ---------------------------------------- #
# Create baseline demographic object: ----
# ---------------------------------------- #
# Demographic object:
demog <- CreateDemoObj(theta = theta, beta = beta, N = c(F = 60, M = 50),
                       omega = 60, alpha = 8, propM = 0.523, yearIni = 2021,
                       R0 = 2, qMax = 0.3, maxImPr = 0.8, immDur = 3)

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
    (1 + exp(-0.1 * (demog$x - demog$settings["omega"] * 0.4)))
  
  # Replace immunity duration:
  demogn$disease$immDur <- covidVars$immDur[vv]
  
  # Mortality probability:
  demogn$disease$qx <- covidVars$maxQx[vv] / 
    (1 + exp(-0.2 * (demog$x - demog$settings["omega"] * 0.4)))
  
  # Run multiple models in parallel:
  outParal <- ParalProjSIRS(nsim = 2000, ncpus = 4, tFin = 10, 
                            demog = demogn, infIni = infIni)
  
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
