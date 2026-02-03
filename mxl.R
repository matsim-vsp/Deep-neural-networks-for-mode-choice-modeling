# ################################################################# #
#### LOAD LIBRARY AND DEFINE CORE SETTINGS                       ####
# ################################################################# #

### Clear memory
rm(list = ls())

### Load Apollo library
library(apollo)
library(tidyverse)
library(dplyr)
library(patchwork)
library(networkD3)
library(sf) #=> geography
library(parallel)
library(stringr)
library(ggalluvial)
library(reshape2)

### Set working directory (only works in RStudio)
apollo_setWorkDir()

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelDescr      = "Very simple mixed logit model based on CR mode choice MATSim data",
  indivID         = "ID", 
  outputDirectory = "output",
  modelName = "simpleMXL",
  weights = "weight",
  nCores = detectCores() - 1
)

# ################################################################# #
#### LOAD DATA AND APPLY ANY TRANSFORMATIONS                     ####
# ################################################################# #

### Loading data from package
### if data is to be loaded from a file (e.g. called data.csv), 
### the code would be: database = read.csv("data.csv",header=TRUE)

database = read_csv("/Users/gregorr/Documents/work/respos/shared-svn/projects/matsim-berlin/data/SrV/2018/converted/trip-choices.csv", comment = "#")  %>%
  rename(ID = person) %>%
  rename(util_money_inp = util_money)


### for data dictionary, use ?apollo_modeChoiceData

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta=c(b_tt_car  = 0,
              b_tt_pt  = 0,
              b_tt_walk = 0,
              b_tt_bike = 0,
              b_tt_ride = 0,
              mu_asc_car = 0,
              mu_asc_pt = 0,
              mu_asc_bike = 0,
              mu_asc_ride = 0,
              sigma_asc_car = 0,
              sigma_asc_pt = 0,
              sigma_asc_bike = 0,
              sigma_asc_ride = 0,
              asc_walk = 0
            )

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c("asc_walk")



# ################################################################# #
#### DEFINE RANDOM COMPONENTS                                    ####
# ################################################################# #

### Set parameters for generating draws
apollo_draws = list(
  interDrawsType = "halton",
  interNDraws    = 800,
  interUnifDraws = c(),
  interNormDraws = c("draws_asc_car","draws_asc_pt","draws_asc_bike","draws_asc_ride"),
  intraDrawsType = "halton",
  intraNDraws    = 0,
  intraUnifDraws = c(),
  intraNormDraws = c()
)

### Create random parameters
apollo_randCoeff = function(apollo_beta, apollo_inputs){
  randcoeff = list()
  #randcoeff[["asc_car"]] = -exp( mu_log_b_tt + sigma_log_b_tt * draws_tt )
  randcoeff[["asc_car"]] = mu_asc_car + draws_asc_car * sigma_asc_car 
  randcoeff[["asc_bike"]] = mu_asc_bike + draws_asc_bike * sigma_asc_bike
  randcoeff[["asc_pt"]] = mu_asc_pt + draws_asc_pt * sigma_asc_pt
  randcoeff[["asc_ride"]] = mu_asc_ride + draws_asc_ride * sigma_asc_ride
  return(randcoeff)
}


# ################################################################# #
#### GROUP AND VALIDATE INPUTS                                   ####
# ################################################################# #

apollo_inputs = apollo_validateInputs()

# ################################################################# #
#### DEFINE MODEL AND LIKELIHOOD FUNCTION                        ####
# ################################################################# #

apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  ### List of utilities: these must use the same names as in mnl_settings, order is irrelevant
  V = list()
  V[["car"]]  = asc_car  + b_tt_car  * car_hours
  V[["pt"]]  = asc_pt  + b_tt_pt  * pt_hours 
  V[["bike"]]  = asc_bike  + b_tt_bike * bike_hours 
  V[["walk"]] = asc_walk + b_tt_walk * walk_hours
  V[["ride"]] = asc_ride + b_tt_ride * ride_hours
  
  ### Define settings for MNL model component
  mnl_settings = list(
    alternatives  = c(walk=1, pt=2, car=3, bike=4, ride=5), 
    avail         = list(walk = walk_valid, pt= pt_valid, car = car_valid, bike = bike_valid, ride = ride_valid), 
    choiceVar     = choice,
    utilities     = V
  )
  
  ### Compute probabilities using MNL model
  P[["model"]] = apollo_mnl(mnl_settings, functionality)
  
  ### Take product across observation for same individual
  P = apollo_panelProd(P, apollo_inputs, functionality)
  
  ## Use weighting 
  if (functionality == "estimate") {
      P= apollo_weighting(P,
                          apollo_inputs,
                          functionality)
  }
  
  ### Average across inter-individual draws
  P = apollo_avgInterDraws(P, apollo_inputs, functionality)
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# ################################################################# #
#### MODEL ESTIMATION                                            ####
# ################################################################# #


estimate_settings = list(
  maxIterations = 400,           # Set maximum iterations
  maxFunctionEvaluations = 86400, # Set maximum function evaluations
  algorithm = "BFGS"              # Switch to BFGS for more stable convergence
)

model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, estimate_settings)

apollo_prediction(model, apollo_probabilities, apollo_inputs )

P_predicted = apollo_probabilities(model$estimate, apollo_inputs, functionality = "prediction")

# Write predictions

df <- bind_rows(P_predicted) %>%
  add_column(ID = database$ID, trip_n = database$trip_n + 1, choice = database$choice, .before = 0)

write_csv(df, "../mxl_output.csv")

# Step 1: Extract the predicted probabilities for each alternative
prob_walk = P_predicted$model$walk
prob_pt = P_predicted$model$pt
prob_car = P_predicted$model$car
prob_bike = P_predicted$model$bike
prob_ride = P_predicted$model$ride

# Step 2: Compute the mean probability (market share) for each alternative
market_share_walk = mean(prob_walk)
market_share_pt = mean(prob_pt)
market_share_car = mean(prob_car)
market_share_bike = mean(prob_bike)
market_share_ride = mean(prob_ride)

# Step 3: Display the predicted market shares
predicted_shares = c(
  Walk = market_share_walk,
  Public_Transport = market_share_pt,
  Car = market_share_car,
  Bike = market_share_bike,
  Ride = market_share_ride
)

# Print the predicted market shares
print(predicted_shares)

# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO SCREEN)                               ----
# ----------------------------------------------------------------- #

apollo_modelOutput(model)

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO FILE, using model name)               ----
# ----------------------------------------------------------------- #

apollo_saveOutput(model)
