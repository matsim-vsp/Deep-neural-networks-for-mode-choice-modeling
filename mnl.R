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
library(stringr)
library(ggalluvial)
library(reshape2)

### Set working directory (only works in RStudio)
apollo_setWorkDir()

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelDescr      = "Simple MNL model based on CR mode choice MATSim data",
  indivID         = "ID",
  outputDirectory = "../../shared-svn/projects/matsim-berlin/data/SrV/2018/halil/output/",
  modelName = "CRspecifcation",
  weights = "weight"
)

# ################################################################# #
#### LOAD DATA AND APPLY ANY TRANSFORMATIONS                     ####
# ################################################################# #

### Loading data from package
### if data is to be loaded from a file (e.g. called data.csv),
### the code would be: database = read.csv("data.csv",header=TRUE)


##we need to split the data set for the dnn model
#testData <- read_csv("../../shared-svn/projects/matsim-berlin/data/SrV/2018/halil/Train-Test Data/train_data.csv",  comment = "#")
fullData <- read_csv(
  "../../shared-svn/projects/matsim-berlin/data/SrV/2018/halil/Train-Test Data/train_dataAllAttributes.csv",
  comment = "#"
)

#testData <- testData %>%
#  left_join(
#    fullData %>% select(person, trip_n, util_money),
#    by = c("person", "trip_n")
#  )


database <- fullData  %>%
  rename(ID = person) %>%
  rename(util_money_inp = util_money)


distances <- database %>% group_by(ID) %>%
  summarise(total_dist = sum(beelineDist))

# Costs calculation based on length of trip (related to total daily length)

database <- left_join(database, distances, by = "ID") %>%
  mutate(
    # Distance weight (same as Biogeme)
    distanceWeight = beelineDist / total_dist,
    distanceWeight = replace_na(distanceWeight, 1),
    # ---- KM COSTS ----
    car_dist_costs  = car_km  * -0.149,
    ride_dist_costs = ride_km * -0.149,
    # ---- FIXED COSTS (ONLY IF MODE IS USED!) ----
    car_fixed_costs = (car_valid == 1) * distanceWeight * -14.30,
    pt_fixed_costs  = (pt_valid  == 1) * distanceWeight * -3
  ) %>%
  # Disable the panel effect (keep as is for MNL)
  mutate(ID = row_number())

##check if there are any na in the database
sum(is.na(database))

### for data dictionary, use ?apollo_modeChoiceData

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta = c(
  asc_car   = 0,
  asc_bike   = 0,
  asc_pt   = 0,
  asc_walk  = 0,
  asc_ride = 0,
  b_tt = 6.88,
  #              b_tt_car  = 0,
  #              b_tt_pt  = 0,
  #              b_tt_walk = 0,
  #              b_tt_bike = 0,
  #              b_tt_ride = 0,
  util_money = 1,
  fixed_costs_perception = 0.3,
  exp_income = 0
)

#b_line_switches <- 1
### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c("asc_walk")

# ################################################################# #
#### GROUP AND VALIDATE INPUTS                                   ####
# ################################################################# #

apollo_inputs = apollo_validateInputs()

# Read mean income from csv, this value is put directly in the code
# mean_income <- read.table(file = "../trip-choices.csv", header = F,nrows = 1, comment.char = "", sep = ":")

# ################################################################# #
#### DEFINE MODEL AND LIKELIHOOD FUNCTION                        ####
# ################################################################# #

apollo_probabilities = function(apollo_beta,
                                apollo_inputs,
                                functionality = "estimate") {
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  # price * util_money * (income/mean(income))^b_elasticitypriceincome
  mean_income <- 1942
  
  ### List of utilities: these must use the same names as in mnl_settings, order is irrelevant
  V = list()
  V[["car"]]  = asc_car - b_tt  * car_hours + (car_fixed_costs * fixed_costs_perception + car_dist_costs) * util_money * (mean_income / income) ** exp_income
  V[["pt"]]  = asc_pt - b_tt * pt_hours + pt_fixed_costs  * fixed_costs_perception * util_money * (mean_income / income) ** exp_income - pt_switches
  V[["bike"]]  = asc_bike  - b_tt * bike_hours
  V[["walk"]] = asc_walk - b_tt * walk_hours
  # Ride has double the time costs because of driver + passenger
  V[["ride"]] = asc_ride - b_tt * ride_hours * 2 + ride_dist_costs * util_money * (mean_income / income) ** exp_income
  
  ### Define settings for MNL model component
  mnl_settings = list(
    alternatives  = c(
      walk = 1,
      pt = 2,
      car = 3,
      bike = 4,
      ride = 5
    ),
    avail         = list(
      walk = walk_valid,
      pt = pt_valid,
      car = car_valid,
      bike = bike_valid,
      ride = ride_valid
    ),
    choiceVar     = choice,
    utilities     = V
  )
  
  ### Compute probabilities using MNL model
  P[["model"]] = apollo_mnl(mnl_settings, functionality)
  
  ### Take product across observation for same individual
  # Currently not used
  #P = apollo_panelProd(P, apollo_inputs, functionality)
  
  ## Use weighting
  # For now only when estimating
  #if (functionality == "estimate") {
  P = apollo_weighting(P, apollo_inputs, functionality)
  #}
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# ################################################################# #
#### MODEL ESTIMATION                                            ####
# ################################################################# #

estimate_settings <- list(
  estimationRoutine = "BFGS",
  constraints = c(
    'util_money > 0.3 - 1e-10',
    'util_money < 1.5 + 1e-10',
    'fixed_costs_perception > 0 - 1e-10',
    'fixed_costs_perception < 1 + 1e-10',
    'exp_income > 0 - 1e-10',
    'exp_income < 1.5  + 1e-10',
    'b_tt > 0 + 1e-10',
    'b_tt < 15 + 1e-10'
  )
)

model = apollo_estimate(apollo_beta,
                        apollo_fixed,
                        apollo_probabilities,
                        apollo_inputs,
                        estimate_settings)

##extracting estimated parameters
beta_hat <- model$estimate

apollo_prediction(model, apollo_probabilities, apollo_inputs)

P_predicted = apollo_probabilities(model$estimate, apollo_inputs, functionality = "prediction")

# Write predictions

df <- bind_rows(P_predicted) %>%
  add_column(
    ID = database$ID,
    trip_n = database$trip_n + 1,
    choice = database$choice,
    .before = 0
  )

write_csv(
  df,
  "../../shared-svn/projects/matsim-berlin/data/SrV/2018/halil/output/mnl_CR_specification.csv"
)

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

###### use test data set
testData <- read_csv(
  "../../shared-svn/projects/matsim-berlin/data/SrV/2018/halil/Train-Test Data/test_dataAllAttributes.csv",
  comment = "#"
)

##apply same transformation 
test_database <- testData %>%
  rename(ID = person) %>%
  rename(util_money_inp = util_money)

distances_test <- test_database %>%
  group_by(ID) %>%
  summarise(total_dist = sum(beelineDist))


test_database <- left_join(test_database, distances_test, by = "ID") %>%
  mutate(
    distanceWeight = beelineDist / total_dist,
    distanceWeight = replace_na(distanceWeight, 1),
    car_dist_costs  = car_km  * -0.149,
    ride_dist_costs = ride_km * -0.149,
    car_fixed_costs = (car_valid == 1) * distanceWeight * -14.30,
    pt_fixed_costs  = (pt_valid  == 1) * distanceWeight * -3
  ) %>%
  mutate(ID = row_number())

##new test inputs
apollo_inputs_test <- apollo_validateInputs(database = test_database)

P_test <- apollo_probabilities(
  beta_hat,
  apollo_inputs_test,
  functionality = "prediction"
)


predictions_test <- data.frame(
  ID      = test_database$ID,
  trip_n  = test_database$trip_n,
  choice  = test_database$choice,
  prob_walk = P_test$model$walk,
  prob_pt   = P_test$model$pt,
  prob_car  = P_test$model$car,
  prob_bike = P_test$model$bike,
  prob_ride = P_test$model$ride
)

write_csv(
  predictions_test,
  "../../shared-svn/projects/matsim-berlin/data/SrV/2018/halil/output/test_predictions.csv"
)






