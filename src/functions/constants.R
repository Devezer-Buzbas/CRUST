################
##
## @description Define constants
##
## @param None
##
## @return None
##
## @lastChange 2018-03-22
##
## @changes
##   Change CARA to BO, NELL to MAVE, ROB to TESS [2018-03-22]
##   Included BIC output parameter [2018-03-22]
##   Included parameters output constants [2017-06-17]
##   Included the O_BETA1_TRUE output field constant [2017-06-08]
##   Change BOB to CARA, NEL to NELL, RAY to REY [2017-10-16]
##
################

## Types of agents
REY <- 1
TESS <- 2
BO <- 3
MAVE <- 4

## Model comparison
TSTATISTICS <- 1
RSQ <- 2
ARSQ <- 3
AIC <- 4
BIC <- 5

## Model output fields
O_NUM_FIELDS <- 22
O_STRATEGY <- 1
O_SELECTED_MODEL <- 2
O_SELECTED_TRUE_MODEL <- 3
O_SELECTED_MODEL_DISTANCE <- 4
O_INITIAL_GLOBAL_MODEL <- 5
O_INITIAL_GLOBAL_TRUE_MODEL <- 6
O_INITIAL_GLOBAL_MODEL_DISTANCE <- 7
O_FINAL_GLOBAL_MODEL <- 8
O_FINAL_GLOBAL_TRUE_MODEL <- 9
O_FINAL_GLOBAL_MODEL_DISTANCE <- 10
O_NUM_PREDICTORS <- 11
O_SAMPLE_SIZE <- 12
O_BETA1_TRUE <- 13
O_BETA1_ESTIMATE <- 14
O_BETA1_ERROR <- 15
O_TSTATISTICS <- 16
O_RSQ <- 17
O_ARSQ <- 18
O_AIC <- 19
O_BIC <- 20
O_REPLICATED <- 21
O_PREDICTORS <- 22

## File output headers
OUTPUT_HEADER <- c("replica",
                   "timestep",
                   "strategy",
                   "selected_model",
                   "selected_true_model",
                   "selected_model_distance",
                   "initial_global_model",
                   "initial_global_true_model",
                   "initial_global_model_distance",
                   "final_global_model",
                   "final_global_true_model",
                   "final_global_model_distance",
                   "num_predictors",
                   "sample_size",
                   "beta1_true",
                   "beta1_estimate",
                   "beta1_error",
                   "tStatistics",
                   "RSQ",
                   "ARSQ",
                   "AIC",
                   "BIC",
                   "replicated",
                   "predictors")

## Parameter output
P_TMODEL <- 1
P_K <- 2
P_SAMPLE_SIZE <- 3
P_SIGMA <- 4
P_CORRELATION <- 5
P_AGENT_WEIGHTS <- 6
P_TRUE_BETAS <- 7
P_XSET <- 8
