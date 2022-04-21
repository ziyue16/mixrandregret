library(haven)
library(apollo)
library(readxl)
#### LOAD DATA ####

database =as.data.frame( read_excel("data.xls"))

# ################################################################# #
#### LOAD LIBRARY AND DEFINE CORE SETTINGS                       ####
# ################################################################# #

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName       = "MMNL_preference_space",
  modelDescr      = "Mixed random regret logit model on simulated data",
  indivID         = "id_ind",  
  mixing          = TRUE,
  nCores          = 1,
  outputDirectory = NULL
)


### Vector of parameters, including any that are kept fixed in estimation
apollo_beta = c(mu_x1    = 1   ,
                sd_x1    = 0.25,
                mu_x2    = 2   ,
                sd_x2    = 0.25,
                b_fix  = 0.01  )

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c()




# ################################################################# #
#### DEFINE RANDOM COMPONENTS                                    ####
# ################################################################# #

### Set parameters for generating draws
apollo_draws = list(
  interDrawsType = "halton",
  interNDraws    = 1000,
  interUnifDraws = c(),
  interNormDraws = c("draws_x1","draws_x2"),
  intraDrawsType = "halton",
  intraNDraws    = 0,
  intraUnifDraws = c(),
  intraNormDraws = c()
)


### Create random parameters
apollo_randCoeff = function(apollo_beta, apollo_inputs){
  randcoeff = list()
  
  randcoeff[["b_x1"]] =  mu_x1 + sd_x1 * draws_x1 
  randcoeff[["b_x2"]] =  mu_x2 + sd_x2 * draws_x2 

  return(randcoeff)
}


apollo_inputs = apollo_validateInputs()



# ################################################################# #
#### DEFINE MODEL AND LIKELIHOOD FUNCTION                        ####
# ################################################################# #

apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  ### Function initialisation: do not change the following three commands
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### x1 differences
# dx1_1 = ( x11 -  x11 ) 
  dx1_2_1 = ( x12 -  x11 ) 
  dx1_3_1 = ( x13 -  x11 ) 

  dx1_1_2 = ( x11 -  x12 ) 
# dx1_2_2 = ( x12 -  x12 ) 
  dx1_3_2 = ( x13 -  x12 ) 

  dx1_1_3 = ( x11 -  x13 ) 
  dx1_2_3 = ( x12 -  x13 ) 
#  dx1_3_3 = ( x13 -  x13 ) 
  

  ### x2 differences
# dx2_1_1 = ( x21 -  x21 ) 
  dx2_2_1 = ( x22 -  x21 ) 
  dx2_3_1 = ( x23 -  x21 ) 
  
  dx2_1_2 = ( x21 -  x22 ) 
# dx2_2_2 = ( x22 -  x22 ) 
  dx2_3_2 = ( x23 -  x22 ) 
  
  dx2_1_3 = ( x21 -  x23 ) 
  dx2_2_3 = ( x22 -  x23 ) 
# dx2_3_3 = ( x23 -  x23 ) 
  
### x_fix  differences
# dfix_1_1 = ( x_fix1 -  x_fix1 ) 
  dfix_2_1 = ( x_fix2 -  x_fix1 ) 
  dfix_3_1 = ( x_fix3 -  x_fix1 ) 
  
  dfix_1_2 = ( x_fix1 -  x_fix2 ) 
# dfix_2_2 = ( x_fix2 -  x_fix2 ) 
  dfix_3_2 = ( x_fix3 -  x_fix2 ) 
  
  dfix_1_3 = ( x_fix1 -  x_fix3 ) 
  dfix_2_3 = ( x_fix2 -  x_fix3 ) 
# dfix_3_3 = ( x_fix3 -  x_fix3 ) 
  
  
  
  ### Create list of probabilities P
  P = list()
  
  ### List of utilities: these must use the same names as in mnl_settings, order is irrelevant
  V = list()
  V[["alt1"]] = - log(1+exp(b_x1  * dx1_2_1  ))  - log(1+exp(b_x1 * dx1_3_1 )) + 
                - log(1+exp(b_x2  * dx2_2_1  ))  - log(1+exp(b_x2 * dx2_3_1 )) +     
                - log(1+exp(b_fix *dfix_2_1 )) - log(1+exp(b_fix * dfix_3_1 ))  
    
  V[["alt2"]] = - log(1+exp(b_x1  * dx1_1_2  ))  - log(1+exp(b_x1 * dx1_3_2 )) + 
                - log(1+exp(b_x2  * dx2_1_2  ))  - log(1+exp(b_x2 * dx2_3_2 )) +     
                - log(1+exp(b_fix *dfix_1_2 )) - log(1+exp(b_fix * dfix_3_2 ))  
  
  V[["alt3"]] = - log(1+exp(b_x1  * dx1_1_3  ))  - log(1+exp(b_x1 * dx1_2_3 )) + 
                - log(1+exp(b_x2  * dx2_1_3  ))  - log(1+exp(b_x2 * dx2_2_3 )) +     
                - log(1+exp(b_fix *dfix_1_3 )) - log(1+exp(b_fix * dfix_2_3 ))  
  
  
  
  ### Define settings for MNL model component
  mnl_settings = list(
    alternatives  = c(alt1=1, alt2=2, alt3=3),
    avail         = list(alt1=1, alt2=1, alt3=1),
    choiceVar     = choice_wide,
    utilities     = V
  )
  
  ### Compute probabilities using MNL model
  P[["model"]] = apollo_mnl(mnl_settings, functionality)
  
  ### Take product across observation for same individual
  P = apollo_panelProd(P, apollo_inputs, functionality)
  
  ### Average across inter-individual draws
  P = apollo_avgInterDraws(P, apollo_inputs, functionality)
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

model = apollo_estimate(apollo_beta, apollo_fixed,apollo_probabilities, apollo_inputs)

# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO SCREEN)                               ----
# ----------------------------------------------------------------- #
apollo_modelOutput(model)
# 
# 
# Model run by u0133260 using Apollo 0.2.7 on R 4.1.2 for Windows.
# www.ApolloChoiceModelling.com
# 
# Model name                       : MMNL_preference_space
# Model description                : Mixed random regret logit model on simulated data
# Model run at                     : 2022-04-19 17:58:14
# Estimation method                : bfgs
# Model diagnosis                  : successful convergence 
# Number of individuals            : 500
# Number of rows in database       : 3000
# Number of modelled outcomes      : 3000
# 
# Number of cores used             :  10 
# Number of inter-individual draws : 1000 (halton)
# 
# LL(start)                        : -744.14
# LL(0)                            : -3295.84
# LL(C)                            : -3293.8
# LL(final)                        : -714.27
# Rho-square (0)                   :  0.7833 
# Adj.Rho-square (0)               :  0.7818 
# Rho-square (C)                   :  0.7831 
# Adj.Rho-square (C)               :  0.7816 
# AIC                              :  1438.54 
# BIC                              :  1468.57 
# 
# Estimated parameters             :  5
# Time taken (hh:mm:ss)            :  00:03:15.84 
# pre-estimation              :  00:01:8.21 
# estimation                  :  00:00:49.28 
# post-estimation             :  00:01:18.35 
# Iterations                       :  15  
# Min abs eigenvalue of Hessian    :  76.59837 
# 
# Unconstrained optimisation.
# 
# Estimates:
#   Estimate        s.e.   t.rat.(0)    Rob.s.e. Rob.t.rat.(0)
# mu_x1    0.885367     0.05194     17.0473     0.05477       16.1655
# sd_x1    0.359805     0.03572     10.0725     0.04293        8.3811
# mu_x2    1.808286     0.09631     18.7764     0.09757       18.5327
# sd_x2   -0.097858     0.08669     -1.1288     0.11758       -0.8323
# b_fix    0.008325     0.04017      0.2072     0.04067        0.2047


conditionals = apollo_conditionals(model,
                                   apollo_probabilities,
                                   apollo_inputs)


(posterior_x1 <- hist(conditionals$b_x1$post.mean))

(posterior_x2 <-hist(conditionals$b_x2$post.mean))


