#### IV-MBM Monte Carlo ###
# Scott Coffin, Ph.D. #
# 6-27-2024 #

# Load necessary libraries
library(data.table)
library(tidyverse)
library(httk)


########## ----- Check to see if IV-MBM model works as expected with simple non-stochastic examples ---###
### copied from Vf_Honda2019.Rmd script ###

# # Source the provided R script
# source("httk/R/armitage.R")
# 
# armitage.dt <- copy(armitage_input)
# chem_phys <- copy(chem.physical_and_invitro.data)

# # join phys-chem data
# armitage.dt2 <- left_join(armitage.dt %>% select(-MW), 
#                           chem_phys, by = c("casrn" = "CAS")) %>% 
#   dplyr::mutate(CAS = casrn,
#                 chem.name = compound_name,
#                 chem.cas = casrn)
# 
# 
# #define user input for modelling
# armitage.dt3 <- armitage.dt2[,well_number:=384] %>% 
#   .[,option.bottom:=TRUE] %>% 
#   .[,option.plastic:=TRUE] %>% 
#   .[,Tsys:=37] %>% 
#   .[,Tref:=298.15] %>% 
#   .[,FBSf:=0.1] %>% 
#   .[,nomconc:=50]
# 
# ## example ##
# armitage.dt3 <- armitage_estimate_sarea(tcdata = armitage.dt3)
# 
# armitage.dt4 <- subset(armitage.dt3,casrn%in%get_cheminfo())
# 
# armitage_output1 <- armitage_eval(tcdata = armitage.dt4[,ac50:=50])
# armitage_output2 <- armitage_eval(tcdata = armitage.dt3[,ac50:=1])
# armitage_output3 <- armitage_eval(tcdata = armitage.dt3[,ac50:=0.001])

### copied from Armitage.R script ###

# Check to see if we have info on the chemical:
"80-05-7" %in% get_cheminfo()

#We do:
temp <- armitage_eval(casrn.vector = c("80-05-7"#, 
                                       #"81-81-2"
                                       ), this.FBSf = 0.1,
                      this.well_number = 384, nomconc = 10)
print(temp$cfree.invitro)

### does not work ###

# Check to see if we have info on the chemical:
"793-24-8" %in% get_cheminfo()

# Since we don't have any info, let's look up phys-chem from dashboard:
cheminfo <- data.frame(
  Compound="6-PPD",
  CASRN="793-24-8",
  DTXSID="DTXSID9025114",
  logP=4.27,
  logHenry=log10(7.69e-8),
  logWSol=log10(1.58e-4),
  MP=	99.4,
  MW=268.404
  )

# Add the information to HTTK's database:
chem.physical_and_invitro.data <- add_chemtable(
 cheminfo,
 current.table=chem.physical_and_invitro.data,
 data.list=list(
 Compound="Compound",
 CAS="CASRN",
  DTXSID="DTXSID",
  MW="MW",
  logP="logP",
  logHenry="logHenry",
  logWSol="logWSol",
  MP="MP"),
  species="Human",
  reference="CompTox Dashboard 31921")

# Run the Armitage et al. (2014) model:
out <- armitage_eval(
  casrn.vector = "793-24-8",
  this.FBSf = 0.1,
  this.well_number = 384,
  nomconc = 10)

# # attempt on multiple
# chem.physical_and_invitro.data %>% 
#   drop_na(logHenry, logWSol, MP)
# 
# out_mulitple <- armitage_eval(
#   casrn.vector = c("793-24-8", "80-05-7"),
#   this.FBSf = 0.1,
#   this.well_number = 384,
#   nomconc = 10)

### does not work ##

######## -------------- Monte Carlo Simulation ---------------- ##############
# Define the number of simulations
n_sim <- 1000

# Define a function to generate simulated data
generate_simulated_data <- function(n_sim) {
  data <- data.table(
    casrn = "793-24-8",  # Example CAS number: replace with actual
    nomconc = rnorm(n_sim, mean = 10, sd = 2),  # Example nominal concentration
    this.FBSf = runif(n_sim, min = 0.1, max = 0.15),  # Example FBS fraction
    this.Tsys = rnorm(n_sim, mean = 37, sd = 1), #example system temp 
    this.pH = rnorm(n_sim, mean = 7, sd = 0.1), #example system pH
    this.P_cells = rnorm(n_sim, mean = 1, sd = 0.1) 
  )
  return(data)
}

# Generate the simulated data
simulated_data <- generate_simulated_data(n_sim)

# Define a function to run IV-MBM on each simulated dataset
run_simulation <- function(row) {
  result <- armitage_eval(
    chem.cas = row$casrn,
    this.FBSf = row$this.FBSf,
    this.Tsys = row$this.Tsys,
    this.pH = row$this.pH,
    this.P_cells = row$this.P_cells
  )
  
  return(result$cfree.invitro)
}

# Run the Monte Carlo simulations
results <- simulated_data %>%
  rowwise() %>%
  mutate(cfree_invitro = run_simulation(cur_data()))

# Summarize the results
summary_results <- results %>%
  group_by(casrn) %>% 
  summarise(
    mean_cfree = mean(cfree_invitro, na.rm = TRUE),
    sd_cfree = sd(cfree_invitro, na.rm = TRUE),
    median_cfree = median(cfree_invitro, na.rm = TRUE),
    ci_lower = quantile(cfree_invitro, 0.025, na.rm = TRUE),
    ci_upper = quantile(cfree_invitro, 0.975, na.rm = TRUE)
  )

print(summary_results)

# Plot the distribution of cfree_invitro
ggplot(results, aes(x = cfree_invitro)) +
  geom_histogram(fill = "blue", alpha = 0.7) +
  labs(title = "Distribution of Free Concentration in Vitro (cfree.invitro)",
       x = "Free Concentration (cfree.invitro)",
       y = "Frequency") +
  theme_minimal()
