library(ABSEIR)
library(splines)

TBdat <- read.csv("./Data/fitting_data")
true_counts <- TBdat$TUBERCULOSISCURRENTQUARTER
counts <- rep(NA, length = 6*4)
## Leave out last year
counts[1:(nrow(TBdat) - 3)] <- TBdat$TUBERCULOSISCURRENTQUARTER[1:(nrow(TBdat) - 3)]

data_model = DataModel(counts,
                       type = "identity",
                       compartment="I_star",
                       cumulative=FALSE)

data_model2 = DataModel(counts,
                       type = "identity",
                       compartment="R_star",
                       cumulative=FALSE)

time_idx <- 1:length(counts)
X <- cbind(1, bs(time_idx,3))
exposure_model_1 = ExposureModel(X,
                                 nTpt = length(counts),
                                 nLoc = 1,
                                 betaPriorPrecision = 0.5,
                                 betaPriorMean = 0)

X_seasonal <- cbind(1, 
                    sin(2*pi*time_idx/4), 
                    cos(2*pi*time_idx/4), 
                    cos(2*pi*time_idx/4)*sin(2*pi*time_idx/4))
exposure_model_2 = ExposureModel(X_seasonal,
                                 nTpt = length(counts),
                                 nLoc = 1,
                                 betaPriorPrecision = 0.5,
                                 betaPriorMean = 0)

# intervention_term = cumsum(Kikwit1995$Date >  as.Date("05-09-1995", "%m-%d-%Y"))
# exposure_model_2 = ExposureModel(cbind(1,intervention_term),
#                                  nTpt = nrow(Kikwit1995),
#                                  nLoc = 1,
#                                  betaPriorPrecision = 0.5,
#                                  betaPriorMean = 0)


## https://www.atsjournals.org/doi/10.1164/rccm.200409-1200OC?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed
## # of recurrences/person-years
reinfect_rate = (48+13+47)/(1794 + 466 + 47)
## crudely convert to person-quarters
reinfect_rate = reinfect_rate/4
pt_out <- poisson.test(48+13+47, 4*(1794 + 466 + 47))
beta_rs_mean <- -log(1-pt_out$estimate)
beta_rs_ci <- -log(1 - pt_out$conf.int)
## 80% chance between
f <- function(log_sd, mu, ci){
  (qnorm(.1, mu, exp(log_sd)) - ci[1])^2 +  
    (qnorm(.9, mu, exp(log_sd)) - ci[2])^2
}
opt_out <- optim(1, 
                 f, 
                 mu = beta_rs_mean, 
                 ci = beta_rs_ci,
                 method = 'Brent',
                 lower = -10, 
                 upper = 10)
beta_rs_precision <- 1/(exp(opt_out$par)^2)

reinfection_model = ReinfectionModel("SEIRS",
                                     X_prs = rep(1, length(counts)),
                                     priorPrecision = beta_rs_precision,
                                     priorMean = beta_rs_mean
                                     )

distance_model = DistanceModel(list(matrix(0)))

### estimates of initial values
### 2013 3/100,000 = 30/million
### https://www.cdc.gov/tb/statistics/tbcases.htm
initial_value_container = InitialValueContainer(S0=300e6,
                                                E0=1000,
                                                I0=9000/4,
                                                R0=3*9000/4)

### Formula: 1 - exp(-rate*time) = quantile
### priors: E -> I 90% chance of moving by 9 months (3 quarters)
### I -> R 90% chance of dying or being treated within 1 year (4 quarters)
transition_priors = ExponentialTransitionPriors(p_ei = 1 - exp(log(1-0.9)/3), 
                                                p_ir= 1 - exp(log(1-0.9)/4),
                                                p_ei_ess = 100,
                                                p_ir_ess = 10)

sampling_control = SamplingControl(seed = 123123, 
                                   n_cores = 16,
                                   algorithm="Beaumont2009",
                                   list(batch_size = 10000,
                                        epochs = 1e6,
                                        max_batches = 1000,
                                        shrinkage = 0.99,
                                        multivariate_perturbation=FALSE,
                                        keep_compartments = TRUE
                                   ) 
)

runtime1 = system.time(result1 <- SpatialSEIRModel(data_model,
                                                   exposure_model_1,
                                                   reinfection_model,
                                                   distance_model,
                                                   transition_priors,
                                                   initial_value_container,
                                                   sampling_control,
                                                   samples = 1000,
                                                   verbose = 0))
saveRDS(result1, file = "./ModelCache/nonlinearfit.RDS")

# Reasonable intensity
runtime2 = system.time(result2 <- SpatialSEIRModel(data_model,
                                                   exposure_model_2,
                                                   reinfection_model,
                                                   distance_model,
                                                   transition_priors,
                                                   initial_value_container,
                                                   sampling_control,
                                                   samples = 1000,
                                                   verbose = 0))
saveRDS(result2, file = "./ModelCache/seasonalfit.RDS")

# runtime2 = system.time(result2 <- SpatialSEIRModel(data_model2,
#                                                    exposure_model_1,
#                                                    reinfection_model,
#                                                    distance_model,
#                                                    transition_priors,
#                                                    initial_value_container,
#                                                    sampling_control,
#                                                    samples = 100,
#                                                    verbose = 0))
# saveRDS(result2, file = "model2fit.RDS")


beta_idx <- grepl("Beta_SE", colnames(result1$param.samples))
eta <- X %*% t(result1$param.samples[,beta_idx])
plot(apply(eta, 1, mean))

simulations1 <- epidemic.simulations(result1, replicates = 50)
simulations2 <- epidemic.simulations(result2, replicates = 50)

saveRDS(simulations1, "./ModelCache/sim1.RDS")
saveRDS(simulations2, "./ModelCache/sim2.RDS")

bf <- compareModels(list(result1, result2))
saveRDS(bf, "./ModelCache/bf.RDS")

plotPosteriorPredictive = function(simulations, obs_count, main, xlab = "Time Index", compartment = "I_star")
{
  allSimulatedCompartment = sapply(simulations$simulationResults, function(x){x[[compartment]]})
  
  lowerQuantile = apply(allSimulatedCompartment, 1, quantile, probs = c(0.025))
  posteriorMean = apply(allSimulatedCompartment, 1, mean)
  upperQuantile = apply(allSimulatedCompartment, 1, quantile, probs = c(0.975))
  
  
  plot(obs_count, ylim = c(0, max(obs_count, na.rm = TRUE)*2),
       xlab = xlab, ylab = "New Cases", main = main)
  lines(upperQuantile, lty = 2, col = "blue")
  lines(lowerQuantile, lty = 2, col = "blue")
  lines(posteriorMean, lty = 1, col = "blue")
  
  legend(x = 100, y = 12, legend = c("Mean", "95% CI", "Observed"), lty = c(1,2,0), 
         pch = c(NA,NA,1), col = c("blue", "blue", "black"), cex = 1)
}
plotPosteriorPredictive(simulations1, counts, "Posterior Predictive Distribution I model", "I_star")
plotPosteriorPredictive(simulations2, counts, "Posterior Predictive Distribution R model", "R_star")
test_data <- read.csv("./Data/test_data.csv")
points(29:32, test_data$count)

Year_cum <- function(simulations, obs_count, main, xlab = "Time Index", compartment = "I_star")
{
  # browser()
  allSimulatedCompartment = sapply(simulations$simulationResults, function(x){
    matrix(colSums(matrix(x[[compartment]], nrow = 4)))
    })
  
  lowerQuantile = apply(allSimulatedCompartment, 1, quantile, probs = c(0.025))
  posteriorMean = apply(allSimulatedCompartment, 1, mean)
  upperQuantile = apply(allSimulatedCompartment, 1, quantile, probs = c(0.975))
  
  
  plot(obs_count, ylim = c(0, max(obs_count, na.rm = TRUE)*2),
       xlab = xlab, ylab = "New Cases", main = main)
  lines(upperQuantile, lty = 2, col = "blue")
  lines(lowerQuantile, lty = 2, col = "blue")
  lines(posteriorMean, lty = 1, col = "blue")
  
  legend(x = 100, y = 12, legend = c("Mean", "95% CI", "Observed"), lty = c(1,2,0), 
         pch = c(NA,NA,1), col = c("blue", "blue", "black"), cex = 1)
}

# Year_cum(simulations1, colSums(matrix(counts, nrow = 4), na.rm = TRUE), "Posterior Predictive Distribution I model", "I_star")




simulations1.R0 <- ComputeR0(simulations1, cores = 16)
simulations2.R0 <- ComputeR0(simulations2, cores = 16)
saveRDS(simulations1.R0, "./ModelCache/R0_mod1.RDS")
saveRDS(simulations2.R0, "./ModelCache/R0_mod2.RDS")

plotR0 = function(simulations, main)
{
  allSimulatedEA_R0 = sapply(simulations$simulationResults, function(x){x$R_EA})
  plot(apply(allSimulatedEA_R0, 1, mean), type = "l", ylim = c(0, 3), lwd =2,
       ylab = "Reproductive Number", main = main)
  lines(apply(allSimulatedEA_R0, 1, mean), lwd = 2, lty = 2, col = "blue")
  lines(apply(allSimulatedEA_R0, 1, quantile, probs = c(0.1)), lwd = 2, lty = 2, col = "blue")
  lines(apply(allSimulatedEA_R0, 1,  quantile, probs = c(0.9)), lwd = 2, lty = 2, col = "blue")
}
plotR0(simulations1.R0, "Model 1: EA-R(t)")
