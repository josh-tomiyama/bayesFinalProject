library(ABSEIR)
library(splines)
library(clipr)
library(dplyr)
library(kableExtra)
library(knitr)
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
plotR0 = function(simulations, main)
{
  allSimulatedEA_R0 = sapply(simulations$simulationResults, function(x){x$R_EA})
  plot(apply(allSimulatedEA_R0, 1, mean), type = "l", ylim = c(0, 3), lwd =2,
       ylab = "Reproductive Number", main = main)
  lines(apply(allSimulatedEA_R0, 1, mean), lwd = 2, lty = 2, col = "blue")
  lines(apply(allSimulatedEA_R0, 1, quantile, probs = c(0.1)), 
        lwd = 2, lty = 2, col = "blue")
  lines(apply(allSimulatedEA_R0, 1,  quantile, probs = c(0.9)), 
        lwd = 2, lty = 2, col = "blue")
}

ciPlot <- function(eta, alpha = 0.5, ...){
  eta_lwr = apply(eta, 1, quantile, prob = .025)
  eta_upr = apply(eta, 1, quantile, prob = .975)
  eta_lwr2 = apply(eta, 1, quantile, prob = 1/6)
  eta_upr2 = apply(eta, 1, quantile, prob = 1 - 1/6)
  x_idx <- 1:length(eta_lwr)
  plot(apply(eta, 1, mean), 
       type = 'l', 
       ylim = c(min(eta_lwr), max(eta_upr)),
       ...)
  polygon(c(x_idx, rev(x_idx)),
          c(eta_lwr2, rev(eta_upr2)),
          col = adjustcolor('gray', alpha.f = alpha),
          border = NA)
  lines(eta_lwr, lty = 2)
  lines(eta_upr, lty = 2)
}

df2Clip <- function(df, ...){
  library(clipr)
  df <- as.data.frame(df)
  kable(df, ...) %>%
    kable_styling(full_width = FALSE) %>%
    write_clip()
}



#-------------------------------------------------------------------------------

TBdat <- read.csv("./Data/fitting_data")
true_counts <- TBdat$TUBERCULOSISCURRENTQUARTER

time_idx <- 1:length(true_counts)
X <- cbind(1, bs(time_idx,3))
X_seasonal <- cbind(1, 
                    sin(2*pi*time_idx/4), 
                    cos(2*pi*time_idx/4), 
                    cos(2*pi*time_idx/4)*sin(2*pi*time_idx/4))
result1 <- readRDS("./ModelCache/nonlinearfit.RDS")
result2 <- readRDS("./ModelCache/seasonalfit.RDS")

#beta---------------------------------------------------------------------------
beta_idx <- grepl("Beta_SE", colnames(result1$param.samples))
eta <- X %*% t(result1$param.samples[,beta_idx])

beta_idx_seasonal <- grepl("Beta_SE", colnames(result2$param.samples))
eta_seasonal <- X_seasonal %*% t(result2$param.samples[,beta_idx_seasonal])
plot(apply(eta_seasonal, 1, mean), 
     type = 'l', 
     main = 'Seasonal Model Linear Predictor',
     ylab = 'XB',
     xlab = 'Time Index')




jpeg("./PresentationImages/ExposuresPlots.jpg",
     quality = 100)
par(mfrow = c(2,1))
ciPlot(eta, 
       main = 'Nonlinear Model Linear Predictor',
       ylab = 'XB',
       xlab = '')
# lines(eta_lwr2, lty = 2)
# lines(eta_upr2, lty = 2)

ciPlot(eta_seasonal, 
       main = 'Seasonal Model Linear Predictor',
       ylab = 'XB',
       xlab = 'Time Index')
dev.off()

#bf---------------------------------------------------------------------------

tedf <- data.frame('Terminating Epsilon' = c(result1$current_eps, 
                                     result2$current_eps))
rownames(tedf) <- c("Nonlinear", "Seasonal")
kable(tedf, digits = 3, caption = "Terminating Epsilon of models",
      format = 'latex', booktabs = T) %>%
  kable_styling(full_width = FALSE) %>%
  write_clip()

#transition---------------------------------------------------------------------

gamma_seasonal_id <- grepl("gamma_", colnames(result2$param.samples))
betars_seasonal_id <- grepl("Beta_RS", colnames(result2$param.samples))
gamma_params <- result2$param.samples[,gamma_seasonal_id]
beta_rs <- result2$param.samples[,betars_seasonal_id]
pi_trans <- cbind(1 - exp(-gamma_params),'beta_rs' = 1 - exp(-exp(beta_rs)))

apply(pi_trans, 2, function(x)(c(mean(x), sd(x), 
                                 quantile(x, probs = c(0.025, 0.975))
                                 )
                               )
      ) %>% t()

df2Clip(sres$parameterEstimates, digits = 3, caption = "Transition Parameter Estimates",
        format = 'latex', booktabs = T)

sres$parameterEstimates %>%
  kable(digits = 3, 
        caption = "Transition Parameters of Seasonal Model",
        format = 'latex', booktabs = T) %>%
  kable_styling(full_width = FALSE)
plot(apply(eta_seasonal, 1, mean), 
     type = 'l', 
     main = 'Seasonal Model Linear Predictor',
     ylab = 'XB',
     xlab = 'Time Index')

simulations1 <- readRDS("./ModelCache/sim1.RDS")
simulations2 <- readRDS("./ModelCache/sim2.RDS")

bf <- readRDS("./ModelCache/bf.RDS")
bf

##ppred------------------------------------------------------------------------

jpeg("./PresentationImages/PostPredPlots.jpg",
     quality = 100)
par(mfrow = c(2,1))
plotPosteriorPredictive(simulations1, true_counts, 
                        "Posterior Predictive Distribution nonlinear", "I_star")
plotPosteriorPredictive(simulations2, true_counts, 
                        "Posterior Predictive Distribution seasonal", "I_star")
dev.off()

#r0-----------------------------------------------------------------------------
simulations1.R0 <- readRDS("./ModelCache/R0_mod1.RDS")
simulations2.R0 <- readRDS("./ModelCache/R0_mod2.RDS")

plotR0(simulations1.R0, "Model 1: EA-R(t)")
plotR0(simulations2.R0, "Model 2: EA-R(t)")



#Probability peak weak is 22----------------------------------------------------
allSimulatedCompartment = sapply(simulations2$simulationResults, function(x){x[['I_star']]})
pkweek = apply(allSimulatedCompartment, 2, function(x){max(x[21:23]) == x[22]})
mean(pkweek)
