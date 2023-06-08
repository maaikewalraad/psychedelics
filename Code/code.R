rm(list = ls())
require(data.table)
require(magrittr) # special pipe operators
require(psych) # describe()
require(dplyr) # data wrangling
require(knitr) 
require(ggpubr)

set.seed(7079540) # for reproducibility

# Import data
dt <- read.csv("../Data/AHRI_DATASET_PM_MANUSCRIPT_DATA.csv")
dt %<>% as.data.table
md <- read.csv("../Data/AHRI_DATASET_PM_MANUSCRIPT_CODEBOOK.csv") # metadata
md %<>% as.data.table

# Data wrangling ----------------------------------------------------------
dt %<>% select(PHQ9_SCORE, PM1_DIAG_CONDITION, C_DP, C_TOTAL, CCI_SCORE) # keep main variables only
str(dt)

# Inclusion criteria include:
dt %<>% filter(C_DP == 1) # experienced depression last year
dt %<>% filter(PM1_DIAG_CONDITION != -99) # keep category 0 and 1 only (no need to use C_DP since these people are *diagnosed* with some psychiatric disorder)
dt$PM1_DIAG_CONDITION <- factor(dt$PM1_DIAG_CONDITION, levels = c(0, 1))
md[VAR == "PHQ9_GE10"]$Description
#dt %<>% filter(PHQ9_GE10 == 1) # "while most patients (88%) with major depression had scores of 10 or greater."
dt %<>% filter(PHQ9_SCORE >= 15) # most patients with PHQ9 score of 15 or higher get treatment

# Descriptive statistics -----------------------------------------------------
describeBy(dt$PHQ9_SCORE, dt$PM1_DIAG_CONDITION)

# Fit standard linear model -----------------------------------------------
mod <- lm(PHQ9_SCORE ~ PM1_DIAG_CONDITION, data = dt)
summary(mod)

# Assumptions -------------------------------------------------------------

# Linearity
plot(x = dt$C_TOTAL, y = dt$PHQ9_SCORE)

# homoscedasticity (constant residual variance)
plot(mod)[1] 

# normality in the outcome
dt %$% hist(PHQ9_SCORE) # DEPENDENT VARIABLE  

# normality in residuals
plot(mod)[2] # QQ-plot shows divergence from normality at the tails

# Gibbs sampling ----------------------------------------------------------

# First step: specify priors for each parameter. 

# # Uninformative priors (compare with standard lm)
# # b0 
# mu00 <- 1 
# zeta00 <- 1.0E4  
# # b1 
# mu10 <- 1
# zeta10 <- 1.0E4
# nu10 <- 100 # prior df's -> bigger = wider tails
# # variance (σ^2) 
# a_0 = 1 
# b_0 = 1
# sig2_0 = 1/rgamma(1, shape = a_0, rate = b_0) 

# Informative priors. These are based on the historical data (see manuscript)
# b0 
mu00 <- 0 # Q. make this based on data or whatever
zeta00 <- 1.0E4 
# b1 
mu10 <- -11.7 # prior mean
zeta10 <- 3^2 # prior variance 
nu10 <- 25 # prior df's -> bigger = wider tails
# variance (σ^2) 
a_0 = 1.0E4 
b_0 = 1.0E4 # large values for scale and shape denote strong belief 

# Second step: the Gibbs algorithm

post_b1 <- function(b0, b1, sig2, y, x1) { # MH function for b1
  I <- -(b1^2) * ((sum(x1^2))/(2*sig2)) + b1 * ((sum(x1*(y - b0))) / (sig2))
  II <- (1 + ((b1 - mu10)^2)/(nu10*zeta10))^((nu10+1)/2)
  r_post <- exp(I) * II
  return(r_post) 
}

gibbs.chains <- function(b0, b1, sig2, y, x1) {
  warmup <- 1000
  H <- 25000 + warmup # number of samples drawn 
  n_par <- 3
  simulated <- matrix(0, nrow = H, ncol = n_par)
  colnames(simulated) <- c("b0","b1", "sig2") 
  simulated[1, ] <- c(b0, b1, sig2)
  
  for (h in 2:H) { # determine full conditional distribution for...

    # b0:
    mu0 <- (sum(y-b1*x1)/sig2 + (mu00/zeta00))/((nrow(dt) / sig2) + (1 / zeta00))  
    sd0 <- sqrt(1 / ((nrow(dt) / sig2) + (1 / zeta00)))
    b0 <- rnorm(1, mu0, sd0) # sample from distribution
    
    # b1 (incl. MH, Non-conjugate prior):
    beta1_c <- b1 
    beta1_n <- rnorm(1, mean = 0, sd = 3)   # Sample from proposal density
    r_post <- post_b1(b0, beta1_n, sig2, y, x1)/post_b1(b0, beta1_c, sig2, y, x1) # Ratio of posterior distributions
    r_prop <- beta1_c / beta1_n   # Ratio of proposal distribution
    r <- r_post * r_prop # Acceptance ratio r
    u <- runif(1, min = 0, max = 1)   # 'Sample a probability' from the uniform distribution
    ifelse(u <= abs(r), b1 <- beta1_n, b1 <- beta1_c) # retain of reject newly sampled value
    
    # # b1 (exc. MH)
    # mu1 <- (sum(x1*(y-b0)) / sig2 + (mu10/zeta10)) / ((sum(x1^2) / sig2) + (1 / zeta10))
    # sd1 <- sqrt(1 / ((sum(x1^2) / sig2) + (1 / zeta10)))
    # b1 <- rnorm(1, mu1, sd1)
    
    # variance: 
    a_1 <- (nrow(dt) / 2) + a_0 # posterior shape
    S <- sum(y - (b0 + b1*x1))^2 # RSS
    b_1 <- (S / 2) + b_0 # posterior scale 
    sig2 <- 1 / rgamma(1, a_1, b_1)
    
    simulated[h, ] <- c(b0, b1, sig2)  # Save them all 
    
  }
  
  return(simulated)
  
}

# Third: obtain the Gibbs estimates

dt$PM1_DIAG_CONDITION <- ifelse(dt$PM1_DIAG_CONDITION == 1, 1, 0) # make numeric
# Choose starting values for each of the parameters
chain1 <- gibbs.chains(b0 = 10, b1 = 0.3, sig2 = 0.5, y = dt$PHQ9_SCORE, x1 = dt$PM1_DIAG_CONDITION)
chain2 <- gibbs.chains(b0 = 2, b1 = 3, sig2 = 5, y = dt$PHQ9_SCORE, x1 = dt$PM1_DIAG_CONDITION) 


# Estimates ---------------------------------------------------------------

chain1 <- as.data.table(chain1[1001:nrow(chain1), ]) # remove warmup
chain2 <- as.data.table(chain2[1001:nrow(chain2), ]) # and transform into datatable

gibbs_stats <- function(chain) { # requires a datatable, with columns as parameters
  
  perc <- c(0.025, 0.975)
  x <- chain[, lapply(.SD, mean)] # EAP
  y <- chain[, lapply(.SD, sd)] # SD
  z <- chain[, lapply(.SD, function(x) quantile(x, perc))] # CI
  estimates <- cbind(t(x), t(y), t(z))
  rownames(estimates) <- c("Intercept", "PM", "Variance")
  colnames(estimates) <- c("Beta", "SD", "2.5% ", "97.5%")
  return(estimates)
  
}

gibbs_stats(chain1) %>% write.csv(., file = "../Output/Estimates_chain1.csv", row.names = F)
gibbs_stats(chain2) %>% write.csv(., file = "../Output/Estimates_chain2.csv", row.names = F)

# Assess convergence ------------------------------------------------------

# Traceplots
traceplot <- function(chain1, chain2) {   # only handles 2 chains (per parameter)
  
  par <- colnames(chain1)
  H <- nrow(chain1)
  
  for (i in 1:length(par)) {
    
    ggplot() + 
      geom_line(chain1, mapping = aes_string(x = I(1:H), y = par[i]), color = "#CC3399", linewidth = 0.4) +  # alpha = 0.8
      geom_line(chain2, mapping = aes_string(x = I(1:H), y = par[i]), color = "#0000FF", linewidth = 0.4) +  # alpha = 0.6
      ggtitle("Convergence plot") +
      ylab("Sampled value") +
      xlab("Number of iterations") +
      theme(
        plot.title = element_text(color = "black", size = 11, face = "bold.italic", hjust = 0.45),
        axis.title.x = element_text(color = "#333333", size = 10, face = "bold"),
        axis.title.y = element_text(color = "#333333", size = 10, face = "bold"))

  ggsave(paste0("../Output/traceplot_", par[i], ".png"), width = 8, height = 4)
  }
  
  
}

traceplot(chain1, chain2)


# Autocorrelation plots
autocors <- function(chain) {    # only handles 2 chains (per parameter)
  
  lags <- 40 # number of lags
  par <- colnames(chain)
  
  autocor <- matrix(0, nrow = lags, ncol = length(par))
  for (i in 1:length(par)) {

      for (j in 1:lags) {
        autocor[j, i] <- cor(chain[1:I(nrow(chain)-j), ..i], chain[(j+1):nrow(chain), ..i])
      }
  } 
  
  autocor %<>% as.data.table()
  setnames(autocor, c("V1", "V2", "V3"), par)
  return(autocor)

}

ac1 <- autocors(chain1)
ac2 <- autocors(chain2) 

autocorplot <- function(ac1, ac2) {
  
  par <- colnames(ac1)
  lags <- nrow(ac1)
  
  for (i in 1:length(par)) {
  
    ggplot() +
    geom_bar(ac1, mapping = aes_string(x = I(1:lags), y = par[i]), col = "#CC3399", width = 0.1, alpha = 0.4, stat = "identity") +
    geom_bar(ac2, mapping = aes_string(x = I(1:lags), y = par[i]), col = "#FF99FF", width = 0.1, alpha = 0.4, stat = "identity") +
    ggtitle("Chain 1") +
    xlab("Lag") +
    ylab("Autocorrelation") +
    ylim(0,1) + # fix y axis to compare among different parameters
    theme(
      plot.title = element_text(color = "grey44", size = 11, face = "bold.italic", hjust = 0.45),
      axis.title.x = element_text(color = "#333333", size = 10, face = "bold"),
      axis.title.y = element_text(color = "#333333", size = 10, face = "bold"))
  
  ggsave(paste0("../Output/acplot_", par[i], ".png"), width = 8, height = 4)
  
  }
  
}

autocorplot(ac1, ac2)

