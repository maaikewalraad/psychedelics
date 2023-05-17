rm(list = ls())
require(data.table)
require(magrittr) # special pipe operators
require(psych) # describe()
require(dplyr) # data wrangling
require(knitr) 
require(ggpubr)

set.seed(7079540) # for reproducibility

# Import data
dt <- read.csv("Data/AHRI_DATASET_PM_MANUSCRIPT_DATA.csv")
dt %<>% as.data.table
md <- read.csv("Data/AHRI_DATASET_PM_MANUSCRIPT_CODEBOOK.csv") # metadata
md %<>% as.data.table

# Data wrangling ----------------------------------------------------------

dt %<>% select(PHQ9_SCORE, PM1_DIAG_CONDITION, C_DP) # keep used variables only
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

# linearity
plot(dt$PM1_DIAG_CONDITION, dt$PHQ9_SCORE)
abline(mod) # add fitted regression line to scatterplot

# homoscedasticity (constant residual variance)
plot(mod)[1] 

# normality in the outcome
dt %$% hist(PHQ9_SCORE) # DEPENDENT VARIABLE 
dt %>% filter(PM1_DIAG_CONDITION == 0) %$% hist(PHQ9_SCORE) 
dt %>% filter(PM1_DIAG_CONDITION == 1) %$% hist(PHQ9_SCORE) 

# normality in residuals
plot(mod)[2] # QQ-plot shows divergence from normality at the tails

# independence of observations (explain analytically)


# Gibbs sampling ----------------------------------------------------------

# First step: specify priors for each parameter. These are based on the historical data (see manuscript)
# b0 priors
mu00 <- 1 
zeta00 <- 1.0E4  

# b1 priors
mu10 <- 1
zeta10 <- 1.0E4
nu10 <- 100 # prior df's -> bigger = wider tails

# variance (σ^2) priors
a_0 = 1 
b_0 = 1
sig2_0 = 1/rgamma(1, shape = a_0, rate = b_0) 

# Second step: make a Gibbs function

post_b1 <- function(b0, b1, sig2, y, x1) { # MH function for b1
  I <- -(b1^2) * ((sum(x1^2))/(2*sig2)) + b1 * ((sum(x1*(y - b0))) / (sig2))
  II <- (1 + ((b1 - mu10)^2)/(nu10*zeta10))^((nu10+1)/2)
  r_post <- exp(I) * II
  return(r_post) 
}

gibbs.chains <- function(b0, b1, sig2, y, x1) {
  set.seed(7079540)
  warmup <- 1000
  H <- 25000 + warmup # number of samples drawn 
  n_par <- 3
  simulated <- matrix(0, nrow = H, ncol = n_par)
  colnames(simulated) <- c("b0","b1", "sig2") 
  simulated[1, ] <- c(b0, b1, sig2)
  
  for (h in 2:H) {
    # determine full conditional value for...
    # b0:
    mu0 <- (sum(y-b1*x1)/sig2 + (mu00/zeta00))/((nrow(dt) / sig2) + (1 / zeta00))  
    sd0 <- sqrt(1 / ((nrow(dt) / sig2) + (1 / zeta00)))
    b0 <- rnorm(1, mu0, sd0)
    
    # b1 (incl. MH, Non-conjugate prior):
    beta1_c <- b1 
    beta1_n <- rnorm(1, mean = 0, sd = 1)   # Sample from proposal density
    r_post <- post_b1(b0, beta1_n, sig2, y, x1)/post_b1(b0, beta1_c, sig2, y, x1) # Ratio of posterior
    r_prop <- beta1_c / beta1_n   # Ratio of proposal distribution
    r <- r_post * r_prop # Acceptance ratio r
    u <- runif(1, min = 0, max = 1)   # 'Sample a probability' from the uniform distribution
    ifelse(u <= abs(r), b1 <- beta1_n, b1 <- beta1_c) 
    
    # # b1 (exc. MH)
    # mu1 <- (sum(x1*(y-b0)) / sig2 + (mu10/zeta10)) / ((sum(x1^2) / sig2) + (1 / zeta10))
    # sd1 <- sqrt(1 / ((sum(x1^2) / sig2) + (1 / zeta10)))
    # b1 <- rnorm(1, mu1, sd1)
    
    # var: 
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
chain1 <- gibbs.chains(b0 = 0.2, b1 = 0.3, sig2 = 0.5, y = dt$PHQ9_SCORE, x1 = dt$PM1_DIAG_CONDITION)
chain2 <- gibbs.chains(b0 = 2, b1 = 3, sig2 = 5, y = dt$PHQ9_SCORE, x1 = dt$PM1_DIAG_CONDITION) 


# Assess convergence ------------------------------------------------------

# Subset chains
chain1 %<>% as.data.table()
chain2 %<>% as.data.table()
parameters <- colnames(chain1)  # saving for plot

# b0:
chains_b0 <- cbind(chain1[1001:nrow(chain1), 1], chain2[1001:nrow(chain2), 1])  # extract b0 & remove the warm up 
colnames(chains_b0) <- c("chain1", "chain2") 
# b1:
chains_b1 <- cbind(chain1[1001:nrow(chain1), 2], chain2[1001:nrow(chain2), 2])   # extract b1 & remove the warm up 
colnames(chains_b1) <- c("chain1", "chain2")
# var:
chains_var <- cbind(chain1[1001:nrow(chain1), 3], chain2[1001:nrow(chain2), 3])   # extract variance & remove the warm up 
colnames(chains_var) <- c("chain1", "chain2")

# Traceplots
traceplot <- function(chains) {   # only handles 2 chains (per parameter)
  
  x <- 0
  colnames <- NULL;
  while (x < ncol(chains))  {  
    x <- x+1  
    colnames <- c(colnames, paste0("chain", x, collapse = "," ))
  }
  
  traceplot <- chains %>% 
    ggplot(aes(x = 1:nrow(chains))) + 
    geom_line(aes(y = chain1), color = "#CC3399", linetype = "dashed", size = 0.4, alpha = 0.8) + 
    geom_line(aes(y = chain2), color= "#339966", linetype = "twodash", size = 0.4, alpha = 0.6) +
    ggtitle("Convergence plot") +
    ylab("Sampled value") +
    xlab("Number of iterations") +
    theme(
      plot.title = element_text(color = "black", size = 11, face = "bold.italic", hjust = 0.45),
      axis.title.x = element_text(color = "#333333", size = 10, face = "bold"),
      axis.title.y = element_text(color = "#333333", size = 10, face = "bold"))
  
  return(traceplot)
  
}

traceplots <- ggarrange(traceplot(chains_b0), traceplot(chains_b1), traceplot(chains_var),
                        labels = parameters,
                        ncol = 2, nrow = 2)

annotate_figure(traceplots,
                top = text_grob("Traceplots for the parameters (2 chains each) \n", color = "black", size = 18, face = "bold.italic"))

#### Alternatively
mcmcplot(mcmcout = chains_b0)

# Autocorrelation plots
autocorplot <- function(chains) { 
  # 'chains' object should be a dataframe with samples in rows and different chains in separate columnns. Handles only 2 chains.
  lags <- 40 # number of lags
  x <- 0
  colnames <- NULL;
  while (x < ncol(chains))  {  
    x <- x+1  
    colnames <- c(colnames, paste0("chain", x, collapse = "," ))
  }
  autocor <- matrix(0, nrow = 40, ncol = ncol(chains))
  
  for (chain in 1:ncol(chains)) {
    
    for (i in 1:lags) {
      autocor[i, chain] <- cor(chains[1:I(nrow(chains)-i), chain], chains[(i+1):nrow(chains), chain])
    }
    
    autocor <- as.data.frame(autocor)
    colnames(autocor) <- colnames
    
  }
  plot1 <- autocor %>%
    ggplot(aes(x = 1:lags)) +
    geom_bar(aes(y = chain1), col = "#CC3399", fill = "#FFFFFF", alpha = 0.4, stat = "identity") +
    ggtitle("Chain 1") +
    xlab("Lag") +
    ylab("Autocorrelation") +
    theme(
      plot.title = element_text(color = "grey44", size = 11, face = "bold.italic", hjust = 0.45),
      axis.title.x = element_text(color = "#333333", size = 10, face = "bold"),
      axis.title.y = element_text(color = "#333333", size = 10, face = "bold"))
  
  plot2 <- autocor %>% 
    ggplot(aes(x = 1:lags)) +
    geom_bar(aes(y = chain2), col =  "#339966", fill = "#FFFFFF", alpha = 0.9, stat = "identity") +
    ggtitle("Chain 2") +
    xlab("Lag") + 
    ylab("Autocorrelation") +
    theme(
      plot.title = element_text(color = "grey44", size = 11, face = "bold.italic", hjust = 0.45),
      axis.title.x = element_text(color = "#333333", size = 10, face = "bold"),
      axis.title.y = element_text(color = "#333333", size = 10, face = "bold"))
  
  corplots <- ggarrange(plot1, plot2)
  
  return(corplots)
  
}

autocorplots <- ggarrange(autocorplot(chains_b0), autocorplot(chains_b1), autocorplot(chains_var),
                          labels = parameters,
                          ncol = 2, nrow = 2) 

annotate_figure(autocorplots,
                top = text_grob("Autocorrelation plots for the parameters (2 chains each) \n", color = "black", size = 18, face = "bold.italic"))


# Obtain parameter estimates, credible intervals and interpretation -------

# Parameter estimates (first chain):
perc <- c(0.025, 0.975)
bhat0_1 <- c(mean(chains_b0[, 1]), sd(chains_b0[, 1]), quantile(chains_b0[, 1], perc))
bhat1_1 <- c(mean(chains_b1[, 1]), sd(chains_b1[, 1]), quantile(chains_b1[, 1], perc))
estimates_c1 <- rbind(bhat0_1, bhat1_1) # also add variance mean
rownames(estimates_c1) <- c("Intercept", "PM")
colnames(estimates_c1) <- c("Beta", "SD", "2.5% ", "97.5%")
estimates_c1

# Parameter estimates (second chain):
bhat0_2 <- c(mean(chains_b0[, 2]), sd(chains_b0[, 2]), quantile(chains_b0[, 2], perc))
bhat1_2 <- c(mean(chains_b1[, 2]), sd(chains_b1[, 2]), quantile(chains_b1[, 2], perc))
estimates_c2 <- rbind(bhat0_2, bhat1_2)
rownames(estimates_c2) <- c("Intercept", "PM")
colnames(estimates_c2) <- c("Beta", "SD", "2.5% ", "97.5%")
estimates_c2

