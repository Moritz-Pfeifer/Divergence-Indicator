tex_table <- TRUE
if(tex_table) {
  library(dplyr)
  library(knitr)
  library(kableExtra)
}
library(ggplot2)
#############################################
### Install and load 'rbsvar' from Github ###
#############################################

if(!("rbsvar" %in% installed.packages())) {
  if(!("devtools" %in% installed.packages())) {
    install.packages("devtools")
  }
  devtools::install_github("jetroant/rbsvar")
  library(rbsvar)
} else {
  library(rbsvar)
}

working_directory <- dirname(rstudioapi::documentPath())
setwd(working_directory)

#################
### Load Data ###
#################

if(!("rmatio" %in% installed.packages())) install.packages("rmatio")

# Read the .mat file

file_path <- "Divergence_Data_M.mat"
divergence <- rmatio::read.mat(file_path)

##########################
### Plot data ###
##########################

y <- divergence$data
colnames(y) <- unlist(divergence$varNames)
y_ts <- ts(y, start = c(1999,1), frequency = 12)
plot(y_ts)

##########################
### Estimate the model ###
##########################

# init_rbsvar() initializes the model and chooses efficient initial values for
# when sampling from the posterior distribution of the parameters in the next step.
svar_divergence <- init_rbsvar(y = y,
                          lags = 12,
                          constant = T,
                          type = "svar",
                          shrinkage = Inf
)

# est_rbsvar() estimates the model by drawing n(=2) chains from the posterior
# distribution of the parameters, starting from values provided by init_rbsvar().
# Chains are of length K*N (=10*N) and they are thinned by a factor of K(=10),
# thus the final sample consists of n*N (=2*N) non-independent draws. The chains
# are drawn with differential evolution Monte Carlo algorithm.
#
# Even though the algorithm is parallelized, this might take some time.
N <- 1000000 # 1000000
burn <- 850000 # OK until 600000
output_divergence<- est_rbsvar(svar_divergence,
                           N = N,
                           n = 2,
                           K = 10,
                           parallel_chains = TRUE,
                           parallel_likelihood = TRUE)

# Save the posterior sample (Or load by commenting/uncommenting the below lines)
saveRDS(output_divergence, "output_divergence.rds")
output_divergence <- readRDS("output_divergence.rds")

# convergence() computes the R_hat convergence statistic of Gelman et al. (2013)
# with respect to every parameter in the model and plots the results.
# Values below the threshold of 1.1 can be interpreted as implying convergence of the chains.
convergence_divergence <- convergence(output_divergence, burn=burn)

# Plot for the distribution of marginal R_hat
ggplot(convergence_divergence, aes(x = R_hat)) +
  geom_histogram(bins = 40, fill = "grey", color = "black") +
  labs(title = NULL, x = "R_hat", y = "Frequency") +
  theme_minimal() +
  theme(panel.grid = element_blank()) + # Remove the grid
  geom_vline(xintercept = 1.1, linetype = "dashed", color = "black", size = 1) # Add vertical line at 1.1

# Density trace plot
plot_densities(output_divergence)

# "Least converged" parameters
plot_chains(output_divergence, variable = rownames(convergence_divergence[order(convergence_divergence[,1], decreasing = T),])[1], burn = burn)
plot_chains(output_divergence, variable = rownames(convergence_divergence[order(convergence_divergence[,1], decreasing = T),])[2], burn = burn)
plot_chains(output_divergence, variable = rownames(convergence_divergence[order(convergence_divergence[,1], decreasing = T),])[3], burn = burn)


# Some other parameters
plot_chains(output_divergence,
            "a60",
            breaks = 30,
            burn = burn)

plot_chains(output_divergence,
            "b9",
            breaks = 30,
            burn = burn)

plot_chains(output_divergence,
            "sgt6",
            breaks = 30,
            burn = burn)

plot_chains(output_divergence,
            "a26",
            breaks = 30,
            burn = burn)

####################################
### Tail exponents of the shocks ### (Highly dependent on convergence)
####################################

post <- post_sample(output_divergence, burn)

# Sgt posterior means
sgt_postmean <- matrix(apply(post$s[,grep("sgt", colnames(post$s))], 2, mean), ncol = 3)
colnames(sgt_postmean) <- c("skew", "log_p", "log_q")
sgt_postmean

# Tail exponent marginal posteriors
par(mfrow = c(1,3))
post <- post_sample(output_divergence, burn)
te_mat <- matrix(NA, ncol = ncol(y), nrow(post$s))
for(i in 1:ncol(y)) {
  log_p <- post$s[,which(colnames(post$s) == "sgt3") + i]
  log_q <- post$s[,which(colnames(post$s) == "sgt6") + i]
  tail_exponent <- exp(log_p + log_q)
  hist(tail_exponent, breaks = 25, probability = T, col = "peachpuff", xlab = "Tail exponent",
       main = paste0("Shock ", i, " / Pr[variance exists] = ", round(100*mean(tail_exponent > 2),1), "%"))
}
par(mfrow = c(1,1))
# dev.off()

##############################
### Skewness of the shocks ### (Skew-parameter (support from -1 to 1), NOT skewness in moment sense)
##############################

post_skew <- post$s[,grep("sgt", colnames(post$s))][,1:3]
par(mfrow = c(1,3))
for(i in 1:3) {
  hist(post_skew[,i], breaks = 20, probability = T, col = "peachpuff", xlab = "lambda",
       main = paste0("Shock ", i))
}
par(mfrow = c(1,1))

##################################
### Impulse response functions ### (These take some time to compute due to lack of parallelization)
##################################


# Compute impulse responses of all the shocks to all variables
irf_obj_all <- irf(model = svar_divergence,
                   output = output_divergence,
                   burn = burn,
                   N = 10000,
                   horizon = 60,
                   cumulate = c())

# Compute impulse responses of the supposed monetary policy (balance sheet) shock to all variables
irf_obj_mon_bs <- irf(model = svar_divergence,
                   output = output_divergence,
                   burn = burn,
#                  N = 10000,
                   N = NA,  # Changed N to NA as requested by authors
                   horizon = 60,
                   cumulate = c(),
                   shocks = 1)

# Compute impulse responses of the supposed monetary policy (shadow rate) shock to all variables
irf_obj_mon_r <- irf(model = svar_divergence,
                   output = output_divergence,
                   burn = burn,
                   N = 10000,
                   horizon = 60,
                   cumulate = c(),
                   shocks = 2)

# Compute impulse responses of the supposed asymmetric shock to all variables
irf_obj_asym <- irf(model = svar_divergence,
                   output = output_divergence,
                   burn = burn,
                   N = 10000,
                   horizon = 60,
                   cumulate = c(),
                   shocks = 3)

# Save the impulse response objects (Or load)
saveRDS(irf_obj_all, "/Users/moritzpfeifer/Desktop/Cycles/replication_files_all/empirical_application_divergence/irf/irf_obj_all.rds")
saveRDS(irf_obj_mon_bs, "/Users/moritzpfeifer/Desktop/Cycles/replication_files_all/empirical_application_divergence/irf/irf_obj_mon_bs.rds")
saveRDS(irf_obj_mon_r, "/Users/moritzpfeifer/Desktop/Cycles/replication_files_all/empirical_application_divergence/irf/irf_obj_mon_r.rds")
saveRDS(irf_obj_asym, "/Users/moritzpfeifer/Desktop/Cycles/replication_files_all/empirical_application_divergence/irf/irf_obj_asym.rds")
irf_obj_all <- readRDS("/Users/moritzpfeifer/Desktop/Cycles/replication_files_all/empirical_application_divergence/irf/irf_obj_all.rds")
irf_obj_mon_bs <- readRDS("/Users/moritzpfeifer/Desktop/Cycles/replication_files_all/empirical_application_divergence/irf/irf_obj_mon_bs.rds")
irf_obj_mon_r <- readRDS("/Users/moritzpfeifer/Desktop/Cycles/replication_files_all/empirical_application_divergence/irf/irf_obj_mon_r.rds")
irf_obj_asym <- readRDS("/Users/moritzpfeifer/Desktop/Cycles/replication_files_all/empirical_application_divergence/irf/irf_obj_asym.rds")


# Plot the impulse response functions and optionally save them as pdf-files.
#     The Shadow Rate rate is plotted in basis points, whereas the other variables
#     are logarithmic time-series

irf_plot(irf_obj_all,
         varnames = c("Money", "Rate", "Divergece"),
         leg = FALSE,
         color="#156082",
         mar = c(2,2,2,1),
         normalize = c(rep(100, 3), 1))

irf_plot(irf_obj_mon_bs,
         varnames = c("Money", "Rate", "Divergece"),
         leg = FALSE,
         color="#156082",
         mar = c(2,2,2,1),
         mfrow = c(2,3),
         normalize = c(rep(100, 3), 1))

irf_plot(irf_obj_mon_r,
         varnames = c("Money", "Rate", "Divergece"),
         leg = FALSE,
         color="#156082",
         mar = c(2,2,2,1),
         mfrow = c(2,3),
         normalize = c(rep(100, 3), 1))

irf_plot(irf_obj_asym,
         varnames = c("Money", "Rate", "Divergece"),
         leg = FALSE,
         color="#156082",
         mar = c(2,2,2,1),
         mfrow = c(2,3),
         normalize = c(rep(100, 3), 1))

### A positive (expansionary) UMP-shock conditional on one sign restriction
### (divergence being negative in first 12 months) ###

f <- function(x) if(sum(x[0:12] > 0) == 0) return(1) else return(0)
t <- apply(irf_obj_mon_bs[[1]][3,,], 2, f)
irf_obj_mon_cond <- irf_obj_mon_bs
irf_obj_mon_cond[[1]] <- irf_obj_mon_bs[[1]][,,which(t == 1)]
t1 <- t  # Save t for decomp object

irf_plot(irf_obj_mon_cond,
         varnames = c("Money", "Rate", "Divergece"),
         leg = FALSE,
         color="#156082",
         mar = c(2,2,2,1),
         mfrow = c(2,3),
         normalize = c(rep(100, 3), 1))

### A positive (contractionary) monetary policy shock conditional on one sign restriction
### Rate being positive in the first 12 months ###

f <- function(x) if(sum(x[0:12] > 0) == 0) return(1) else return(0)
t <- apply(irf_obj_mon_r[[2]][2,,], 2, f)
irf_obj_mon_cond <- irf_obj_mon_r
irf_obj_mon_cond[[2]] <- irf_obj_mon_r[[2]][,,which(t == 0)]

irf_plot(irf_obj_mon_cond,
         varnames = c("Money", "Rate", "Divergece"),
         leg = FALSE,
         mar = c(2,2,2,1),
         color="#156082",
         mfrow = c(2,3),
         normalize = c(rep(100, 3), 1))

### An asymmetric-shock conditional on one sign restriction ###
### Interest Rate being negative in the first 24 month ###

f <- function(x) if(sum(x[0:24] > 0) == 0) return(1) else return(0)
t <- apply(irf_obj_asym[[3]][2,,], 2, f)
irf_obj_mon_cond <- irf_obj_asym
irf_obj_mon_cond[[3]] <- irf_obj_asym[[3]][,,which(t == 1)]

irf_plot(irf_obj_mon_cond,
         varnames = c("Money", "Rate", "Divergece"),
         leg = FALSE,
         mar = c(2,2,2,1),
         color="#156082",
         mfrow = c(2,3),
         normalize = c(rep(100, 3), 1))



#######################################
### Historical shock decompositions ###
#######################################

# Compute historical shock decompositions
decomp_obj <- rbsvar:::shock_decomp(model = svar_divergence,
                           output = output_divergence,
                           burn = burn,
                           N = 10000)

#-------------Balance Sheet--------------------#
# October 2012 (UMP) shock decomposition (ESM)
decomp_oct_2012 <- rbsvar:::shock_decomp_pick(decomp_obj = decomp_obj,
                                             show_date = c(2012,10),
                                             start_date = c(1999,1), # 12 lags
                                             freq = 12)

# March 2015 (UMP) shock decomposition (EAPP)
decomp_mar_2015 <- rbsvar:::shock_decomp_pick(decomp_obj = decomp_obj,
                                             show_date = c(2015,3),
                                             start_date = c(1999,1), # 12 lags
                                             freq = 12)

# March 2020 (UMP) shock decomposition (PEPP)
decomp_mar_2020 <- rbsvar:::shock_decomp_pick(decomp_obj = decomp_obj,
                                             show_date = c(2020,3),
                                             start_date = c(1999,1), # 12 lags
                                             freq = 12)

#-------------Interest Sheet--------------------#
# October 2008 (rate) shock decomposition (neg/expansioanry)
decomp_oct_2008 <- rbsvar:::shock_decomp_pick(decomp_obj = decomp_obj,
                                              show_date = c(2008,10),
                                              start_date = c(1999,1),
                                              freq = 12)

# September 2022 (rate) shock decomposition (pos/contractionary)
decomp_sep_2022 <- rbsvar:::shock_decomp_pick(decomp_obj = decomp_obj,
                                              show_date = c(2022,9),
                                              start_date = c(1999,1),
                                              freq = 12)


#-------------Divergence--------------------#
# Jan 2000 (symmetric) shock decomposition (QE)
decomp_jan_2000 <- rbsvar:::shock_decomp_pick(decomp_obj = decomp_obj,
                                              show_date = c(2000,1),
                                              start_date = c(1999,1), # 12 lags
                                              freq = 12)

# October 2009 (asymmetric) shock decomposition (Euro debt crisis)
decomp_oct_2009 <- rbsvar:::shock_decomp_pick(decomp_obj = decomp_obj,
                                              show_date = c(2009,10),
                                              start_date = c(1999,1), # 12 lags
                                              freq = 12)

# July 2021 (asymmetric) shock decomposition (Inflation/end of Covid-19)
decomp_jul_2021 <- rbsvar:::shock_decomp_pick(decomp_obj = decomp_obj,
                                              show_date = c(2021,7),
                                              start_date = c(1999,1), # 12 lags
                                              freq = 12)



#--------------------Plot Marginal Poster Distributions------------------------#

# Positive UMP monetary policy shocks (ranges)
par(mfrow = c(1,3))
par(mar = c(2,2,2,1))

# Positive UMP monetary policy shocks (points)
hist(decomp_oct_2012$E[,1], main = "October 2012 (ESM)",
     col = "darkolivegreen3", probability = T, breaks = 20)
abline(v = 0, lty = 2, lwd = 2)

hist(decomp_mar_2015$E[,1], main = "March 2015 (EAPP)",
     col = "darkolivegreen3", probability = T, breaks = 20)
abline(v = 0, lty = 2, lwd = 2)

hist(decomp_mar_2020$E[,1], main = "March 2020 (PEPP)",
     col = "darkolivegreen3", probability = T, breaks = 20)
abline(v = 0, lty = 2, lwd = 2)


# Expansionary/contractionary monetary policy shocks:
hist(decomp_oct_2008$E[,2], main = "October 2008",
     col = "darkolivegreen3", probability = T, breaks = 10)
abline(v = 0, lty = 2, lwd = 2)
hist(decomp_sep_2022$E[,2], main = "September 2022",
     col = "salmon", probability = T)
abline(v = 0, lty = 2, lwd = 2)


# Asymmetric/symmetric shocks point
par(mfrow = c(1,3))
par(mar = c(2,2,2,1))
hist(decomp_jan_2000$E[,3], main = "January 2000",
     col = "salmon", probability = T, breaks = 20)
abline(v = 0, lty = 2, lwd = 2)

hist(decomp_oct_2009$E[,3], main = "October 2009",
     col = "darkolivegreen3", probability = T, breaks = 20)
abline(v = 0, lty = 2, lwd = 2)

hist(decomp_jul_2021$E[,3], main = "July 2021",
     col = "darkolivegreen3", probability = T, breaks = 20)
abline(v = 0, lty = 2, lwd = 2)


#############################################
### Narrative signs - shock probabilities ###
#############################################

narrative_sign_probs <- function(decomp_obj,
                                 start_date,
                                 freq,
                                 dates,
                                 signs) {

  N <- dim(decomp_obj$E)[1]
  m <- dim(decomp_obj$E)[3]
  probs <- matrix(NA, nrow = length(dates), ncol = m)
  total_probs <- rep(NA, m)
  for(i in 1:m) {
    for(j in 1:length(dates)) {
      picked <- rbsvar:::shock_decomp_pick(decomp_obj, show_date = dates[[j]], start_date = start_date, freq = freq)
      if(signs[j] == 1) agree_dummy <- picked$E[,i] > 0
      if(signs[j] == -1) agree_dummy <- picked$E[,i] < 0
      if(j == 1) {
        agree_dummies <- agree_dummy
      } else {
        agree_dummies <- cbind(agree_dummies, agree_dummy)
      }
    }
    probs[,i] <- apply(agree_dummies, 2, mean)
    total_probs[i] <- mean(apply(agree_dummies, 1, function(x) ifelse(sum(x) == length(x), TRUE, FALSE)))
  }
  ret <- rbind(probs, total_probs)
  rownames(ret) <- c(1:length(dates), "Total")
  ret
}

date_list <- list(c(2000,1),
                  c(2009,10),
                  c(2021,7))

date_list <- list(c(2012,10),
                  c(2015,3),
                  c(2020,3))

date_list <- list(c(2008,10),
                  c(2022,9))

probs <- narrative_sign_probs(decomp_obj = decomp_obj,
                              start_date = c(1999,1),
                              freq = 12,
                              dates = date_list,
                              signs = c(1,1, 1) #rep(1, length(date_list))
                              )
for(i in 1:length(date_list)) {
  rownames(probs)[i] <- paste(date_list[[i]][1], date_list[[i]][2], sep = "/")
}
colnames(probs) <- paste("Shock", 1:ncol(y))

# Fraction of the posterior satisfying the narrative sign restrictions for every shock
round(probs, 2)
probs
saveRDS(probs, "narrative_probs.rds")

# Tex table of the above
if(tex_table) {
  rownames(probs)[nrow(probs)] <- " "
  knitr::kable(round(probs, 2),
               "latex", booktabs = T, align = "r",
               caption = "Fraction of the posterior satisfying the narrative sign restrictions for every shock") %>%
    kable_styling(latex_options = c("hold_position")) %>%
    kableExtra::pack_rows(index = c("Contractionary shocks" = 4, "Expansionary shocks" = 4, "Total" = 1)) %>%
    kable_styling(position = "center") %>%
    kable_styling(font_size = 12)
}


########################################################
### Narrative signs - shock probabilities - two-ways ###
########################################################

shock_decomp2 <- function(model,
                          output,
                          burn = 0,
                          N = 1000,
                          verbose = TRUE) {

  if(model$type != "svar") stop("model$type must be 'svar'")
  m <- ncol(model$y)
  p <- model$lags
  yy <- model$xy$yy
  xx <- model$xy$xx
  n <- output$args$n
  m0 <- output$args$m0
  chains <- output$chains
  b_indices <- (model$cpp_args$first_b + 1):model$cpp_args$first_sgt
  a_indices <- 1:model$cpp_args$first_b

  #Collect A and B sample
  burn <- burn + 1
  b_list <- list()
  a_list <- list()
  for(i in 1:n) {

    #Pick one chain
    row_indices <- seq(from = m0 + i - n, by = n, length.out = output$args$N + 1)
    one_chain <- chains[row_indices,][-c(1:burn),]

    #Pick A
    a_list[[i]] <- one_chain[,a_indices]

    #Pick B
    b_list[[i]] <- one_chain[,b_indices]
  }
  for(i in 1:n) {
    if(i == 1) {
      full_b <- b_list[[i]]
      full_a <- a_list[[i]]
    } else {
      full_b <- rbind(full_b, b_list[[i]])
      full_a <- rbind(full_a, a_list[[i]])
    }
  }

  #Sample N draws and compute B_inv
  if(is.na(N)) {
    N <- nrow(full_a)
    rows <- 1:N
  } else {
    rows <- sample.int(nrow(full_a), N, replace = TRUE)
  }

  ret_E <- array(NA, dim = c(N, nrow(yy), m))
  ret_U <- list()
  for(i in 1:m) ret_U[[i]] <- array(NA, dim = c(N, nrow(yy), m))

  full_a <- full_a[rows,]
  full_b <- full_b[rows,]
  full_binv <- matrix(NA, ncol = m^2, nrow = nrow(full_b))
  for(j in 1:nrow(full_b)) {
    B <- diag(m)
    B[,] <- full_b[j,]
    full_binv[j,] <- c(solve(B))
  }

  if(verbose == TRUE) print(paste0("Computing shock decompositions..."))
  if(verbose == TRUE) pb <- txtProgressBar(min = 0, max = N, style = 3)
  for(i in 1:N) {

    B <- diag(m)
    B_inv <- diag(m)
    B[,] <- full_b[i,]
    B_inv[,] <- full_binv[i,]
    A <- matrix(full_a[i,], ncol = m)
    U <- yy - xx %*% A
    E <- U %*% t(B)
    ret_E[i,,] <- E
    for(j in 1:m) {
      ret_U[[j]][i,,] <- matrix(c(t(E)) * B_inv[j,], byrow = TRUE, ncol = m)
    }
    if(verbose == TRUE) setTxtProgressBar(pb, i)
  }
  if(verbose == TRUE) close(pb)

  ret <- list("E" = ret_E,
              "U" = ret_U,
              "par" = list("a" = full_a, "b" = full_b, "b_inv" = full_binv))
  ret
}


decomp_obj_2 <- shock_decomp2(model = svar_divergence,
                              output = output_divergence,
                              burn = burn,
#                              N = 10000
                               N = NA) # Set N = NA as requested by authors


probs2 <- narrative_sign_probs(decomp_obj = decomp_obj_2,
                               start_date = c(1999,1),
                               freq = 12,
                               dates = date_list,
                               signs = c(1,1,1) #rep(1, length(date_list))
)
for(i in 1:length(date_list)) {
  rownames(probs2)[i] <- paste(date_list[[i]][1], date_list[[i]][2], sep = "/")
}
colnames(probs2) <- paste("Shock", 1:ncol(y))

# Fraction of the posterior satisfying the narrative sign restrictions for every shock
round(probs2, 2)
probs2

# Compute shock decompositions conditional on negative divergence
### (divergence being negative in first 12 months) ###

decomp_obj_2_res <- decomp_obj_2
decomp_obj_2_res$E <- decomp_obj_2_res$E[which(t1 == 1),,]

probs2_res <- narrative_sign_probs(decomp_obj = decomp_obj_2_res,
                               start_date = c(1999,1),
                               freq = 12,
                               dates = date_list,
                               signs = c(1,1,1) #rep(1, length(date_list))
)
for(i in 1:length(date_list)) {
  rownames(probs2_res)[i] <- paste(date_list[[i]][1], date_list[[i]][2], sep = "/")
}
colnames(probs2_res) <- paste("Shock", 1:ncol(y))

# Fraction of the posterior satisfying the narrative sign restrictions for negative divergence shock
round(probs2_res, 2)
probs2_res

########################
### Ranges of Shocks ###
########################

# Let's create a new function for ranges
shock_decomp_pick_new <- function (decomp_obj, start_show_date, end_show_date, start_date, freq)
{
  row_indices <- ts(1:dim(decomp_obj$E)[2], start = start_date, frequency = freq)
  row_index <- window(row_indices, start = start_show_date, end = end_show_date)
  E <- decomp_obj$E[, row_index, ]
  m <- dim(decomp_obj$U[[1]])[3]
  U <- list()
  for (i in 1:m) {
    U[[i]] <- decomp_obj$U[[i]][, row_index, ]
  }
  ret <- list(E = E, U = U)
  ret
}

#------------------------------Compute Ranges---------------------------------#

#-------------Balance Sheet--------------------#
# October 2012 (UMP) shock decomposition (ESM)
ESM <- shock_decomp_pick_new(decomp_obj = decomp_obj,
                             start_show_date = c(2012,9),
                             end_show_date = c(2012,11),
                             start_date = c(1999,1), # 12 lags
                             freq = 12)


# March 2015 (UMP) shock decomposition (EAPP)
EAPP <- shock_decomp_pick_new(decomp_obj = decomp_obj,
                              start_show_date = c(2015,3),
                              end_show_date = c(2015,5),
                              start_date = c(1999,1), # 12 lags
                              freq = 12)

# March 2020 (UMP) shock decomposition (PEPP)
PEPP <- shock_decomp_pick_new(decomp_obj = decomp_obj,
                              start_show_date = c(2012,9),
                              end_show_date = c(2012,11),
                              start_date = c(1999,1), # 12 lags
                              freq = 12)


#-------------Divergence--------------------#
# EURO launch
EURO <- shock_decomp_pick_new(decomp_obj = decomp_obj,
                              start_show_date = c(1999,10),
                              end_show_date = c(2000,3),
                              start_date = c(1999,1), # 12 lags
                              freq = 12)


# Sovereign Debt Crisis
SDC <- shock_decomp_pick_new(decomp_obj = decomp_obj,
                             start_show_date = c(2009,10),
                             end_show_date = c(2010,3),
                             start_date = c(1999,1), # 12 lags
                             freq = 12)

# Post Covid-19 inflation
Covid19 <- shock_decomp_pick_new(decomp_obj = decomp_obj,
                                 start_show_date = c(2021,7),
                                 end_show_date = c(2021,12),
                                 start_date = c(1999,1), # 12 lags
                                 freq = 12)


# UMP Shock ranges
hist(apply(ESM$E[, , 1], 1, mean), main = "ESM 2012",
     col = "darkolivegreen3", probability = T, breaks = 20)
abline(v = 0, lty = 2, lwd = 2)

hist(apply(EAPP$E[, , 1], 1, mean), main = "EAPP 2015",
     col = "darkolivegreen3", probability = T, breaks = 20)
abline(v = 0, lty = 2, lwd = 2)

hist(apply(PEPP$E[, , 1], 1, mean), main = "ESM 2020",
     col = "darkolivegreen3", probability = T, breaks = 20)
abline(v = 0, lty = 2, lwd = 2)


# Asymmetric/symmetric shocks ranges
hist(apply(EURO$E[, , 3], 1, mean), main = "Euro 1999/2000",
     col = "salmon", probability = T, breaks = 20)
abline(v = 0, lty = 2, lwd = 2)

hist(apply(SDC$E[, , 3], 1, mean), main = "EAPP 2015",
     col = "darkolivegreen3", probability = T, breaks = 20)
abline(v = 0, lty = 2, lwd = 2)

hist(apply(Covid19$E[, , 3], 1, mean), main = "Post-Covid Inflation",
     col = "darkolivegreen3", probability = T, breaks = 20)
abline(v = 0, lty = 2, lwd = 2)


# Let's create ranges
start_year <- 2012
end_year <- 2013
date_list <- list()

# Loop through each year and month
for (year in start_year:end_year) {
  for (month in 1:12) {
    # Append each year and month as a vector to the list
    date_list[[length(date_list) + 1]] <- c(year, month)
  }
}
date_list

#
# Shock decompositions conditional upon negative divergence shock restrictions
#

#-------------Balance Sheet--------------------#
# October 2012 (UMP) shock decomposition (ESM)
decomp_oct_2012_res <- rbsvar:::shock_decomp_pick(decomp_obj = decomp_obj_2_res,
                                              show_date = c(2012,10),
                                              start_date = c(1999,1), # 12 lags
                                              freq = 12)

# March 2015 (UMP) shock decomposition (EAPP)
decomp_mar_2015_res <- rbsvar:::shock_decomp_pick(decomp_obj = decomp_obj_2_res,
                                              show_date = c(2015,3),
                                              start_date = c(1999,1), # 12 lags
                                              freq = 12)

# March 2020 (UMP) shock decomposition (PEPP)
decomp_mar_2020_res <- rbsvar:::shock_decomp_pick(decomp_obj = decomp_obj_2_res,
                                              show_date = c(2020,3),
                                              start_date = c(1999,1), # 12 lags
                                              freq = 12)


par(mfrow = c(1,3))
par(mar = c(2,2,2,1))

# Positive UMP monetary policy shocks (points)
hist(decomp_oct_2012_res$E[,1], main = "October 2012 (ESM)",
     col = "darkolivegreen3", probability = T, breaks = 20)
abline(v = 0, lty = 2, lwd = 2)

hist(decomp_mar_2015_res$E[,1], main = "March 2015 (EAPP)",
     col = "darkolivegreen3", probability = T, breaks = 20)
abline(v = 0, lty = 2, lwd = 2)

hist(decomp_mar_2020_res$E[,1], main = "March 2020 (PEPP)",
     col = "darkolivegreen3", probability = T, breaks = 20)
abline(v = 0, lty = 2, lwd = 2)



#-------------Balance Sheet--------------------#
# October 2012 (UMP) shock decomposition (ESM)
ESM_res <- shock_decomp_pick_new(decomp_obj = decomp_obj_2_res,
                             start_show_date = c(2012,9),
                             end_show_date = c(2012,11),
                             start_date = c(1999,1), # 12 lags
                             freq = 12)


# March 2015 (UMP) shock decomposition (EAPP)
EAPP_res <- shock_decomp_pick_new(decomp_obj = decomp_obj_2_res,
                              start_show_date = c(2015,3),
                              end_show_date = c(2015,5),
                              start_date = c(1999,1), # 12 lags
                              freq = 12)

# March 2020 (UMP) shock decomposition (PEPP)
PEPP_res <- shock_decomp_pick_new(decomp_obj = decomp_obj_2_res,
                              start_show_date = c(2012,9),
                              end_show_date = c(2012,11),
                              start_date = c(1999,1), # 12 lags
                              freq = 12)

# UMP Shock ranges
hist(apply(ESM_res$E[, , 1], 1, mean), main = "ESM 2012",
     col = "darkolivegreen3", probability = T, breaks = 20)
abline(v = 0, lty = 2, lwd = 2)

hist(apply(EAPP_res$E[, , 1], 1, mean), main = "EAPP 2015",
     col = "darkolivegreen3", probability = T, breaks = 20)
abline(v = 0, lty = 2, lwd = 2)

hist(apply(PEPP_res$E[, , 1], 1, mean), main = "ESM 2020",
     col = "darkolivegreen3", probability = T, breaks = 20)
abline(v = 0, lty = 2, lwd = 2)




###############################
### Shock volatility figure ###
###############################

mp_shock <- decomp_obj$E[,,1] # 1 = monetary policy shock; 3 = asymmetric shock
#mp_quant <- t(apply(mp_shock, 2, function(x) quantile(x, probs = c(0.05, 0.5, 0.95))))
#dim(mp_shock)

par(mar = c(3, 3, 2, 1))
col <- "tomato"
probs <- c(0.76, 0.00, 0.01)
mean_volatility <- ts(apply(mp_shock, 2, mean), start = attributes(y_ts)$tsp[1] + 1, frequency = attributes(y_ts)$tsp[3])
#y_labels <- as.expression(sapply(1:ncol(mean_volatility), function(x) bquote(sigma[.(x)])))
sub_shock <- apply(mp_shock, 2, function(x) quantile(x, probs = probs))
plot(mean_volatility, main = "", ylab = "", lwd = 1,
     ylim = c(-20, 20)
     #ylim = c(min(sub_shock) - 0.1, max(sub_shock) + 0.1),
     )#xaxt = "n")
#axis(1, at = seq(1960, 2010, by = 10), labels = rep("", 6))
grid()
fanplot::fan(data = sub_shock, data.type = "values", probs = probs,
             start = attributes(y_ts)$tsp[1] + 1, frequency = attributes(y_ts)$tsp[3],
             fan.col = colorRampPalette(c(col, "white")),
             rlab = NULL, ln = NULL)
#lines(mean_volatility, lwd = 1, lty = 2)
# dev.off()

if((i-1) %% rows == 0) {


  if(rows == 1) {
    plot(mean_volatility, main = "Shock Volatility", ylab = y_labels[i], lwd = 1, ylim = c(0, max(sub_volatility) + 0.1))
  } else {
    plot(mean_volatility, main = "Shock Volatility", ylab = "", lwd = 1, ylim = c(0, max(sub_volatility) + 0.1),
         xaxt = "n")
    axis(1, at = seq(1960, 2010, by = 10), labels = rep("", 6))
  }
} else {
  if(i %% rows == 0) par(mar = c(3, 4.5, 2.5, 1)) else par(mar = c(0, 4.5, 2.5, 1))
  if(i %% rows == 0) {
    plot(mean_volatility, main = "", ylab = y_label, lwd = 1, ylim = c(0, max(sub_volatility) + 0.1))
  } else {
    plot(mean_volatility, main = "", ylab = y_label, lwd = 1, ylim = c(0, max(sub_volatility) + 0.1), xaxt = "n")
    axis(1, at = seq(1950, 2000, by = 10), labels = rep("", 6))
  }
}
grid()
fanplot::fan(data = sub_volatility, data.type = "values", probs = probs,
             start = attributes(y_ts)$tsp[1] + 1, frequency = attributes(y_ts)$tsp[3],
             fan.col = colorRampPalette(c(col, "white")),
             rlab = NULL, ln = NULL)
lines(mean_volatility, lwd = 1, lty = 2)

###########################
### IRFs with aer-bands ###
###########################

irfs_aer <- list("narrative" = list(), "signres" = list())
for(i in 1:6) {
  irfs_aer$narrative[[i]] <- rmatio::read.mat(paste0("./aer_data/irfs_narrative_", i, ".mat"))[[1]]
  irfs_aer$signres[[i]] <- rmatio::read.mat(paste0("./aer_data/irfs_signres_", i, ".mat"))[[1]]
}

irf_plot_with_bands <- function (irf_obj, varnames = NULL, probs = c(0.05, 0.16, 0.84, 0.95),
                                 mar = c(2, 2, 2, 1), mfrow = NULL, color = "tomato",
                                 leg = TRUE, plot_median = FALSE, normalize = NULL,
                                 band_list = NULL,
                                 band_lty = c(2, 3), band_lwd = c(1, 1), band_col = c(1, 1)) {
  if (max(probs) > 0.99 | min(probs) < 0.01) stop("Values of 'probs' need to be between 0.01 and 0.99.")
  probs <- c(probs[1] - 0.01, probs, probs[length(probs)] + 0.01)
  m <- length(irf_obj)
  for(i in 1:length(irf_obj)) if (!is.null(nrow(irf_obj[[i]]))) m <- nrow(irf_obj[[i]])
  if(is.null(varnames)) varnames <- paste0("var", 1:m)
  if(is.null(mfrow)) mfrow <- c(m, m)
  par(mar = mar)
  par(mfrow = mfrow)
  indexmat <- matrix(1:m^2, ncol = m)
  row <- 0
  col <- 0
  count <- 0
  for(fig_index in 1:m^2) {
    if ((fig_index - 1)%%m == 0) {
      row <- row + 1
      col <- 1
    } else {
      col <- col + 1
    }
    if(length(irf_obj) < row) next
    if(is.null(irf_obj[[row]])) next
    count <- count + 1
    sub_irfs <- t(irf_obj[[row]][col, , ])
    if(!is.null(normalize)) sub_irfs <- sub_irfs * normalize[col]
    median_sub_irfs <- ts(apply(sub_irfs, 2, median), start = 0)
    quant <- function(column) quantile(column, probs = probs)
    quantiles_sub_irfs <- apply(sub_irfs, 2, quant)
    if(!is.null(band_list)) {
      ylims <- c(min(c(quantiles_sub_irfs, band_list$narrative[[count]])), # band_list$signres[[count]])),
                 max(c(quantiles_sub_irfs, band_list$narrative[[count]]))) #, band_list$signres[[count]])))
    } else {
      ylims <- c(min(quantiles_sub_irfs), max(quantiles_sub_irfs))
    }
    plot(median_sub_irfs, lwd = 2, lty = 2, col = color,
         ylab = "", xlab = "", main = paste0("Shock ", row," on ", varnames[col]),
         ylim = ylims)
    grid()
    fanplot::fan(data = quantiles_sub_irfs, data.type = "values",
                 probs = probs, start = 0, fan.col = colorRampPalette(c(color, "white")),
                 rlab = NULL, ln = NULL)
    abline(h = 0, lwd = 2, lty = 2)
    if(plot_median) lines(median_sub_irfs, lwd = 2, lty = 1)
    post_mass <- (max(probs[-length(probs)]) - min(probs[-1])) * 100
    if(col == 1 & row == 1 & leg == TRUE) legend("topleft", c(paste0(post_mass, "% of post. prob. mass")), lwd = 0, bty = "n", col = "tomato")
    if(!is.null(band_list)) {
      lines(ts(band_list$narrative[[count]][,1]), lty = band_lty[1], lwd = band_lwd[1], col = band_col[1])
      lines(ts(band_list$narrative[[count]][,3]), lty = band_lty[1], lwd = band_lwd[1], col = band_col[1])
      #lines(ts(band_list$signres[[count]][,1]), lty = band_lty[2], lwd = band_lwd[2], col = band_col[2])
      #lines(ts(band_list$signres[[count]][,3]), lty = band_lty[2], lwd = band_lwd[2], col = band_col[2])
    }
  }
  par(mfrow = c(1, 1))
}

irf_plot_with_bands(irf_obj_mon,
                    varnames = c("GDP", "Deflator", "CPI", "TotRes", "BorRes", "Fed"),
                    leg = FALSE,
                    mar = c(2,2,2,1),
                    mfrow = c(2,3),
                    normalize = c(rep(100, 6), 1),
                    band_list = irfs_aer,
                    band_lty = c(1, 1), band_lwd = c(2, 2), band_col = c("steelblue", "darkgrey"))
#dev.off()

################################################
### Marginal posteriors of on impact effects ###
################################################

N_post <- nrow(output_uhlig$chains)
post_binv <- output_uhlig$chains[(N_post - (N - burn) * output_uhlig$args$n + 1):N_post, grep("b", colnames(output_uhlig$chains))]
post_b <- post_binv
for(i in 1:nrow(post_b)) {
  post_b[i,] <- c(solve(matrix(post_binv[i,], ncol = 6)))
}

#pdf("./empirical_application_Uhlig/on_impact_all.pdf", width = 11, height = 9)
imat <- t(matrix(1:(6*6), ncol = 6))
par(mfrow = c(6, 6))
par(mar = c(2, 4, 3.5, 1))
row <- 1
col <- 0
for(i in 1:(6*6)) {
  col <- col + 1
  if(col > 6) {
    col <- 1
    row <- row + 1
  }
  hist(post_b[,imat[i]] * 100, probability = TRUE,
       main = paste0("Shock ", row, " on ", c("GDP", "Deflator", "CPI", "TotRes", "BorRes", "Fed")[col]),
       ylab = ifelse(col == 1, "Density", ""),
       xlab = "")
}
par(mfrow = c(6, 6))
#dev.off()






