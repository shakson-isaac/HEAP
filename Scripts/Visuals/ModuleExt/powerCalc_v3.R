library(pwr)
library(ggplot2)

# # Assuming 10 covariates (e.g., age, sex, PCs), 1 SNP as predictors
# R2 <- 0.01
# ans1 <- pwr.f2.test(u = 11,         # 10 covariates + 1 SNP
#             v = NULL,       # Unknown, to be computed
#             f2 = R2/(1-R2),    # Calculated effect size
#             sig.level = 5e-8, # Genome-wide significance
#             power = 0.8)  
# plot(ans1)


#Power for GWAS/EWAS type studies using R2 of the specific exposure
library(viridis)
plot_power <- function(R2_end = 0.1, N_samp = 10000, N_cov = 20, 
                       sig = 5e-8, plot_title = "Power Analysis",
                       R2_start = 0.01,
                       R2_iter = 0.01){
  # Define parameters
  effect_sizes <- seq(R2_start, R2_end, by = R2_iter) # Varying R2
  sample_sizes <- seq(100, N_samp, by = 50)   # Sample sizes
  covariates <- N_cov                         # Number of predictors (e.g., Covariates + Main Predictor)
  
  # Store parameters
  power_results <- expand.grid(effect_size = effect_sizes, sample_size = sample_sizes)
  power_results$power <- NA
  
  # Calculate power for each combo
  for (i in 1:nrow(power_results)) {
    R2 <- power_results$effect_size[i]
    n <- power_results$sample_size[i]
    v <- n - covariates - 1
    
    # Calculate power
    power_results$power[i] <- pwr.f2.test(u = covariates, v = v, f2 = R2/(1-R2), sig.level = sig)$power
  }

  # Plot power vs. sample size (for diff. effect size)
  gg1 <- ggplot(power_results, aes(x = sample_size, y = power, color = as.factor(effect_size))) +
    geom_line(linewidth = 1) +
    labs(title = paste0(plot_title),
         x = "Sample Size",
         y = "Power",
         color = bquote("Effect Size"~R^2)) +
    theme_minimal() +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red", linewidth = 1) + # Reference line for 80% power
  theme_minimal(base_size = 14) + # Use a minimal theme with a larger base text size
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"), # Centered bold title
      #axis.title = element_text(face = "bold"),   # Bold axis titles
      legend.position = "bottom",                 # Move legend to bottom
      #legend.title = element_text(face = "bold"), # Bold legend title
      #legend.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5), # Styled legend box
      legend.key.size = unit(1, "cm")             # Adjust legend key size
    ) +
    # Add text annotation for significance level and number of covariates
    annotate("text", x = max(sample_sizes) * 0.6, y = 0.1, 
             label = paste("p-value threshold:", sig, "\nNumber of covariates:", N_cov), 
             hjust = 0, vjust = 0, size = 3, color = "black", 
             fontface = "italic")
  
  # Plot power vs. sample size for different effect sizes
  # gg1 <- ggplot(power_results, aes(x = sample_size, y = power, color = effect_size)) +
  #   geom_line(linewidth = 1.2) +  # Use linewidth instead of size
  #   #geom_point(size = 2) +        # Point size adjustment
  #   scale_color_viridis_c() +     # Use viridis color scale
  #   labs(title = paste0(plot_title, ": pval thresh (", sig, ") #covar (", N_cov, ")"),
  #        x = "Sample Size",
  #        y = "Power",
  #        color = bquote("Effect Size"~R^2)) +
  #   geom_hline(yintercept = 0.8, linetype = "dashed", color = "red", linewidth = 1) + # Reference line for 80% power
  #   theme_minimal(base_size = 14) + # Use a minimal theme with a larger base text size
  #   theme(
  #     plot.title = element_text(hjust = 0.5, face = "bold"), # Centered bold title
  #     axis.title = element_text(face = "bold"),   # Bold axis titles
  #     legend.position = "bottom",                 # Move legend to bottom
  #     legend.title = element_text(face = "bold"), # Bold legend title
  #     legend.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5), # Styled legend box
  #     legend.key.size = unit(1, "cm")             # Adjust legend key size
  #   )
  return(gg1)
  
}
plot_power(plot_title = "Power Analysis GWAS variants")

#Under bonferroni correction:
#0.05/(3000*250) ~ 7e-8
plot_power(N_cov = 30, sig = round(7e-8, digits = 10), plot_title = "Power Analysis Exposomic Architecture \n of the Proteome")




plot_power(N_cov = 10, sig = round(7e-8, digits = 10), plot_title = "Power Analysis Exposomic Architecture \n of the Proteome")


plot_power(N_cov = 10, 
           N_samp = 20000,
           sig = round(7e-8, digits = 10), 
           plot_title = "Power Analysis Exposomic Architecture \n of the Proteome")




plot_power(R2_end = 0.1,
           N_cov = 30, 
           N_samp = 50000,
           sig = round(7e-8, digits = 10), 
           plot_title = "Power Analysis Exposomic Architecture \n of the Proteome",
           R2_start = 1e-3,
           R2_iter = 0.005)


plot_powerv2 <- function(effect_sizes, N_samp = 10000, N_cov = 20, 
                       sig = 5e-8, plot_title = "Power Analysis"){
  # Define parameters:
  #User defined: effect_sizes - Varying R2
  sample_sizes <- seq(100, N_samp, by = 50)   # Sample sizes
  covariates <- N_cov                         # Number of predictors (e.g., Covariates + Main Predictor)
  
  # Store parameters
  power_results <- expand.grid(effect_size = effect_sizes, sample_size = sample_sizes)
  power_results$power <- NA
  
  # Calculate power for each combo
  for (i in 1:nrow(power_results)) {
    R2 <- power_results$effect_size[i]
    n <- power_results$sample_size[i]
    v <- n - covariates - 1
    
    # Calculate power
    power_results$power[i] <- pwr.f2.test(u = covariates, v = v, f2 = R2/(1-R2), sig.level = sig)$power
  }
  
  # Plot power vs. sample size (for diff. effect size)
  gg1 <- ggplot(power_results, aes(x = sample_size, y = power, color = as.factor(effect_size))) +
    geom_line(linewidth = 1) +
    labs(title = paste0(plot_title),
         x = "Sample Size",
         y = "Power",
         color = bquote("Effect Size"~R^2)) +
    theme_minimal() +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red", linewidth = 1) + # Reference line for 80% power
    theme_minimal(base_size = 14) + # Use a minimal theme with a larger base text size
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"), # Centered bold title
      #axis.title = element_text(face = "bold"),   # Bold axis titles
      legend.position = "bottom",                 # Move legend to bottom
      #legend.title = element_text(face = "bold"), # Bold legend title
      #legend.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5), # Styled legend box
      legend.key.size = unit(1, "cm")             # Adjust legend key size
    ) +
    # Add text annotation for significance level and number of covariates
    annotate("text", x = max(sample_sizes) * 0.6, y = 0.2, 
             label = paste("p-value threshold:", sig, "\nNumber of covariates:", N_cov), 
             hjust = 0, vjust = 0, size = 3, color = "black", 
             fontface = "italic")
  
  # Plot power vs. sample size for different effect sizes
  # gg1 <- ggplot(power_results, aes(x = sample_size, y = power, color = effect_size)) +
  #   geom_line(linewidth = 1.2) +  # Use linewidth instead of size
  #   #geom_point(size = 2) +        # Point size adjustment
  #   scale_color_viridis_c() +     # Use viridis color scale
  #   labs(title = paste0(plot_title, ": pval thresh (", sig, ") #covar (", N_cov, ")"),
  #        x = "Sample Size",
  #        y = "Power",
  #        color = bquote("Effect Size"~R^2)) +
  #   geom_hline(yintercept = 0.8, linetype = "dashed", color = "red", linewidth = 1) + # Reference line for 80% power
  #   theme_minimal(base_size = 14) + # Use a minimal theme with a larger base text size
  #   theme(
  #     plot.title = element_text(hjust = 0.5, face = "bold"), # Centered bold title
  #     axis.title = element_text(face = "bold"),   # Bold axis titles
  #     legend.position = "bottom",                 # Move legend to bottom
  #     legend.title = element_text(face = "bold"), # Bold legend title
  #     legend.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5), # Styled legend box
  #     legend.key.size = unit(1, "cm")             # Adjust legend key size
  #   )
  return(gg1)
  
}
plot_powerv2(effect_sizes = c(1e-3, 2.5e-3, 5e-3,
                              7.5e-3, 1e-2, 2.5e-2, 
                              5e-2, 7.5e-2, 1e-1),
           N_samp = 50000,
           N_cov = 30,
           sig = round(7e-8, digits = 10), 
           plot_title = "Power Analysis of HEAP Associations")

#'*SAVE THE PLOT LATERRRRR*
#'*SAVE THE PLOT!!!*










#'*GxE signal not well powered - not sure how tiny is ok*
plot_power(R2_end = 0.01, 
           N_samp = 50000,
           N_cov = 30, 
           sig = round(0.05/(3000*250), digits = 10), 
           plot_title = "Power Analysis GxE Architecture \n of the Proteome",
           R2_start = 1e-3,
           R2_iter = 1e-3)

#Power for Polygenic GxE signals: just base it off R2. 



#Power in Longitudinal Framework:
# Load the longpower package
#library(longpower)
# Load necessary libraries
library(lme4)  # For linear mixed-effects models
library(simr)  # For power analysis through simulation

# Define parameters
n_subjects <- 30         # Number of subjects
n_timepoints <- 2        # Number of time points (repeated measures)
effect_size <- 0.5       # Cohen's d for the fixed effect
sigma_within <- 1        # Within-subject standard deviation
sigma_between <- 2       # Between-subject standard deviation
alpha <- 0.05            # Significance level

# Simulate data
set.seed(123)
data <- expand.grid(
  subject = 1:n_subjects,
  time = 1:n_timepoints
)
data$response <- with(data, rnorm(n = n_subjects * n_timepoints,
                                  mean = effect_size * time,
                                  sd = sigma_within))

# Fit a linear mixed-effects model
model <- lmer(response ~ time + (1 | subject), data = data)

# Power analysis using simulation
power_result <- powerSim(model, nsim = 10)  # Run 1000 simulations

# Display the results
summary(power_result)

# Extend the model to allow for power analysis
model_power <- extend(model, along = "subject", n = seq(10, 50, by = 10))  # Vary number of subjects

# Perform power analysis using powerCurve
power_curve <- powerCurve(model_power, nsim = 100)
plot(power_curve)

#Good review of lmer and power simulations here:
#https://humburg.github.io/Power-Analysis/simr_power_analysis.html
#use LRT to get pvals.




