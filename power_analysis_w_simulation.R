# load library
```{r}
library(tidyverse)
```

# function likelihood ratio test
```{r}
computeLRT <- function(observed_counts, expected_proportions_null) {
  # Calculate total number of observations
  n <- sum(observed_counts)
  
  # Calculate the log-likelihood for the null model 
  logL_null <- sum(dmultinom(observed_counts, size = n, prob = expected_proportions_null, log = TRUE))
  
  # Calculate the log-likelihood for the alternative model (MLE)
  observed_proportions_alternative <- observed_counts / n
  logL_alternative <- sum(dmultinom(observed_counts, size = n, prob = observed_proportions_alternative, log = TRUE))
  
  # Calculate the test statistic
  test_statistic <- 2 * (logL_alternative - logL_null)
  
  # Degrees of freedom and P-value
  df <- length(observed_counts) - 1 
  p_value <- pchisq(test_statistic, df, lower.tail = FALSE)

  return(list(test_statistic = test_statistic, p_value = p_value))
}

# Test
observed_counts <- c(50, 25, 15, 5, 5) 
expected_proportions_null <- c(0.2, 0.2, 0.2, 0.2, 0.2) 

lrt_result <- computeLRT(observed_counts, expected_proportions_null)
print(lrt_result)

```

# function computing KL divergence
```{r}
calculate_KL_divergence <- function(P, Q) {

  # KL divergence
  KL_divergence <- sum(P * log(P / Q))
  
  return(KL_divergence)
}

# Test
P <- c(0.2, 0.4, 0.1, 0.2, 0.1) 
Q <- c(0.1, 0.3, 0.3, 0.1, 0.2)
KL_result <- calculate_KL_divergence(P, Q)

print(KL_result)
```

# simulation examples
## example 1 VAMP3 n_variant = 98                                                                                                  
```{r}
set.seed(1234)
i <- 4
n_sample <- sum(file_for_logit_fullset2_poisson[i,2:6]) 
p_null <- unlist(file_for_logit_fullset2_poisson[i,7:11])
alpha <- 0.05
n_simulation <- 1000

power_result <- c()
KL_divergence_result <- c()
index <- 1
for (random_index in sample(1:nrow(file_for_logit_fullset2_poisson), 100, replace = FALSE))
{
  p_alt <- unlist(file_for_logit_fullset2_poisson[random_index,12:16])
  KL_result <- calculate_KL_divergence(p_alt, p_null)
  
  simulated_data_matrix <- rmultinom(n_simulation, n_sample, p_alt)
  successes <- 0
  for (j in 1:n_simulation)
  {
    simulated_data <- unlist(simulated_data_matrix[,j])
    lrt_result <- computeLRT(simulated_data, p_null)
    lrt_pvalue <- lrt_result$p_value
    if (lrt_pvalue < alpha)
    {
      successes <- successes + 1
    }
  }
  power <- successes/n_simulation
  power_result <- c(power_result, power)
  KL_divergence_result <- c(KL_divergence_result, KL_result)
  index <- index + 1
}

```




