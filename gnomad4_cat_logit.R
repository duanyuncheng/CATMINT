# load library
```{r}
library(tidyverse)
#library(brglm2)
#library(VGAM)
library(nnet)
```

```{r}
#library(devtools)
devtools::install_local("C:/Users/yuncheng/Downloads/ACAT-master.zip")
library(ACAT)
```

# pseudocounts with nnet (high/medium/low confidence)
## excluding zero expected counts
```{r}
# filter on the expected counts
# ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_summary file from C:/Users/yuncheng/Documents/R/Joint/summary_of_gnomad4_exome_5_classes_LRT_result.Rmd
ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset <- ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset_poisson <- ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_summary %>% filter(., sum_pr_non_functional > 0 & sum_pr_OP > 0 & sum_pr_no_probably_damaging_missense > 0 & sum_pr_missense_probably_damaging > 0 & sum_pr_lof > 0)

# add pseudo counts to observation
ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_non_functional <- ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_non_functional + ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_non_functional_to_all

ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_OP <- ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_OP + ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_OP_to_all

ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_no_probably_damaging_missense <- ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_no_probably_damaging_missense + ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_no_probably_damaging_to_all

ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_missense_probably_damaging <- ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_missense_probably_damaging + ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_missense_probably_damaging_to_all

ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_lof <- ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_lof + ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_lof_to_all

ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_all <- ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_all + 1

#observation ratio recalculate
ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_non_functional_to_all <-  ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_non_functional/ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_all

ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_OP_to_all <- ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_OP/ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_all

ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_no_probably_damaging_to_all <- ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_no_probably_damaging_missense/ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_all

ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_missense_probably_damaging_to_all <-
ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_missense_probably_damaging/ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_all

ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_lof_to_all <-
ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_lof/ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$num_all

```

```{r}
ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_alpha_COC <- log(ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_OP_to_all/ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_non_functional_to_all)

ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_alpha_LDM <- log(ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_no_probably_damaging_to_all/ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_non_functional_to_all)

ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_alpha_PDM <- log(ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_missense_probably_damaging_to_all/ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_non_functional_to_all)

ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_alpha_LOF <- log(ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_lof_to_all/ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_non_functional_to_all)



ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_alpha_COC <- log(ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_OP_to_all/ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_non_functional_to_all)

ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_alpha_LDM <- log(ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_no_probably_damaging_to_all/ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_non_functional_to_all)

ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_alpha_PDM <- log(ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_missense_probably_damaging_to_all/ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_non_functional_to_all)

ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_alpha_LOF <- log(ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_lof_to_all/ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_non_functional_to_all)
```

```{r}
ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$logit_delta_COC <- ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_alpha_COC - ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_alpha_COC

ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$logit_delta_LDM <- ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_alpha_LDM - ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_alpha_LDM

ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$logit_delta_PDM <- ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_alpha_PDM - ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_alpha_PDM

ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$logit_delta_LOF <- ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$obs_alpha_LOF - ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset$sim_alpha_LOF
```

## run baseline category logit with nnet
### file prepare
```{r}
gnomad4_exome_5_classes_file_for_logit_fullset <- ensembl_sim_obs_for_multinonmial_gnomad4_exome_5_classes_logit_fullset %>% select(., gene_name, num_non_functional, num_OP, num_no_probably_damaging_missense, num_missense_probably_damaging, num_lof, sim_non_functional_to_all, sim_OP_to_all, sim_no_probably_damaging_to_all, sim_missense_probably_damaging_to_all, sim_lof_to_all, obs_non_functional_to_all, obs_OP_to_all, obs_no_probably_damaging_to_all, obs_missense_probably_damaging_to_all, obs_lof_to_all, sim_alpha_COC, sim_alpha_LDM, sim_alpha_PDM, sim_alpha_LOF, label)
colnames(gnomad4_exome_5_classes_file_for_logit_fullset) <- c("gene_name","X1","X2","X3","X4","X5", "pi_01","pi_02","pi_03","pi_04","pi_05","pi_11","pi_12","pi_13","pi_14","pi_15","alpha_2","alpha_3","alpha_4","alpha_5", "label")
head(gnomad4_exome_5_classes_file_for_logit_fullset)
```

### run logit and get estimates and z-scores for four categories
```{r}
gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset <- matrix(NA, nrow(gnomad4_exome_5_classes_file_for_logit_fullset), 10)
colnames(gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset) <- c("gene_name", "delta_COC","delta_BDM","delta_PDM","delta_LOF","zscore_COC", "zscore_BDM", "zscore_PDM", "zscore_LOF", "label")
for (i in 1:nrow(gnomad4_exome_5_classes_file_for_logit_fullset))
{
  gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset[i,1] <- gnomad4_exome_5_classes_file_for_logit_fullset$gene_name[i]
  gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset[i,10] <- gnomad4_exome_5_classes_file_for_logit_fullset$label[i]
  data_mat <- matrix(unlist(gnomad4_exome_5_classes_file_for_logit_fullset[i,c("X2","X3","X4","X5","X1")]), ncol = 5)
  colnames(data_mat) <- c("cat1", "cat2", "cat3", "cat4", "baseline")
  alpha <- unlist(gnomad4_exome_5_classes_file_for_logit_fullset[i,c("alpha_2","alpha_3","alpha_4","alpha_5")])
  offset_mat <- matrix(alpha, nrow = length(alpha))
  
  # Convert the matrix to a multinom dataframe
  data_multinom <- data.frame(
  category = factor(colnames(data_mat)),
  count = data_mat[1,])
  data_multinom$category <- relevel(data_multinom$category, ref = "baseline")
  
  # Fit the model with nnet
  fit_multinom <- multinom(category ~ 1, data = data_multinom, weights = data_multinom$count)
  deltas <- coefficients(fit_multinom) - offset_mat
  gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset[i,2] <- deltas[1]
  gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset[i,3] <- deltas[2]
  gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset[i,4] <- deltas[3]
  gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset[i,5] <- deltas[4]
  
  #  Compute Zscores
  zscores <- deltas/sqrt(diag(vcov(fit_multinom)))
  gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset[i,6] <- zscores[1]
  gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset[i,7] <- zscores[2]
  gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset[i,8] <- zscores[3]
  gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset[i,9] <- zscores[4]
}
```

```{r}
df_gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset <- data.frame(gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset)

df_gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset$delta_COC <- as.numeric(df_gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset$delta_COC)

df_gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset$delta_BDM <- as.numeric(df_gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset$delta_BDM)

df_gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset$delta_PDM <- as.numeric(df_gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset$delta_PDM)

df_gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset$delta_LOF <- as.numeric(df_gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset$delta_LOF)

df_gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset$zscore_COC <- as.numeric(df_gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset$zscore_COC)

df_gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset$zscore_BDM <- as.numeric(df_gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset$zscore_BDM)

df_gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset$zscore_PDM <- as.numeric(df_gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset$zscore_PDM)

df_gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset$zscore_LOF <- as.numeric(df_gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset$zscore_LOF)
```

### run logit and sum of zscore_PDM and zscore_LOF
```{r}
gnomad4_exome_5_classes_file_for_nf_enrichment <- matrix(NA, nrow(gnomad4_exome_5_classes_file_for_logit_fullset), 7)
colnames(gnomad4_exome_5_classes_file_for_nf_enrichment) <- c("gene_name", "zscore_NF", "pvalue_NF", "zscore_PDM", "zscore_LOF", "Var_Y", "label")
for (i in 1:nrow(gnomad4_exome_5_classes_file_for_logit_fullset))
{
  gnomad4_exome_5_classes_file_for_nf_enrichment[i,1] <- gnomad4_exome_5_classes_file_for_logit_fullset$gene_name[i]
  gnomad4_exome_5_classes_file_for_nf_enrichment[i,7] <- gnomad4_exome_5_classes_file_for_logit_fullset$label[i]
  data_mat <- matrix(unlist(gnomad4_exome_5_classes_file_for_logit_fullset[i,c("X2","X3","X4","X5","X1")]), ncol = 5)
  colnames(data_mat) <- c("cat1", "cat2", "cat3", "cat4", "baseline")
  alpha <- unlist(gnomad4_exome_5_classes_file_for_logit_fullset[i,c("alpha_2","alpha_3","alpha_4","alpha_5")])
  offset_mat <- matrix(alpha, nrow = length(alpha))
  
  # Convert the matrix to a multinom dataframe
  data_multinom <- data.frame(
  category = factor(colnames(data_mat)),
  count = data_mat[1,])
  data_multinom$category <- relevel(data_multinom$category, ref = "baseline")
  
  # Fit the model with nnet
  fit_multinom <- multinom(category ~ 1, data = data_multinom, weights = data_multinom$count)
  deltas <- coefficients(fit_multinom) - offset_mat
  zscores <- deltas/sqrt(diag(vcov(fit_multinom)))
  zscore3 <- zscores[3]
  zscore4 <- zscores[4]
  
  Y <- zscore3 + zscore4
  Var_Y <- 2 + 2*vcov(fit_multinom)[3,4]/sqrt(vcov(fit_multinom)[3,3]*vcov(fit_multinom)[4,4])

  
  gnomad4_exome_5_classes_file_for_nf_enrichment[i,2] <- Y/sqrt(Var_Y)
  gnomad4_exome_5_classes_file_for_nf_enrichment[i,3] <- pnorm(Y/sqrt(Var_Y))
  gnomad4_exome_5_classes_file_for_nf_enrichment[i,4] <- zscore3
  gnomad4_exome_5_classes_file_for_nf_enrichment[i,5] <- zscore4
  gnomad4_exome_5_classes_file_for_nf_enrichment[i,6] <- Var_Y
}
```

```{r}
df_gnomad4_exome_5_classes_file_for_nf_enrichment <- data.frame(gnomad4_exome_5_classes_file_for_nf_enrichment)
df_gnomad4_exome_5_classes_file_for_nf_enrichment$zscore_NF <- as.numeric(df_gnomad4_exome_5_classes_file_for_nf_enrichment$zscore_NF)
df_gnomad4_exome_5_classes_file_for_nf_enrichment$pvalue_NF <- as.numeric(df_gnomad4_exome_5_classes_file_for_nf_enrichment$pvalue_NF)
df_gnomad4_exome_5_classes_file_for_nf_enrichment$zscore_PDM <- as.numeric(df_gnomad4_exome_5_classes_file_for_nf_enrichment$zscore_PDM)
df_gnomad4_exome_5_classes_file_for_nf_enrichment$zscore_LOF <- as.numeric(df_gnomad4_exome_5_classes_file_for_nf_enrichment$zscore_LOF)
df_gnomad4_exome_5_classes_file_for_nf_enrichment$Var_Y <- as.numeric(df_gnomad4_exome_5_classes_file_for_nf_enrichment$Var_Y)
```

### Cauchy Test for Multiple P-values
```{r}
df_gnomad4_exome_5_classes_file_pvalue_combination <- df_gnomad4_exome_5_classes_file_for_four_logit_zscore_fullset

df_gnomad4_exome_5_classes_file_pvalue_combination$pvalue_PDM <- 1-pnorm(df_gnomad4_exome_5_classes_file_pvalue_combination$zscore_PDM)

df_gnomad4_exome_5_classes_file_pvalue_combination$pvalue_LOF <- 1 - pnorm(df_gnomad4_exome_5_classes_file_pvalue_combination$zscore_LOF)

df_gnomad4_exome_5_classes_file_pvalue_combination$pvalue_PDMLOF <- 1 - df_gnomad4_exome_5_classes_file_for_nf_enrichment$pvalue_NF
```

```{r}
input_matrix_for_ACAT <- as.matrix(df_gnomad4_exome_5_classes_file_pvalue_combination[,c("pvalue_PDM", "pvalue_LOF", "pvalue_PDMLOF")])
input_matrix_for_ACAT <- t(input_matrix_for_ACAT)
```

```{r}
df_gnomad4_exome_5_classes_file_pvalue_combination$pvalue_combined <- ACAT(input_matrix_for_ACAT)
```





