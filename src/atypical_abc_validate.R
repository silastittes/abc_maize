library(abc)
library(tidyverse)
library(ggridges)
library(patchwork)
library(cowplot)
library(magrittr)
theme_set(theme_cowplot(font_size = 20))
par(ask=F)

n_ind <- 26
fold_bins <- (1:(n_ind/2))
to_keep <- c(fold_bins, fold_bins + n_ind*1, fold_bins + n_ind*2)
#var0 <- function(x) var(x) > 0


full_df <- 
  vroom::vroom(
    "NAM_DFE_ABC/pop_sfs.txt", delim = " ",
    col_names = FALSE
  ) %>%
  dplyr::select(-ncol(.)) %>% 
  dplyr::select(all_of(to_keep)) %>% 
  mutate(sm  = rowSums(.)) %>%
  rowwise() %>%
  mutate_all(.funs = ~ .x / sm) %>%
  dplyr::select(-sm) %>% 
  ungroup()


atypical_df <- 
  vroom::vroom(
    "NAM_DFE_ABC/atypical_pop_sfs.txt", delim = " ",
    col_names = FALSE
  ) %>%
  dplyr::select(-ncol(.)) %>% 
  dplyr::select(all_of(to_keep)) %>% 
  mutate(sm  = rowSums(.)) %>%
  rowwise() %>%
  mutate_all(.funs = ~ .x / sm) %>%
  dplyr::select(-sm) %>% 
  ungroup()

pcx <- prcomp(full_df, scale. = F, center = F)
plot(pcx)
#plot(pcx$x[,1], pcx$x[,2], pch = ".")
varz <- cumsum((pcx$sdev ^ 2) / sum(pcx$sdev ^ 2)) <= 0.99
mx_col <- max(c(2,length(varz[varz])))
pc_df <- data.frame(pcx$x) %>% 
  dplyr::select(1:mx_col)


params <- 
  c("Na", "N0", "Nb", "B_t", "mu_sv",
    "sfs1_shape", "sfs1_mean", "sfs2_shape", "sfs2_mean")

param_df_full <- read_delim(
  "NAM_DFE_ABC/pop.txt",
  col_names = params,
  delim = " "
) %>% 
  dplyr::select(-c(Na, B_t))


nrow(param_df_full)
q <- 0.3
param_df_full %>% 
  filter(mu_sv < quantile(mu_sv, q), sfs2_mean < quantile(sfs2_mean, q)) %>% 
  nrow()


atypical_param_df <- 
  read_delim(
    "NAM_DFE_ABC/atypical_pop.txt",
    col_names = c("ngenes", "Na", "N0", "Nb", "B_t", "mu_sv",
                  "sfs1_shape", "sfs1_mean", "sfs2_shape", "sfs2_mean"),
    delim = " ", comment = "#"
  ) %>% 
  dplyr::select(-c(Na, B_t))


itts <- nrow(atypical_df)
p <- progress_estimated(itts)



set.seed(123)
val_df <- 
  1:itts %>% 
  map_df(function(target_idx){
    print(p$tick())
    
    # atypical_param_df %>% 
    #   mutate(idx = 1:n()) %>% 
    # filter(sfs2_mean == 0 | sfs1_mean == 0)
    
    #raw
    target_stats <- atypical_df[target_idx, ]
    sumstat_df <- full_df
    
    target_df <- atypical_param_df[target_idx, ]
    param_df <-  param_df_full
    
    
    res <- abc(target=target_stats,
               param=param_df,
               sumstat=sumstat_df,
               tol=0.002,
               transf=c("none"),
               method = "neuralnet", 
               sizenet = 2
    )
    
    
    posts_df <- data.frame(res$adj.values)
    
    names(posts_df) %>% 
      map_df(function(c_name){
        prop_gt <- mean(posts_df[[c_name]] >= target_df[[c_name]])
        creds <- quantile(posts_df[[c_name]], c(0,0.025,0.975,1))
        w.in_cred <- as.numeric(creds[2] < target_df[[c_name]] & creds[3] > target_df[[c_name]])
        w.in_prior <- as.numeric(creds[1] < target_df[[c_name]] & creds[4] > target_df[[c_name]])
        diff_sc <- mean( (posts_df[[c_name]] - target_df[[c_name]]) / (target_df[[c_name]] + 1) )
        var_sc <- sd(posts_df[[c_name]] / sd(param_df[[c_name]]))
        mean_diff <- mean(sqrt((posts_df[[c_name]] -  target_df[[c_name]])^2) /( target_df[[c_name]] + 1))
        mean_post <- mean(posts_df[[c_name]])
        ngenes <- target_df$ngenes
        tibble(mean_diff, parameter = c_name, w.in_cred, w.in_prior, prop_gt, diff_sc, var_sc, mean_post)
      }) %>% 
      bind_cols(. , target_df)
  })


make_plots <- function(var, val, logp = TRUE){
  
  var <- enquo(var)
  
  if(logp){
    dff <- 
      val_df %>% 
      filter(sfs2_mean != 0, sfs1_mean !=0) %>% 
      mutate(ngenes = log10(ngenes))  
    median_ngenes <- log10(300)
    low <-  log10(88)
    high <- log10(720)
    xtitle <- "log10(number of genes)"
  } else {
    dff <- val_df
    median_ngenes <- 300
    low <-  88
    high <- 720
    xtitle <- "number of genes"
  }
  dff %>% 
    filter(parameter == val) %>%
    ggplot(aes(ngenes, !!var)) +
    geom_jitter(width = 0, height = 0.0) +
    geom_smooth(se = F, colour = "black", alpha = 0.7) +
    geom_vline(xintercept = median_ngenes, lty = 1) +
    geom_vline(xintercept = low, lty = 2) +
    geom_vline(xintercept = high, lty = 2) +
    xlab(xtitle)
}


val_df %<>% mutate(is_zero = factor(if_else(sfs1_mean == 0 | sfs2_mean == 0, true = "zero", false = "not_zero")))

param <- "sfs1_mean"
make_plots(w.in_cred, param) +
  ggtitle(param) +
  make_plots(mean_diff, param) +
  make_plots(log(var_sc), param) +
  ggsave(filename = str_glue("NAM_DFE_ABC/figures/atypical_{param}.png"))

param <- "sfs2_mean"
make_plots(w.in_cred, param) +
  ggtitle(param)+
  make_plots(mean_diff, param) +
  make_plots(log(var_sc), param) +
  ggsave(filename = str_glue("NAM_DFE_ABC/figures/atypical_{param}.png"))



val_df %>% 
  filter(parameter %in% c("sfs1_mean", "sfs2_mean"),
         mu_sv > -INF) %>%
  mutate(var = mean_post) %>%
  ggplot(aes(parameter, -var)) +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.05) +
  facet_grid(~is_zero)

val_df %>% 
  filter(parameter %in% c("sfs1_mean", "sfs2_mean"), 
         sfs2_mean == 0, sfs1_mean == 0,
         mu_sv > 1e-8) %>%
  mutate(var = mean_post) %>% 
  group_by(parameter) %>% 
  summarise(mean(mean_post))


val_df %>% 
  filter(parameter %in% c("sfs2_mean")) %>%
  mutate(var = mean_post) %>%
  ggplot(aes(mu_sv, -var, colour = is_zero)) +
  geom_point() +
  geom_smooth()


val_df %>% 
  filter(parameter %in% c("sfs2_mean", "mu_sv")) %>%
  group_by(parameter, is_zero) %>% 
  mutate(idx = 1:n()) %>% 
  select(idx, parameter, mean_post) %>% 
  pivot_wider(id_cols = c(idx, is_zero), names_from = "parameter", values_from = mean_post) %>% 
  ggplot(aes(mu_sv, -sfs2_mean, colour = is_zero)) +
  geom_point() +
  geom_smooth()



val_df %>% 
  filter(parameter %in% c("sfs2_mean", "mu_sv")) %>%
  group_by(parameter, is_zero) %>% 
  mutate(idx = 1:n()) %>% 
  select(idx, parameter, mean_post) %>% 
  pivot_wider(id_cols = c(idx, is_zero), names_from = "parameter", values_from = mean_post) %>% 
  ggplot(aes(sfs2_mean, fill = is_zero)) +
  geom_density()


val_df %>% 
  filter(parameter %in% c("mu_sv")) %>%
  mutate(var = mean_post) %>%
  ggplot(aes(mu_sv, var, colour = is_zero)) +
  geom_point() +
  geom_smooth()


val_df %>% 
  filter(parameter %in% c("sfs2_mean")) %>%
  ggplot(aes(log10(ngenes), log(var_sc), colour = is_zero)) +
  geom_point() +
  geom_smooth()


val_df %>% 
  filter(parameter %in% c("sfs2_mean")) %>%
  group_by(is_zero) %>% 
  summarise(mean(w.in_cred))
  
val_df %>% 
  filter(parameter %in% c("sfs2_mean")) %>%
  ggplot(aes(mu_sv, w.in_prior, colour = is_zero)) +
  geom_jitter(width = 0, height = 0.05) +
  geom_smooth()

val_df %>% filter(parameter %in% c("sfs2_mean"), sfs2_mean == 0) %>% pull(ngenes)

val_df %>% 
  filter(parameter %in% c("sfs2_mean"), sfs2_mean == 0) %>% 
  lm(-mean_post ~ N0 + Nb + mu_sv + sfs2_shape + sfs1_shape, data = .) %>% 
  summary()


val_df %>% 
  filter(parameter %in% c("sfs2_mean"), sfs2_mean != 0) %>% 
  lm(-mean_post ~ N0 + Nb + mu_sv + sfs2_shape + sfs1_shape, data = .) %>% 
  summary()


