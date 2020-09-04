library(abc)
library(tidyverse)
library(ggridges)
library(patchwork)
library(cowplot)
theme_set(theme_cowplot(font_size = 25))
par(ask=F)

n_ind <- 26
fold_bins <- (1:(n_ind/2))
to_keep <- c(fold_bins, fold_bins + n_ind*1, fold_bins + n_ind*2)
#var0 <- function(x) var(x) > 0

abc_val <- function(sfs_file, params_file){
  full_df <- vroom::vroom(
    #"NAM_DFE_ABC/pop_sfs_fixedMu.txt", delim = " ",
    sfs_file, delim = " ",
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
    #"NAM_DFE_ABC/pop_fixedMu.txt",
    params_file,
    col_names = params,
    delim = " "
  ) %>% 
    dplyr::select(-c(Na, B_t))
  
  
  idx_0 <- which(param_df_full$sfs1_mean != 0)
  
  itts <- 100
  p <- progress_estimated(itts)
  
  set.seed(123)
  
  val_df <- 
    sample(seq_len(nrow(full_df)), itts) %>% 
    map_df(function(target_idx){
      print(p$tick())
      
      #target_idx <- 73
      #raw
      target_stats <- full_df[target_idx, ]
      sumstat_df <- full_df[-target_idx, ]
      
      target_df <- param_df_full[target_idx, ]
      param_df <-  param_df_full[-target_idx, ]
      
      
      res <- abc(target=target_stats,
                 #param=select(param_df[idx_0, ], -mu_sv),
                 #sumstat=sumstat_df[idx_0, ],
                 param=select(param_df, -mu_sv),
                 sumstat=sumstat_df,
                 #param= param_df,
                 #sumstat=sumstat_df,
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
          #ngenes <- target_df$ngenes
          tibble(mean_diff, parameter = c_name, w.in_cred, w.in_prior, prop_gt, diff_sc, var_sc, mean_post)
        }) %>% 
        bind_cols(. , target_df)
    })
  
  val_df <- val_df %>% 
    mutate(is_zero = ifelse(sfs1_mean == 0, "zero", "not_zero"))

  val_df
}



viz_val <- function(val_df){
  bns <- 30
  pgt <- val_df %>% 
    ggplot(aes(x = prop_gt, y = parameter, fill = parameter, colour =  parameter)) +
    geom_density_ridges(stat = "binline", bins = bns, scale = 1, rel_min_height = 0.01) +
    geom_vline(xintercept = 0.5, lty = 2) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme(legend.position = "n") 
  
  cred <- val_df %>% 
    ggplot(aes(x = w.in_cred, y = parameter, fill = parameter, colour =  parameter)) +
    geom_density_ridges(stat = "binline", bins = 2, scale = 0.5, rel_min_height = 0.01) +
    theme(axis.text.y=element_blank(),
          legend.position = "n") +
    scale_x_continuous(breaks = c(0,1)) +
    ylab("")
  
  wpr <- val_df %>% 
    ggplot(aes(x = w.in_prior, y = parameter, fill = parameter, colour =  parameter)) +
    geom_density_ridges(stat = "binline", bins = 2, scale = 0.5, rel_min_height = 0.01) +
    theme(axis.text.y=element_blank(),
          legend.position = "n") +
    scale_x_continuous(breaks = c(0,1)) +
    ylab("")
  
  
  df_sc <- val_df %>% 
    ggplot(aes(x = diff_sc, y = parameter, fill = parameter, colour =  parameter)) +
    geom_density_ridges(stat = "binline", bins = bns, scale = 1, rel_min_height = 0.01) +
    geom_vline(xintercept = 0, lty = 2) +
    theme(axis.text.y=element_blank(),
          legend.position = "n") +
    ylab("")
  
  v_sc <- 
    val_df %>%
    filter(sfs1_mean != 0) %>% 
    ggplot(aes(x = log(var_sc), y = parameter, fill = parameter, colour =  parameter)) +
    geom_density_ridges(stat = "binline", bins = bns, scale = 1, rel_min_height = 0.001) +
    geom_vline(xintercept = log(1), lty = 2, lwd =1.2) +
    theme(axis.text.y=element_blank(),
          legend.position = "n") +
    ylab("")
  
  pgt + cred + wpr + v_sc + plot_layout(nrow = 1) 
}

val_df08 <- abc_val("NAM_DFE_ABC/pop_sfs_fixedMu_1e-08.txt", "NAM_DFE_ABC/pop_fixedMu_1e-08.txt")
val_df09 <- abc_val("NAM_DFE_ABC/pop_sfs_fixedMu_1e-09.txt", "NAM_DFE_ABC/pop_fixedMu_1e-09.txt")
val_df10 <- abc_val("NAM_DFE_ABC/pop_sfs_fixedMu_1e-10.txt", "NAM_DFE_ABC/pop_fixedMu_1e-10.txt")

viz_val(val_df08)
viz_val(val_df09)
viz_val(val_df10)



write_csv(x = val_df09, path = "val.csv")
val_df <- read_csv("val.csv")

val_df %>% 
  mutate(
    parameter = str_replace(parameter, "sfs1_shape", "shape0"),
    parameter = str_replace(parameter, "sfs2_shape", "shapeSV"),
    parameter = str_replace(parameter, "sfs1_mean", "s0"),
    parameter = str_replace(parameter, "sfs2_mean", "sSV"),
    ) %>% 
viz_val() + ggsave("NAM_DFE_ABC/figures/abc_validate.pdf")

val_df %>% 
  filter(parameter %in% c("sfs1_mean", "sfs2_mean")) %>%
  mutate(var = mean_post) %>%
  ggplot(aes(parameter, -var)) +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.05) +
  facet_grid(~is_zero)

val_df %>% 
  filter(parameter %in% c("sfs1_mean", "sfs2_mean"), sfs2_mean == 0, sfs1_mean == 0) %>%
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
  filter(parameter %in% c("mu_sv")) %>%
  mutate(var = mean_post) %>%
  ggplot(aes(mu_sv, var, colour = is_zero)) +
  geom_point() +
  geom_smooth()


val_df %>% 
  filter(parameter %in% c("sfs2_mean"), sfs2_mean == 0) %>% 
  lm(-mean_post ~ N0 + Nb + mu_sv + sfs2_shape + sfs1_shape, data = .) %>% 
  summary()

val_df %>% 
  filter(parameter %in% c("sfs2_mean"), sfs2_mean != 0) %>% 
  lm(-mean_post ~ N0 + Nb + mu_sv + sfs2_shape + sfs1_shape, data = .) %>% 
  summary()


