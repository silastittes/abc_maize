library(abc)
library(tidyverse)
library(optparse)
par(ask=F)

option_list <-  list(
  make_option(c("-s", "--sim_sfs"), type="character", default=NULL, 
              help="Simulated SFS file. Space delimited", metavar="character"),
  make_option(c("-p", "--sim_param"), type="character", default=NULL, 
              help="Simulated file listing true parameters. Space delimited", metavar="character"),
  make_option(c("-e", "--obs_sfs"), type="character", default=NULL, 
              help="Observed sfs file, each line contains SFS counts for 4fold, 0fold, and sv categories.", metavar="character"),
  make_option(c("-c", "--sfs_count"), type="numeric", default=0, 
              help="Minimum raw count of required for window to be analyzed", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="Name of output file containing all posterior draws.", metavar="character")
)

opt_parser <-  OptionParser(option_list=option_list)
opt <-  parse_args(opt_parser)

n_ind <- 26
fold_bins <- (1:(n_ind/2))
to_keep <- c(fold_bins, fold_bins + n_ind*1, fold_bins + n_ind*2)

#-----------------------
#simulated SFS
#opt$sim_sfs <- "~/Dropbox/NAM_DFE_ABC/pop_sfs.txt" 
#"~/Dropbox/NAM_DFE_ABC/pop_sfs.txt"
full_df <- read_delim(opt$sim_sfs, delim = " ", col_names = FALSE) %>%
  dplyr::select(-ncol(.)) %>% 
  dplyr::select(all_of(to_keep)) %>% 
  mutate(sm  = rowSums(.)) %>%
  ungroup() %>% 
  rowwise() %>% 
  mutate_all(.funs = ~ .x / sm) %>% 
  dplyr::select(-sm)

#
#opt <- NULL
#opt$obs_sfs <- "~/Dropbox/NAM_DFE_ABC/NAM_all.txt"
#opt$sfs_count <- 500
sv_sum <- apply(read_delim(opt$obs_sfs, delim = "\t", col_names = FALSE)[ fold_bins + n_ind*2], 1, sum) > opt$sfs_count
window_targets <- which(sv_sum)

if (sum(sv_sum) > opt$sfs_count){
  #observed nam sfs data
  nam_df <- read_delim(
    opt$obs_sfs, delim = "\t",
    col_names = FALSE
  ) %>%
    dplyr::select(all_of(to_keep)) %>% 
    mutate(sm  = rowSums(.)) %>%
    rowwise() %>% 
    mutate_all(.funs = ~ .x / sm) %>% 
    dplyr::select(-sm) %>% 
    ungroup()
  
  params <- 
    c("Na", "N0", "Nb", "B_t", 
      "sfs1_shape", "sfs1_mean", "sfs2_shape", "sfs2_mean", 
      "p_neutral", "p_sfs1")
  
  
  #-----------------------
  #opt$sim_param <- "~/Dropbox/NAM_DFE_ABC/pop.txt"
  #~/Dropbox/NAM_DFE_ABC/pop.txt
  param_df <- read_delim(
    opt$sim_param,
    col_names = params,
    delim = " "
  ) %>% 
    dplyr::select(-c(Na, B_t))
  
  
  all_abc_post <- 
    window_targets %>% 
    map_df(function(target_idx){
      #raw
      target_stats <- nam_df[target_idx, ]
      res <- abc(target=target_stats,
                 param=param_df,
                 sumstat=full_df,
                 tol=0.05,
                 transf=c("log"),
                 method = "neuralnet",
                 sizenet = 2)
      
      posts_df <- data.frame(res$adj.values) 
      
      mean_post <- summarise_all(posts_df, .funs = median)
      prior_low <- summarise_all(param_df, .funs = min)
      prior_high <- summarise_all(param_df, .funs = max)
      post_in <- mean(mean_post > prior_low & mean_post < prior_high)
      
      if (post_in == 1) {
        posts_df %>% 
          mutate(window = target_idx)
      } else {
        tibble()
      }
    })
  
  #-----------------------
  write_csv(all_abc_post, opt$out)
} else {
  print("No windows has sufficient SV data. Writing empty tibble to file")
  write_csv(tibble(), opt$out)
}

window_targets
unique(all_abc_post$window)

# all_abc_post %>% 
#   ggplot(aes(x = log10(sfs1_mean), group = window)) +
#   geom_density(alpha = 0.1) +
#   geom_density(data = param_df, aes(x = log10(sfs1_mean)), inherit.aes = F, colour = "red")
# 
# all_abc_post %>% 
#   ggplot(aes(x = log10(sfs2_mean), group = window)) +
#   geom_density() +
#   geom_density(data = param_df, aes(x = log10(sfs2_mean)), inherit.aes = F, colour = "red")
# 
# ggplot() +
#   xlim(0, 0.01) +
#   geom_density(data = param_df, aes(x = sfs2_mean), inherit.aes = F, fill = "red") +
#   geom_density(data = all_abc_post, aes(x = sfs2_mean, group = window), colour=alpha("black", 0.05))
# 
# ggplot() +
#   xlim(0, 0.01) +
#   ylim(0, 250) +
#   geom_density(data = param_df, aes(x = sfs2_mean), inherit.aes = F, fill = "red") +
#   geom_density(data = all_abc_post, aes(x = sfs2_mean, group = window), colour=alpha("black", 0.05))
# 
# all_abc_post %>%
#   ggplot(aes(x = log(N0), group = window)) +
#   geom_density() +
#   geom_density(data = param_df, aes(x = log(N0)), inherit.aes = F, colour = "red")
# 
# all_abc_post %>% 
#   ggplot(aes(x = Nb, group = window)) +
#   geom_density() +
#   geom_density(data = param_df, aes(x = Nb), inherit.aes = F, colour = "red")
# 
# 
# all_abc_post %>% 
#   filter(window %in% 1:10) %>% 
#   ggplot(aes(x = sfs2_mean, group = window)) +
#   geom_density() +
#   geom_density(data = param_df, aes(x = sfs2_mean), inherit.aes = F, colour = "red")
# 
# all_abc_post %>% 
#   filter(window %in% 1:10) %>% 
#   ggplot(aes(x = sfs1_mean, group = window)) +
#   geom_density() +
#   geom_density(data = param_df, aes(x = sfs1_mean), inherit.aes = F, colour = "red")
# 
# 
# post_window_mean <- 
#   all_abc_post %>%
#   group_by(window) %>% 
#   mutate(pr_dfe2gtdfe1 = sfs1_mean < sfs2_mean) %>% 
#   summarise_all(.funs = mean)
# 
# length(unique(post_window_mean$window))
# mean(post_window_mean$pr_dfe2gtdfe1 > 0.95)
# hist(post_window_mean$pr_dfe2gtdfe1, breaks = 100)
# hist(log(post_window_mean$sfs1_mean, base = 10), breaks = 100)
# hist(post_window_mean$sfs2_mean, breaks = 100)
# mean(post_window_mean$sfs2_mean)
# mean(post_window_mean$sfs1_mean)
# plot(post_window_mean$sfs1_mean, post_window_mean$sfs2_mean)
# abline(0,1)
# plot(post_window_mean$Nb, post_window_mean$N0)
# hist(post_window_mean$N0, breaks = 100)
# hist(post_window_mean$Nb, breaks = 100)
# hist(post_window_mean$p_sfs1, breaks = 100)
# hist(post_window_mean$p_neutral, breaks = 100)
# plot(post_window_mean$sfs2_mean)
# names(posts_df)
# c_name <- "sfs1_shape"
# 
# den_prior <- density(param_df[[c_name]])
# den_post <- density(posts_df[[c_name]])
# plot(den_prior, 
#      xlim = range(c(range(den_prior$x), range(den_post$x))),
#      ylim = range(c(range(den_prior$y), range(den_post$y))),
#      col = "red"
# )
# abline(v = mean(param_df[[c_name]]), col = "red")
# 
# lines(den_post, col = "blue")
# abline(v = mean(posts_df[[c_name]]), col = "blue")
# diff(quantile(param_df[[c_name]], c(0.025, 0.975))) > diff(quantile(posts_df[[c_name]], c(0.025, 0.975)))
# 
# mean(posts_df[["sfs1_mean"]] < posts_df[["sfs2_mean"]])
# plot(posts_df[["sfs1_mean"]], posts_df[["sfs2_mean"]])
# abline(0,1)
# 
# mean(posts_df[["N0"]] > posts_df[["Nb"]])
# plot(posts_df[["N0"]], posts_df[["Nb"]])
# abline(0,1)
# 
# hist(posts_df[["Nb"]], breaks = 500)
# 
# 
# quantile(posts_df[["sfs1_mean"]], c(0.025, 0.975))
# quantile(posts_df[["sfs2_mean"]], c(0.025, 0.975))
# quantile(param_df[[c_name]], c(0.025, 0.975))
# quantile(posts_df[[c_name]], c(0.025, 0.975))

