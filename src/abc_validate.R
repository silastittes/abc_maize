library(abc)
library(tidyverse)
library(ggridges)
library(patchwork)
library(cowplot)
theme_set(theme_cowplot())
par(ask=F)

n_ind <- 26
fold_bins <- (1:(n_ind/2))
to_keep <- c(fold_bins, fold_bins + n_ind*1, fold_bins + n_ind*2)
#var0 <- function(x) var(x) > 0

full_df <- read_delim(
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


itts <- 100
p <- progress_estimated(itts)

set.seed(123)
val_df <- 
  sample(seq_len(nrow(full_df)), itts) %>% 
  map_df(function(target_idx){
    print(p$tick())
    #PC
    # target_stats <- pc_df[target_idx, ]
    # sumstat_df <- pc_df[-target_idx, ]

    #raw
    target_stats <- full_df[target_idx, ]
    sumstat_df <- full_df[-target_idx, ]
    
    target_df <- param_df_full[target_idx, ]
    param_df <-  param_df_full[-target_idx, ]
    
    
    res <- abc(target=target_stats,
               param=param_df,
               sumstat=sumstat_df,
               tol=0.005,
               transf=c("none"),
               method = "neuralnet",
               sizenet = 2)
    
    posts_df <- data.frame(res$adj.values)
    
    names(posts_df) %>% 
      map_df(function(c_name){
        prop_gt <- mean(posts_df[[c_name]] > target_df[[c_name]])
        creds <- quantile(posts_df[[c_name]], c(0,0.025,0.975,1))
        w.in_cred <- as.numeric(creds[2] < target_df[[c_name]] & creds[3] > target_df[[c_name]])
        w.in_prior <- as.numeric(creds[1] < target_df[[c_name]] & creds[4] > target_df[[c_name]])
        diff_sc <- mean( (posts_df[[c_name]] - target_df[[c_name]]) / target_df[[c_name]] )
        var_sc <- sd(posts_df[[c_name]] / sd(param_df[[c_name]]))
        tibble(parameter = c_name, w.in_cred, w.in_prior, prop_gt, diff_sc, var_sc)
      })
  })


bns <- 30
pgt <- val_df %>% 
  ggplot(aes(x = prop_gt, y = parameter, fill = parameter, colour =  parameter)) +
  geom_density_ridges(stat = "binline", bins = bns, scale = 1, rel_min_height = 0.01) +
  geom_vline(xintercept = 0.5, lty = 2) +
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
  ggplot(aes(x = log(var_sc), y = parameter, fill = parameter, colour =  parameter)) +
  geom_density_ridges(stat = "binline", bins = bns, scale = 1, rel_min_height = 0.001) +
  geom_vline(xintercept = log(1), lty = 2) +
  theme(axis.text.y=element_blank(),
        legend.position = "n") +
  ylab("")


pgt + cred + wpr + v_sc + plot_layout(nrow = 1) +
  ggsave("~/Desktop/abc_val.png")

hist(param_df_full$sfs1_shape)
mean(param_df_full$sfs1_shape)
mean(param_df_full$sfs2_shape)

target_idx <- sample(seq_len(nrow(full_df)), 1)
#PC
target_stats <- pc_df[target_idx, ]
sumstat_df <- pc_df[-target_idx, ]

#raw
target_stats <- full_df[target_idx, ]
sumstat_df <- full_df[-target_idx, ]

target_df <- param_df_full[target_idx, ]
param_df <-  param_df_full[-target_idx, ]


#target_stats[,(39-12):(39)] <- 0
#fake <- c(1919,1254,742,192,441,333,238,542,469,797,508,136,30,0,0,0,0,0,0,0,0,0,0,0,0,0,18,11,16,10,24,8,4,6,14,1,11,5,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
#target_stats <- fake[to_keep] / sum(fake)
res <- abc(target=target_stats,
           param=param_df,
           sumstat=sumstat_df,
           tol=0.005,
           transf=c("none"),
           method = "neuralnet",
           sizenet = 2)

# res <- abc(target=target_stats,
#            param=param_df,
#            sumstat=sumstat_df,
#            tol=0.1,
#            transf=c("none"),
#            method = "loclinear",
#            sizenet = 5)

posts_df <- data.frame(res$adj.values)
c_name <- "sfs2_mean"

den_prior <- density(param_df[[c_name]])
den_post <- density(posts_df[[c_name]])

plot(posts_df[["sfs1_mean"]], posts_df[["sfs2_mean"]])
mean(posts_df[["sfs1_mean"]] > posts_df[["sfs2_mean"]])
abline(0,1)

mean(posts_df[["N0"]] > posts_df[["Nb"]])
plot(posts_df[["N0"]], posts_df[["Nb"]])
abline(0,1)

plot(den_prior, 
     xlim = range(c(range(den_prior$x), range(den_post$x))),
     ylim = range(c(range(den_prior$y), range(den_post$y))),
     col = "red"
)
abline(v = mean(param_df[[c_name]]), col = "red")

lines(den_post, col = "blue")
abline(v = target_df[[c_name]])
abline(v = mean(posts_df[[c_name]]), col = "blue")

log(sd(posts_df[[c_name]] / sd(param_df[[c_name]])))

diff(quantile(param_df[[c_name]], c(0.025, 0.975))) > diff(quantile(posts_df[[c_name]], c(0.025, 0.975)))

quantile(param_df[[c_name]], c(0.025, 0.975))
quantile(posts_df[[c_name]], c(0.025, 0.975))

mean(posts_df[[c_name]] > target_df[[c_name]])
posts_df <- data.frame(res$adj.values)
mean(posts_df[[c_name]] > target_df[[c_name]])

target_df

