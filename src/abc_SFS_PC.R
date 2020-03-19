library(abc)
library(tidyverse)
par(ask=F)

# require(abc.data)
# data(human)
# nrow(stat.3pops.sim)

var0 <- function(x) var(x) > 0

full_df <- read_delim(
  "~/Dropbox/NAM_DFE_ABC/pop_sfs.txt", delim = " ",
  col_names = FALSE
) %>%
  select(- ncol(.)) %>% 
  select_if(var0) %>% 
  rowwise() %>% 
  mutate(sm  = sum(.)) %>%
  ungroup() %>% 
  mutate_all(.funs = ~ .x / sm) %>% 
  select(-sm)

pcx <- prcomp(full_df, scale. = F, center = F)

#plot(pcx)
#plot(pcx$x[,1], pcx$x[,2])
varz <- cumsum((pcx$sdev ^ 2) / sum(pcx$sdev ^ 2)) <= 1
mx_col <- max(c(2,length(varz[varz])))
pc_df <- data.frame(pcx$x) %>% select(1:mx_col)


params <- 
  c("Na", "N0", "Nb", "B_t", 
    "sfs1_shape", "sfs1_mean", "sfs2_shape", "sfs2_mean", 
    "p_neutral", "p_sfs1")

param_df <- read_delim(
  "~/Dropbox/NAM_DFE_ABC/pop.txt",
  col_names = params,
  delim = " "
) %>% 
  select(-c(Na, B_t))


target_idx <- sample(seq_len(nrow(full_df)), 1)

# target_stats <- pc_df[target_idx, ]
# sumstat_df <- pc_df[-target_idx, ]

target_stats <- full_df[target_idx, ]
sumstat_df <- full_df[-target_idx, ]

target_df <- param_df[target_idx, ]
param_df <-  param_df[-target_idx, ]

# res.gfit.bott <- gfit(target=target_stats, 
#                       sumstat=sumstat_df,
#                       statistic=mean, nb.replicate=10)
# summary(res.gfit.bott)$pvalue


dim(target_stats)
dim(param_df)
dim(sumstat_df)

res <- abc(target=target_stats,
           param=param_df,
           sumstat=sumstat_df, 
           tol=0.01, 
           transf=c("log"), 
           method = "neuralnet", 
           sizenet = 5)

# res <- abc(target=pc_target,
#            param=param_df,
#            sumstat=pc_df,
#            tol=0.05, transf=c("log"), method="neuralnet")

posts_df <- data.frame(res$adj.values)
names(posts_df)
c_name <- "N0"

den_prior <- density(param_df[[c_name]])
den_post <- density(posts_df[[c_name]])


plot(den_prior, 
     xlim = range(c(range(den_prior$x), range(den_post$x))),
     ylim = range(c(range(den_prior$y), range(den_post$y))),
     col = "red"
)
abline(v = mean(param_df[[c_name]]), col = "red")

lines(den_post, col = "blue")
abline(v = target_df[[c_name]])
abline(v = mean(posts_df[[c_name]]), col = "blue")
diff(quantile(param_df[[c_name]], c(0.025, 0.975))) > diff(quantile(posts_df[[c_name]], c(0.025, 0.975)))

quantile(param_df[[c_name]], c(0.025, 0.975))
quantile(posts_df[[c_name]], c(0.025, 0.975))

mean(posts_df[[c_name]] > target_df[[c_name]])
posts_df <- data.frame(res$adj.values)
mean(posts_df[[c_name]] > target_df[[c_name]])

target_df

# i <- sample(1:nrow(param_df), 1)
# shp <- param_df$sfs1_shape[i]
# #shp <- 2
# mn <- param_df$sfs1_mean[i]
# dst <- rgamma(
#   n = 100, 
#   shape = shp, 
#   scale = mn/shp)
# hist(
#   dst
# )
# sd(dst)
# mean(dst)
# abline(v = param_df$sfs1_mean[i])
# 
# shp <- 20; mn <- 1e-4
# skrr <- rgamma(n = 100, shape = shp, scale = mn)
# hist(skrr)
# abline(v = shp*mn, col = "blue")
# abline(v = mean(skrr))
# 
# 
# 
# shp <- runif(1, 2, 100)
# mn <- runif(1, 0, 0.01)
# dst <- rgamma(
#   n = 100, 
#   shape = shp, 
#   scale = mn/shp)
# hist(
#   -dst
# )
# sd(dst)
# mean(dst)
# abline(v = param_df$sfs1_mean[i])
