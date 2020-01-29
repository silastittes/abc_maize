library(abc)
library(tidyverse)
par(ask=F)

#set.seed(123)
con <- file("~/Dropbox/NAM_DFE_ABC/pop.txt","r")
first_line <- readLines(con,n=1)
close(con)
col_string <- str_split(first_line, " ")[[1]] %>% str_remove(pattern = "#")

full_df <- read_delim(
  "~/Dropbox/NAM_DFE_ABC/pop.txt", delim = " ",
  col_names = col_string, comment = "#")

pc_sumstat_df <- prcomp(
  select(full_df, 
         -c(mu, Na, N0, Nb, B_t, sfs1_shape, sfs1_mean, sfs2_shape, sfs2_mean, p_neutral, p_sfs1)
         ), 
  center = T, scale. = T)
plot(pc_sumstat_df$x[,1], pc_sumstat_df$x[,2])
(pc_sumstat_df$sdev^2/sum(pc_sumstat_df$sdev^2))

hist(full_df$sfs_pi1)
abline(v = 0.008)
summary(full_df$sfs_pi1)
hist(full_df$sfs_td1)
range(full_df$sfs_td1)
#pairs(select(full_df, -c(mu, Na, N0, Nb, B_t, sfs1_shape, sfs1_mean, sfs2_shape, sfs2_mean, p_neutral, p_sfs1)))
cor(full_df)

target_idx <- sample(seq_len(nrow(full_df)), 1)
target_df <- full_df[target_idx, ]
full_df <- full_df[-target_idx, ]

sumstat_df <- select(full_df, -c(mu, Na, N0, Nb, B_t, sfs1_shape, sfs1_mean, sfs2_shape, sfs2_mean, p_neutral, p_sfs1))
param_df <- select(full_df, c(mu, Na, N0, Nb, B_t, sfs1_shape, sfs1_mean, sfs2_shape, sfs2_mean, p_neutral, p_sfs1))
target_stats <- select(target_df, -c(mu, Na, N0, Nb, B_t, sfs1_shape, sfs1_mean, sfs2_shape, sfs2_mean, p_neutral, p_sfs1))
target_params <- select(target_df, c(mu, Na, N0, Nb, B_t, sfs1_shape, sfs1_mean, sfs2_shape, sfs2_mean, p_neutral, p_sfs1))


res.gfit.bott <- gfit(target=target_stats, 
                   sumstat=sumstat_df,
                   statistic=mean, nb.replicate=10)
summary(res.gfit.bott)$pvalue

res <- abc(target=target_stats,
           param=param_df,
           sumstat=sumstat_df, 
           tol=0.05, transf=c("log"), method="ridge")


posts_df <- data.frame(res$unadj.values)
names(posts_df)
c_name <- "N0"
plot(density(posts_df[[c_name]]))

posts_df[[c_name]]

den_prior <- density(param_df[[c_name]])
den_post <- density(posts_df[[c_name]])
plot(den_prior, 
     xlim = range(c(range(den_prior$x), range(den_post$x))),
     ylim = range(c(range(den_prior$y), range(den_post$y)))
     )

lines(den_post, col = "blue")
abline(v = target_params[[c_name]])
diff(quantile(param_df[[c_name]], c(0.025, 0.975))) > diff(quantile(posts_df[[c_name]], c(0.025, 0.975)))

mean(posts_df[[c_name]] > target_params[[c_name]])
posts_df <- data.frame(res$adj.values)
mean(posts_df[[c_name]] > target_params[[c_name]])

