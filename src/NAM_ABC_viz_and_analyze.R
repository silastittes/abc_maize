library(tidyverse)
library(patchwork)
library(cowplot)
theme_set(theme_cowplot())

params <- 
  c("Na", "N0", "Nb", "B_t", "mu_sv",
    "sfs1_shape", "sfs1_mean", "sfs2_shape", "sfs2_mean")

param_df <- read_delim(
  "NAM_DFE_ABC/pop.txt",
  col_names = params,
  delim = " "
) %>% 
  dplyr::select(-c(Na, B_t))


all_abc_post <-
  bind_rows(
    read_delim(file = "NAM_DFE_ABC/NAM_abc_out/NAM_ABC_SFS-dup_inv_tra_ins_del_df.csv", delim = ",")
  )

win_df <- read_delim(
  "NAM_DFE_ABC/NAM_abc_out/window_stats_20Mb.txt", delim = "\t", skip = 1,
  col_names = c("chrom", "start", "end", "fold0", "fold4", "cds_bps", "genes")) %>% 
  mutate(window = 1:n(),
         length = end - start)
  

# sfs <- read_delim("~/Dropbox/NAM_DFE_ABC/abc_out/NAM_SFS_INS_DEL-gte2Kb_all.txt", "\t", col_names = FALSE)
# sums <- apply(sfs, 1, sum)

all_abc_post %>% 
  group_by(window) %>% 
  summarise()

unique(all_abc_post$sfs)

xx <- seq(0, quantile(all_abc_post$sfs2_mean, 0.99), length.out = 500)
DFE2_vec <- sample(1:nrow(all_abc_post), 1) %>% 
  map(~{
    dgamma(x = xx, shape = all_abc_post$sfs2_shape[.x], scale = all_abc_post$sfs2_mean[.x]/all_abc_post$sfs2_shape[.x])
  })
mxx <- max(map_dbl(DFE2_vec, max))
plot(xx, xx, type = 'n', ylim=c(0, 1))
walk(DFE2_vec, ~{lines(xx, .x/max(.x), col = alpha("black", 0.9))})

xx <- seq(0, quantile(all_abc_post$sfs1_mean, 0.99), length.out = 500)
DFE2_vec <- sample(1:nrow(all_abc_post), 1) %>% 
  map(~{
    dgamma(x = xx, shape = all_abc_post$sfs1_shape[.x], scale = all_abc_post$sfs1_mean[.x]/all_abc_post$sfs1_shape[.x])
  })
mxx <- max(map_dbl(DFE2_vec, max))
#plot(xx, xx, type = 'n', ylim=c(0, mxx))
walk(DFE2_vec, ~{lines(xx, .x/max(.x), col = alpha("blue", 0.9))})



all_abc_post %>%
  ggplot(aes(x = sfs1_mean, group = window)) +
  geom_density(colour=alpha("black", 0.05)) +
  geom_density(data = param_df, aes(x = sfs1_mean), inherit.aes = F, colour = "red") +
  facet_wrap(~sfs)
 
all_abc_post %>%
  ggplot(aes(x = sfs2_mean, group = window)) +
  geom_density(colour=alpha("black", 0.05)) +
  geom_density(data = param_df, aes(x = sfs2_mean), inherit.aes = F, colour = "red") +
  facet_wrap(~sfs)

all_abc_post %>%
  ggplot(aes(x = log10(N0), group = window)) +
  geom_density(colour=alpha("black", 0.05)) +
  geom_density(data = param_df, aes(x = log10(N0)), inherit.aes = F, colour = "red")

all_abc_post %>%
  ggplot(aes(x = N0, group = window)) +
  geom_density(colour=alpha("black", 0.05)) +
  geom_density(data = param_df, aes(x = N0), inherit.aes = F, colour = "red") 


 all_abc_post %>%
  ggplot(aes(x = Nb, group = window)) +
  geom_density(colour=alpha("black", 0.05)) +
  geom_density(data = param_df, aes(x = Nb), inherit.aes = F, colour = "red")
 
 all_abc_post %>%
   ggplot(aes(x = log10(mu_sv), group = window)) +
   geom_density(colour=alpha("black", 0.05)) +
   geom_density(data = param_df, aes(x = log10(mu_sv)), inherit.aes = F, colour = "red")
   
hist(log(all_abc_post$mu_sv) )

post_window_mean <-
  all_abc_post %>%
  group_by(window, sfs) %>%
  mutate(pr_dfe2gtdfe1 = sfs1_mean < sfs2_mean) %>%
  summarise_all(.funs = mean) %>% 
  mutate(idx = 1:n()) %>% 
  full_join(., win_df, by = "window")


post_window_mean %>%
  mutate(dfe_0fold = sfs1_mean*1000, 
         dfe_sv = sfs2_mean*1000) %>% 
  dplyr::select(sfs, dfe_0fold, dfe_sv) %>% 
  pivot_longer(cols = c(dfe_0fold, dfe_sv), names_to = "dfe", values_to = "Ns") %>% 
  ggplot(aes(Ns, dfe)) +
  ggridges::geom_density_ridges(alpha = 0.5, rel_min_height = 0.001) +
  ggsave(filename = "~/Desktop/dfe.png")

post_window_mean %>%
  mutate(dfe_0fold = sfs1_mean, 
         dfe_sv = sfs2_mean) %>% 
  dplyr::select(sfs, dfe_0fold, dfe_sv) %>% 
  pivot_longer(cols = c(dfe_0fold, dfe_sv), names_to = "dfe", values_to = "s") %>% 
  ggplot(aes(dfe, -s)) +
  ylab("mean posterior selection coefficient") +
  xlab("Variant type") +
  scale_x_discrete(labels = c('0-fold','SV')) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.4) +
  ggsave("NAM_DFE_ABC/figures/dfe.pdf")


posts_df %>%
  mutate(dfe_0fold = sfs1_mean, 
         dfe_sv = sfs2_mean) %>% 
  dplyr::select(dfe_0fold, dfe_sv) %>% 
  pivot_longer(cols = c(dfe_0fold, dfe_sv), names_to = "dfe", values_to = "Ns") %>% 
  ggplot(aes(Ns, dfe)) +
  ggridges::geom_density_ridges(alpha = 0.5, rel_min_height = 0.001)

post_window_mean %>%
  mutate(dfe_0fold = sfs1_mean, 
         dfe_sv = sfs2_mean) %>% 
  dplyr::select(sfs, dfe_0fold, dfe_sv) %>% 
  pivot_longer(cols = c(dfe_0fold, dfe_sv), names_to = "dfe", values_to = "s") %>% 
  ggplot(aes(-s, dfe)) +
  ggridges::geom_density_ridges(alpha = 0.5, rel_min_height = 0.001) +
  ggsave(filename = "~/Desktop/dfe.png")

post_window_mean %>%
  ggplot(aes(sfs1_mean, sfs2_mean, colour = sfs)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

post_window_mean %>%
  ggplot(aes(window, sfs2_mean, colour = chrom)) +
  geom_point() +
  ylab("s") +
  ggsave(filename = "~/Desktop/s_window.png")

post_window_mean %>%
  ggplot(aes(window, sfs1_mean, colour = chrom)) +
  geom_point() +
  ggsave(filename = "~/Desktop/s_window.png")

post_window_mean %>%
  ggplot(aes(genes, -sfs1_mean)) +
  geom_point() +
  geom_smooth(se = T) + 
  ylab("0 fold s") +
  xlab("number of genes") +
  ggsave(filename = "~/Desktop/fold0_gene_window.png")

post_window_mean %>%
  ggplot(aes(genes, -sfs2_mean)) +
  geom_point() +
  geom_smooth(se = T) +
  ylab("SV s") +
  xlab("number of genes") +
  ggsave(filename = "~/Desktop/SV_gene_window.png")

summary(lm(sfs2_mean ~ genes, data = post_window_mean))
summary(lm(sfs1_mean ~ genes, data = post_window_mean))

  
length(unique(post_window_mean$window))
mean(post_window_mean$pr_dfe2gtdfe1 >= 0.95)
hist(post_window_mean$pr_dfe2gtdfe1, breaks = 100)
1/mean(1/post_window_mean$pr_dfe2gtdfe1)
exp(mean(log(post_window_mean$pr_dfe2gtdfe1)))
mean(post_window_mean$pr_dfe2gtdfe1)
hist(log10(post_window_mean$sfs1_mean))
hist(log10(post_window_mean$sfs2_mean))
hist(post_window_mean$sfs2_mean)
mean(post_window_mean$sfs2_mean)
mean(post_window_mean$sfs1_mean)
plot(post_window_mean$sfs1_mean, post_window_mean$sfs2_mean)
abline(0,1)
mean(post_window_mean$sfs1_mean < post_window_mean$sfs2_mean)
plot(post_window_mean$Nb, post_window_mean$N0)
hist(log10(post_window_mean$N0), breaks = 100)
hist(log10(post_window_mean$Nb), breaks = 100)
plot(log10(post_window_mean$N0), log10(post_window_mean$Nb))

hist(log10(post_window_mean$N0), breaks = 100)
hist(log10(post_window_mean$Nb), breaks = 100)
plot(post_window_mean$N0, post_window_mean$Nb)
abline(0,1)

plot(post_window_mean$sfs2_mean)
plot(post_window_mean$sfs1_mean)

plot(post_window_mean$sfs2_mean, post_window_mean$sfs1_mean)
abline(0,1)

plot(post_window_mean$sfs2_mean, post_window_mean$mu_sv)

