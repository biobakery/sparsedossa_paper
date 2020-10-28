rm(list = ls())
library(magrittr)
library(ggplot2)

load("results/DA/PA_results.RData")

### power analysis

tb_bench <- 
  tidyr::expand_grid(
    n = seq(100, 1000, by = 100),
    effect_size = c(0.5, 1, 2),
    R = seq(1, 500)
  )

p_power <- tb_bench %>% 
  dplyr::mutate(p = results,
                effect_size = factor(effect_size,
                                     levels = c(0.5, 1, 2))) %>% 
  dplyr::group_by(n, effect_size) %>% 
  dplyr::summarise(power = mean(p < 0.05),
                   sd = sd((p < 0.05) * 1) / sqrt(500)) %>% 
  dplyr::ungroup() %>% 
  ggplot(aes(x = n, y = power, color = effect_size)) +
  geom_point(size = 2) + 
  geom_line() +
  geom_errorbar(aes(ymin = power - sd, ymax = power + sd),
                width = 5) +
  scale_color_discrete(labels = paste0(c(0.5, 1, 2)),
                       name = "Log FC") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_blank())

### reproducible
load("data/mouse/physeq_mouse.RData")
physeq_mouse <- physeq_mouse_filtered %>% 
  phyloseq::transform_sample_counts(function(x) x / sum(x))
tb_filter <- tibble::tibble(
  day = c(0, 1, 1, 5, 5),
  diet = c("Chow", "Meat", "Tuber", "Meat", "Tuber")
)

load("results/Mouse/fitted.RData")
l_sim <- list()
l_df_sim <- list()
for(i_filter in seq_len(nrow(tb_filter))) {
  i_day <- tb_filter$day[i_filter]
  i_diet <- tb_filter$diet[i_filter]
  i_physeq <- physeq_mouse %>% 
    phyloseq::subset_samples(Day.as.listed.in.MS.SI == i_day &
                               diet_overall == i_diet)
  n_sample <- i_physeq %>% phyloseq::nsamples()
  i_df <- i_physeq %>% smar::sample_data2()
  l_df_sim[[i_filter]] <- rbind(i_df, i_df, i_df, i_df, i_df)
  
  l_sim[[i_filter]] <- SparseDOSSA2::SparseDOSSA2(
    l_fit[[i_filter]],
    n_sample = n_sample * 5,
    new_features = FALSE,
    spike_metadata = "none"
  )$simulated_matrices$rel
}

mat_sim <- l_sim %>% purrr::reduce(cbind)
df_sim <- l_df_sim %>% purrr::reduce(rbind)

dist_sim <- vegan::vegdist(t(mat_sim))
ord <- ape::pcoa(dist_sim)
perc <- (ord$values$Relative_eig[1:2] * 100) %>% round(digits = 2)

df_plot <- data.frame(`Axis 1` = ord$vectors[, 1],
                      `Axis 2` = ord$vectors[, 2],
                      check.names = FALSE) %>% 
  cbind(df_sim) %>% 
  dplyr::mutate(shape = ifelse(Day.as.listed.in.MS.SI == 1,
                               "circle",
                               "solid"),
                Date = Day.as.listed.in.MS.SI %>% 
                  as.character() %>% 
                  dplyr::recode_factor(
                    "0" = "Baseline",
                    "1" = "Day 1",
                    "5" = "Day 5"
                  ),
                Diet = diet_overall %>% 
                  factor(levels = c("Chow", "Tuber", "Meat")))

p_mouse <- df_plot %>% 
  ggplot(aes(x = `Axis 1`, y = `Axis 2`, 
             color = Diet)) +
  geom_point(aes(shape = Date), size = 3) +
  scale_color_manual(values = c("Chow" = "grey",
                                "Meat" = "purple",
                                "Tuber" = "brown")) +
  scale_shape_manual(values = c("Baseline" = 16,
                                "Day 1" = 1,
                                "Day 5" = 16)) +
  theme_bw() +
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.background = element_blank(),
        legend.direction = "horizontal") +
  xlab(paste0("Axis 1 (", perc[1], "%)")) +
  ylab(paste0("Axis 2 (", perc[2], "%)"))

cowplot::plot_grid(p_power, p_mouse, labels = c("A", "B"),
                   rel_widths = c(1, 1), nrow = 1) %>% 
  ggsave("figures/Fig4/Fig4.pdf", ., width = 10, height = 5)
