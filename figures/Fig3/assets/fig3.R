rm(list = ls())
library(magrittr)
library(ggplot2)

# order feature by abundance
load("results/fitted/SparseDOSSA2_stool.RData")
load("data/physeqs/stool.RData")

features_all <- phyloseq::taxa_names(physeq_stool_bl)

## p abundance
features_ordered <- physeq_stool_bl %>% 
  phyloseq::otu_table() %>% 
  apply(2, function(x) x > 0) %>% 
  apply(1, mean) %>% 
  {order(-.)} %>% 
  {features_all[.]} 
features_toSpike <- seq(1, 16, 
                        length.out = 16) %>% 
  round() %>% {features_ordered[.]}
  

metadata <- matrix(rep(c(0, 1), each = 500),
                   nrow = 1000, 
                   ncol = 1)
effect_sizes <- rep(c(1, -1), 8)

df_spike <- data.frame(
  metadata_datum = 1,
  feature_spiked = features_toSpike,
  associated_property = "abundance",
  effect_size = effect_sizes
)

set.seed(1)
fit_simulation <- SparseDOSSA2::SparseDOSSA2(
  template = fit_SparseDOSSA2,
  n_sample = 1000, 
  new_features = FALSE,
  spike_metadata = "abundance",
  feature_metadata_spike_df = df_spike,
  metadata_matrix = metadata
)

mat_abd_simulated <- fit_simulation$simulated_data %>% 
  apply(2, SparseDOSSA2:::TSS)

df_effectSize <- features_all %>% 
  purrr::map_dfr(function(i_feature) {
    ind <- mat_abd_simulated[i_feature, ] > 0
    y <- log(mat_abd_simulated[i_feature, ind])
    x <- metadata[ind, 1]
    if(length(unique(x)) <= 1)
      return(
        data.frame(feature = i_feature,
                   estimate = NA,
                   se = NA,
                   p = NA)
      )
    fit_lm <- lm(y ~ x)
    data.frame(feature = i_feature,
               estimate = summary(fit_lm)$coef[2, 1],
               se = summary(fit_lm)$coef[2, 2],
               p = summary(fit_lm)$coef[2, 4])
  }) %>% 
  dplyr::filter(!is.na(p)) %>% 
  dplyr::mutate(q = p.adjust(p, method = "bonf"))

p_abd <- df_effectSize %>% 
  dplyr::mutate(spiked = ifelse(feature %in% features_toSpike,
                                "Spiked",
                                "Null") %>% 
                  factor(levels = c("Spiked", "Null")),
                sig = ifelse(q < 0.05,
                             "Bonf. p < 0.05",
                             "Bonf. p >= 0.05") %>% 
                  factor(levels = c("Bonf. p < 0.05", "Bonf. p >= 0.05"))) %>% 
  dplyr::arrange(spiked, 
                 abs(estimate) * 
                   ifelse(spiked == "Spiked",
                          0,
                          -1)) %>% 
  dplyr::mutate(feature = feature %>% forcats::as_factor()) %>% 
  ggplot(aes(x = feature,
             y = estimate,
             color = spiked,
             alpha = sig)) +
  geom_point(aes(size = spiked)) +
  geom_errorbar(aes(y = estimate,
                    ymax = estimate + 1.96 * se,
                    ymin = estimate - 1.96 *se)) +
  scale_color_manual(values = c("Spiked" = "red",
                                "Null" = "black")) +
  scale_alpha_manual(values = c("Bonf. p < 0.05" = 1, "Bonf. p >= 0.05" = 0.2)) +
  scale_size_manual(values = c("Spiked" = 2, "Null" = 0.5)) +
  scale_y_continuous(breaks = c(-1, 0, 1), 
                     limits = c(-2, 2),
                     oob = scales::oob_keep) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  geom_hline(yintercept = -1, color = "red", linetype = "dashed") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.direction = "horizontal", 
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.title = element_blank()) +
  xlab("Feature") + ylab("Log FC (abundance)")

## p prevalence
features_ordered <- physeq_stool_bl %>% 
  phyloseq::otu_table() %>% 
  apply(2, function(x) x > 0) %>% 
  apply(1, mean) %>% 
  {order(abs(. - 0.5))} %>% 
  {features_all[.]} 
features_toSpike <- seq(1, 16, 
                        length.out = 16) %>% 
  round() %>% {features_ordered[.]}


metadata <- matrix(rep(c(0, 1), each = 500),
                   nrow = 1000, 
                   ncol = 1)
effect_sizes <- rep(c(1, -1), 8)

df_spike <- data.frame(
  metadata_datum = 1,
  feature_spiked = features_toSpike,
  associated_property = "prevalence",
  effect_size = effect_sizes
)

set.seed(1)
fit_simulation <- SparseDOSSA2::SparseDOSSA2(
  template = fit_SparseDOSSA2,
  n_sample = 1000, 
  new_features = FALSE,
  spike_metadata = "presence",
  feature_metadata_spike_df = df_spike,
  metadata_matrix = metadata
)

mat_prev_simulated <- fit_simulation$simulated_data %>% 
  apply(2, SparseDOSSA2:::TSS)

df_effectSize <- features_all %>% 
  purrr::map_dfr(function(i_feature) {
    y <- (mat_prev_simulated[i_feature, ] > 0) * 1
    x <- metadata[, 1]
    fit_glm <- glm(y ~ x, family = "binomial")
    data.frame(feature = i_feature,
               estimate = summary(fit_glm)$coef[2, 1],
               se = summary(fit_glm)$coef[2, 2],
               p = summary(fit_glm)$coef[2, 4])
  }) %>% 
  dplyr::filter(!is.na(p)) %>% 
  dplyr::mutate(q = p.adjust(p, method = "bonf"))

  
p_prev <- df_effectSize %>% 
  dplyr::mutate(spiked = ifelse(feature %in% features_toSpike,
                                "Spiked",
                                "Null") %>% 
                  factor(levels = c("Spiked", "Null")),
                sig = ifelse(q < 0.05,
                             "Bonf. p < 0.05",
                             "Bonf. p >= 0.05") %>% 
                  factor(levels = c("Bonf. p < 0.05", "Bonf. p >= 0.05"))) %>% 
  dplyr::arrange(spiked, 
                 abs(estimate) * 
                   ifelse(spiked == "Spiked",
                          0,
                          -1)) %>% 
  dplyr::mutate(feature = feature %>% forcats::as_factor()) %>% 
  ggplot(aes(x = feature,
             y = estimate,
             color = spiked,
             alpha = sig)) +
  geom_point(aes(size = spiked)) +
  geom_errorbar(aes(y = estimate,
                    ymax = estimate + 1.96 * se,
                    ymin = estimate - 1.96 *se)) +
  scale_color_manual(values = c("Spiked" = "red",
                                "Null" = "black")) +
  scale_alpha_manual(values = c("Bonf. p < 0.05" = 1, "Bonf. p >= 0.05" = 0.2)) +
  scale_size_manual(values = c("Spiked" = 2, "Null" = 0.5)) +
  scale_y_continuous(breaks = c(-1, 0, 1), 
                     limits = c(-2, 2),
                     oob = scales::oob_keep) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  geom_hline(yintercept = -1, color = "red", linetype = "dashed") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  xlab("Feature") + ylab("Log OR (prevalence)")


### feature-feature spike-in
set.seed(1)
load("results/fitted/SparseDOSSA2/fit_3.RData")
features_ordered <- physeq_stool_bl %>% 
  phyloseq::otu_table() %>% 
  apply(2, function(x) x / sum(x)) %>% 
  apply(1, mean) %>% 
  {order(-.)} %>% 
  {features_all[.]} 
features_sub <- seq(1, 10) %>% {features_ordered[.]}

features_toSpike <- features_sub[sample.int(10, size = 4)]

metadata <- matrix(rnorm(200000),
                   nrow = 100000, 
                   ncol = 2)
df_spike <- 
  rbind(data.frame(
    metadata_datum = c(1, 1, 2, 2),
    feature_spiked = features_toSpike,
    associated_property = "abundance",
    effect_size = c(5, 5, 5, -5)),
    data.frame(
      metadata_datum = c(1, 1, 2, 2),
      feature_spiked = features_toSpike,
      associated_property = "prevalence",
      effect_size = c(5, 5, 5, -5)))

fit_simulation <- SparseDOSSA2::SparseDOSSA2(
  template = fit_EM,
  n_sample = 100000, 
  new_features = FALSE,
  spike_metadata = "both",
  feature_metadata_spike_df = df_spike,
  metadata_matrix = metadata
)

cor_a <- fit_simulation$simulated_matrices$a_null[features_sub, ] %>% 
  t() %>% cor(method = "sp")
cor_rel <- fit_simulation$simulated_matrices$a_null %>% 
  {apply(., 2, SparseDOSSA2:::TSS)[features_sub, ]} %>% 
  t() %>% cor(method = "sp")

tb_1 <- list(cor_a, cor_rel) %>% 
  purrr::map2_dfr(c("Abs", "Rel"),
                  function(cor_mat, class)
                    cor_mat %>% 
                    as.data.frame() %>% 
                    tibble::rownames_to_column("Feature1") %>% 
                    tidyr::pivot_longer(cols = dplyr::all_of(features_sub),
                                        names_to = "Feature2",
                                        values_to = "Cor") %>% 
                    dplyr::mutate(Feature1 = factor(Feature1, levels = features_sub),
                                  Feature2 = factor(Feature2, levels = features_sub)) %>% 
                    dplyr::filter(class == "Abs" & as.integer(Feature1) > as.integer(Feature2) |
                                    class == "Rel" & as.integer(Feature1) < as.integer(Feature2))) %>% 
  rbind(tibble::tibble(Feature1 = features_sub,
                       Feature2 = features_sub,
                       Cor = NA)) %>% 
  dplyr::mutate(case = "Null")


cor_a <- fit_simulation$simulated_matrices$a_spiked[features_sub, ] %>%
  t() %>% cor(method = "sp")
cor_rel <- fit_simulation$simulated_matrices$rel[features_sub, ] %>%
  t() %>% cor(method = "sp")
tb_2 <- list(cor_a, cor_rel) %>%
  purrr::map2_dfr(c("Abs", "Rel"),
                  function(cor_mat, class)
                    cor_mat %>%
                    as.data.frame() %>%
                    tibble::rownames_to_column("Feature1") %>%
                    tidyr::pivot_longer(cols = dplyr::all_of(features_sub),
                                        names_to = "Feature2",
                                        values_to = "Cor") %>%
                    dplyr::mutate(Feature1 = factor(Feature1, levels = features_sub),
                                  Feature2 = factor(Feature2, levels = features_sub)) %>%
                    dplyr::filter(class == "Abs" & as.integer(Feature1) > as.integer(Feature2) |
                                    class == "Rel" & as.integer(Feature1) < as.integer(Feature2))) %>%
  rbind(tibble::tibble(Feature1 = features_sub,
                       Feature2 = features_sub,
                       Cor = NA)) %>% 
  dplyr::mutate(case = "Spiked")

p_cor <- rbind(tb_1, tb_2) %>% 
  dplyr::mutate(case = case %>% forcats::as_factor()) %>% 
  dplyr::mutate(
    Assoc. = ifelse((Feature1 == features_toSpike[1] & Feature2 == features_toSpike[2] |
                       Feature2 == features_toSpike[1] & Feature1 == features_toSpike[2] |
                       Feature1 == features_toSpike[3] & Feature2 == features_toSpike[4] |
                       Feature1 == features_toSpike[4] & Feature2 == features_toSpike[3]) &
                      case != "Null",
                    "TP",
                    "Null"
    ) %>% 
      factor(levels = c("TP", "Null"))) %>% 
  ggplot(aes(x = Feature1, y = Feature2, fill = Cor, color = Assoc.)) +
  geom_tile(width = 0.9, height = 0.9, size = 1) +
  theme_bw() +
  scale_fill_gradient2(breaks = seq(-1, 1, length.out = 11),
                       limits = c(-1, 1),
                       low = "blue", mid = "white", high = "red",
                       midpoint = 0, na.value = "grey") +
  scale_color_manual(values = c("TP" = "black",
                                "Null" = "white")) +
  # scale_size_manual(values = c("TP" = 1,
  #                               "Null" = 0)) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  facet_grid(.~case) +
  xlab("Feature1") + ylab("Feature2") +
  ggtitle("Relative abundance Spearman (top left) vs.\nAbsolute abundance Spearman (bottom right)") +
  guides(color = guide_legend(override.aes = list(fill = "red")))

cowplot::plot_grid(p_abd, p_prev, p_cor, labels = c("A", "B", "C"),
                   ncol = 1,
                   rel_heights = c(1, 1, 1.5)) %>% 
  ggsave("figures/Fig3/Fig3.pdf", 
         ., width = 7, height = 9)
