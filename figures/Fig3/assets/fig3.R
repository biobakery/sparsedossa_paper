rm(list = ls())
library(magrittr)
library(ggplot2)
## Stool
# order feature by abundance
load("../sparsedossa_update/SparseDOSSA2/R/sysdata.rda")
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
  template = "Stool",
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
                                "Spiked features",
                                "Null features") %>% 
                  factor(levels = c("Spiked features", "Null features")),
                sig = ifelse(q < 0.05,
                             "Bonf. p < 0.05",
                             "Bonf. p >= 0.05") %>% 
                  factor(levels = c("Bonf. p < 0.05", "Bonf. p >= 0.05"))) %>% 
  dplyr::arrange(spiked, 
                 -abs(estimate)) %>% 
  dplyr::mutate(feature = feature %>% forcats::as_factor()) %>% 
  ggplot(aes(x = feature,
             y = estimate,
             color = spiked,
             alpha = sig)) +
  geom_point(aes(size = spiked)) +
  geom_errorbar(aes(ymax = estimate + 1.96 * se,
                    ymin = estimate - 1.96 *se)) +
  scale_color_manual(values = c("Spiked features" = "red",
                                "Null features" = "black")) +
  scale_alpha_manual(values = c("Bonf. p < 0.05" = 1, "Bonf. p >= 0.05" = 0.2)) +
  scale_size_manual(values = c("Spiked features" = 2, "Null features" = 0.5)) +
  scale_y_continuous(breaks = c(-1, 0, 1), 
                     limits = c(-2, 2),
                     oob = scales::squish) +
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
  xlab("Feature") + ylab("Log FC (non-zero abundance)")
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
  template = "Stool",
  n_sample = 1000, 
  new_features = FALSE,
  spike_metadata = "prevalence",
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
                                "Spiked features",
                                "Null features") %>% 
                  factor(levels = c("Spiked features", "Null features")),
                sig = ifelse(q < 0.05,
                             "Bonf. p < 0.05",
                             "Bonf. p >= 0.05") %>% 
                  factor(levels = c("Bonf. p < 0.05", "Bonf. p >= 0.05"))) %>% 
  dplyr::arrange(spiked, 
                 -abs(estimate)) %>% 
  dplyr::mutate(feature = feature %>% forcats::as_factor()) %>% 
  ggplot(aes(x = feature,
             y = estimate,
             color = spiked,
             alpha = sig)) +
  geom_point(aes(size = spiked)) +
  geom_errorbar(aes(ymax = estimate + 1.96 * se,
                    ymin = estimate - 1.96 *se)) +
  scale_color_manual(values = c("Spiked features" = "red",
                                "Null features" = "black")) +
  scale_alpha_manual(values = c("Bonf. p < 0.05" = 1, "Bonf. p >= 0.05" = 0.2)) +
  scale_size_manual(values = c("Spiked features" = 2, "Null features" = 0.5)) +
  scale_y_continuous(breaks = c(-1, 0, 1), 
                     limits = c(-2, 2),
                     oob = scales::squish) +
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


# Vaginal results for supp figure 3
load("data/physeqs/vaginal.RData")
features_all <- phyloseq::taxa_names(physeq_vaginal_bl)
## p abundance
features_ordered <- physeq_vaginal_bl %>%
  phyloseq::otu_table() %>%
  apply(2, function(x) x > 0) %>%
  apply(1, mean) %>%
  {order(-.)} %>%
  {features_all[.]}
features_toSpike <- seq(1, 6,
                        length.out = 6) %>%
  round() %>% {features_ordered[.]}
metadata <- matrix(rep(c(0, 1), each = 500),
                   nrow = 1000,
                   ncol = 1)
effect_sizes <- rep(c(1, -1), 3)
df_spike <- data.frame(
  metadata_datum = 1,
  feature_spiked = features_toSpike,
  associated_property = "abundance",
  effect_size = effect_sizes
)
set.seed(1)
fit_simulation <- SparseDOSSA2::SparseDOSSA2(
  template = "Vaginal",
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
p_abd_vag <- df_effectSize %>%
  dplyr::mutate(spiked = ifelse(feature %in% features_toSpike,
                                "Spiked features",
                                "Null features") %>%
                  factor(levels = c("Spiked features", "Null features")),
                sig = ifelse(q < 0.05,
                             "Bonf. p < 0.05",
                             "Bonf. p >= 0.05") %>%
                  factor(levels = c("Bonf. p < 0.05", "Bonf. p >= 0.05"))) %>%
  dplyr::arrange(spiked,
                 -abs(estimate)) %>%
  dplyr::mutate(feature = feature %>% forcats::as_factor()) %>%
  ggplot(aes(x = feature,
             y = estimate,
             color = spiked,
             alpha = sig)) +
  geom_point(aes(size = spiked)) +
  geom_errorbar(aes(ymax = estimate + 1.96 * se,
                    ymin = estimate - 1.96 *se)) +
  scale_color_manual(values = c("Spiked features" = "red",
                                "Null features" = "black")) +
  scale_alpha_manual(values = c("Bonf. p < 0.05" = 1, "Bonf. p >= 0.05" = 0.2)) +
  scale_size_manual(values = c("Spiked features" = 2, "Null features" = 0.5)) +
  scale_y_continuous(breaks = c(-1, 0, 1),
                     limits = c(-2, 2),
                     oob = scales::squish) +
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
  xlab("Feature") + ylab("Log FC\n(non-zero abundance)")

## p prevalence
features_ordered <- physeq_vaginal_bl %>%
  phyloseq::otu_table() %>%
  apply(2, function(x) x > 0) %>%
  apply(1, mean) %>%
  {order(abs(. - 0.5))} %>%
  {features_all[.]}
features_toSpike <- seq(1, 6,
                        length.out = 6) %>%
  round() %>% {features_ordered[.]}
metadata <- matrix(rep(c(0, 1), each = 500),
                   nrow = 1000,
                   ncol = 1)
effect_sizes <- rep(c(1, -1), 3)
df_spike <- data.frame(
  metadata_datum = 1,
  feature_spiked = features_toSpike,
  associated_property = "prevalence",
  effect_size = effect_sizes
)
set.seed(1)
fit_simulation <- SparseDOSSA2::SparseDOSSA2(
  template = "Vaginal",
  n_sample = 1000,
  new_features = FALSE,
  spike_metadata = "prevalence",
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
p_prev_vag <- df_effectSize %>%
  dplyr::mutate(spiked = ifelse(feature %in% features_toSpike,
                                "Spiked features",
                                "Null features") %>%
                  factor(levels = c("Spiked features", "Null features")),
                sig = ifelse(q < 0.05,
                             "Bonf. p < 0.05",
                             "Bonf. p >= 0.05") %>%
                  factor(levels = c("Bonf. p < 0.05", "Bonf. p >= 0.05"))) %>%
  dplyr::arrange(spiked,
                 -abs(estimate)) %>%
  dplyr::mutate(feature = feature %>% forcats::as_factor()) %>%
  ggplot(aes(x = feature,
             y = estimate,
             color = spiked,
             alpha = sig)) +
  geom_point(aes(size = spiked)) +
  geom_errorbar(aes(ymax = estimate + 1.96 * se,
                    ymin = estimate - 1.96 *se)) +
  scale_color_manual(values = c("Spiked features" = "red",
                                "Null features" = "black")) +
  scale_alpha_manual(values = c("Bonf. p < 0.05" = 1, "Bonf. p >= 0.05" = 0.2)) +
  scale_size_manual(values = c("Spiked features" = 2, "Null features" = 0.5)) +
  scale_y_continuous(breaks = c(-1, 0, 1),
                     limits = c(-2, 2),
                     oob = scales::squish) +
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
cowplot::plot_grid(p_abd_vag, p_prev_vag,
                   ncol = 1,
                   rel_heights = c(1, 1),
                   align = "hv") %>% 
  ggsave("Supp_Materials/SuppFigures/SuppFig3.pdf", 
         ., width = 7, height = 4)

### feature-feature spike-in
set.seed(1)
load("results/fitted/SparseDOSSA2_full/fit_1.RData")
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
                    tidyr::pivot_longer(cols = dplyr::one_of(features_sub),
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
                    tidyr::pivot_longer(cols = dplyr::one_of(features_sub),
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

tb_plot <- rbind(tb_1, tb_2) %>% 
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
  dplyr::mutate(Assoc. = dplyr::case_when(
    is.na(Cor) ~ " (pos)",
    Cor > 0 ~ " (pos)",
    TRUE ~ " (neg)"
  ) %>% 
    paste0(Assoc., .) %>% 
    factor(levels = c("TP (pos)", "TP (neg)",
                      "Null (pos)", "Null (neg)")))
p_cor <- tb_plot %>% 
  ggplot(aes(x = Feature1, y = Feature2,
             fill = Cor, color = `Assoc.`)) +
  geom_tile(aes(fill = Cor),
            width = 0.9, height = 0.9, size = 1) +
  theme_bw() +
  scale_fill_gradient2(breaks = seq(-1, 1, by = 0.5),
                       limits = c(-1, 1),
                       low = "blue", mid = "white", high = "red",
                       midpoint = 0, na.value = "grey") +
  scale_color_manual(values = c("TP (pos)" = "black",
                                "TP (neg)" = "black",
                                "Null (pos)" = "white",
                                "Null (neg)" = "white"),
                     labels = c("TP", "", "Null", "")) +
  # scale_size_manual(values = c("TP" = 1,
  #                               "Null" = 0)) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  facet_grid(.~case) +
  xlab("Feature1") + ylab("Feature2") +
  ggtitle("Relative abundance Spearman (top left) vs.\nAbsolute abundance Spearman (bottom right)") +
  guides(color = guide_legend(override.aes = list(fill = c("red", "blue", "lightpink", "thistle2"))))
cowplot::plot_grid(p_abd, p_prev, p_cor, labels = c("A", "B", "C"),
                   ncol = 1,
                   rel_heights = c(1, 1, 1.5)) %>% 
  ggsave("figures/Fig3/Fig3.pdf", 
         ., width = 7, height = 10)

# More correlation spike-in cases for SuppFig 4
tb_2$case <- "Spiked, Effect Size = 5"
df_spike <- 
  rbind(data.frame(
    metadata_datum = c(1, 1, 2, 2),
    feature_spiked = features_toSpike,
    associated_property = "abundance",
    effect_size = c(1, 1, 1, -1)),
    data.frame(
      metadata_datum = c(1, 1, 2, 2),
      feature_spiked = features_toSpike,
      associated_property = "prevalence",
      effect_size = c(1, 1, 1, -1)))
fit_simulation <- SparseDOSSA2::SparseDOSSA2(
  template = fit_EM,
  n_sample = 100000, 
  new_features = FALSE,
  spike_metadata = "both",
  feature_metadata_spike_df = df_spike,
  metadata_matrix = metadata
)
cor_a <- fit_simulation$simulated_matrices$a_spiked[features_sub, ] %>%
  t() %>% cor(method = "sp")
cor_rel <- fit_simulation$simulated_matrices$rel[features_sub, ] %>%
  t() %>% cor(method = "sp")
tb_3 <- list(cor_a, cor_rel) %>%
  purrr::map2_dfr(c("Abs", "Rel"),
                  function(cor_mat, class)
                    cor_mat %>%
                    as.data.frame() %>%
                    tibble::rownames_to_column("Feature1") %>%
                    tidyr::pivot_longer(cols = dplyr::one_of(features_sub),
                                        names_to = "Feature2",
                                        values_to = "Cor") %>%
                    dplyr::mutate(Feature1 = factor(Feature1, levels = features_sub),
                                  Feature2 = factor(Feature2, levels = features_sub)) %>%
                    dplyr::filter(class == "Abs" & as.integer(Feature1) > as.integer(Feature2) |
                                    class == "Rel" & as.integer(Feature1) < as.integer(Feature2))) %>%
  rbind(tibble::tibble(Feature1 = features_sub,
                       Feature2 = features_sub,
                       Cor = NA)) %>% 
  dplyr::mutate(case = "Spiked, Effect Size = 1")
df_spike <- 
  rbind(data.frame(
    metadata_datum = c(1, 1, 2, 2),
    feature_spiked = features_toSpike,
    associated_property = "abundance",
    effect_size = c(2, 2, 2, -2)),
    data.frame(
      metadata_datum = c(1, 1, 2, 2),
      feature_spiked = features_toSpike,
      associated_property = "prevalence",
      effect_size = c(2, 2, 2, -2)))
fit_simulation <- SparseDOSSA2::SparseDOSSA2(
  template = fit_EM,
  n_sample = 100000, 
  new_features = FALSE,
  spike_metadata = "both",
  feature_metadata_spike_df = df_spike,
  metadata_matrix = metadata
)
cor_a <- fit_simulation$simulated_matrices$a_spiked[features_sub, ] %>%
  t() %>% cor(method = "sp")
cor_rel <- fit_simulation$simulated_matrices$rel[features_sub, ] %>%
  t() %>% cor(method = "sp")
tb_4 <- list(cor_a, cor_rel) %>%
  purrr::map2_dfr(c("Abs", "Rel"),
                  function(cor_mat, class)
                    cor_mat %>%
                    as.data.frame() %>%
                    tibble::rownames_to_column("Feature1") %>%
                    tidyr::pivot_longer(cols = dplyr::one_of(features_sub),
                                        names_to = "Feature2",
                                        values_to = "Cor") %>%
                    dplyr::mutate(Feature1 = factor(Feature1, levels = features_sub),
                                  Feature2 = factor(Feature2, levels = features_sub)) %>%
                    dplyr::filter(class == "Abs" & as.integer(Feature1) > as.integer(Feature2) |
                                    class == "Rel" & as.integer(Feature1) < as.integer(Feature2))) %>%
  rbind(tibble::tibble(Feature1 = features_sub,
                       Feature2 = features_sub,
                       Cor = NA)) %>% 
  dplyr::mutate(case = "Spiked, Effect Size = 2")

tb_plot1 <- rbind(tb_1, tb_3, tb_4, tb_2) %>% 
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
      factor(levels = c("TP", "Null")))

# vaginal
set.seed(2)
load("data/physeqs/vaginal.RData")
load("results/fitted/SparseDOSSA2_full/fit_8.RData")
features_all <- phyloseq::taxa_names(physeq_vaginal_bl)
features_ordered <- physeq_vaginal_bl %>% 
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
                    tidyr::pivot_longer(cols = dplyr::one_of(features_sub),
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
                    tidyr::pivot_longer(cols = dplyr::one_of(features_sub),
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

# More correlation spike-in cases for SuppFig 4
tb_2$case <- "Spiked, Effect Size = 5"
df_spike <- 
  rbind(data.frame(
    metadata_datum = c(1, 1, 2, 2),
    feature_spiked = features_toSpike,
    associated_property = "abundance",
    effect_size = c(1, 1, 1, -1)),
    data.frame(
      metadata_datum = c(1, 1, 2, 2),
      feature_spiked = features_toSpike,
      associated_property = "prevalence",
      effect_size = c(1, 1, 1, -1)))
fit_simulation <- SparseDOSSA2::SparseDOSSA2(
  template = fit_EM,
  n_sample = 100000, 
  new_features = FALSE,
  spike_metadata = "both",
  feature_metadata_spike_df = df_spike,
  metadata_matrix = metadata
)
cor_a <- fit_simulation$simulated_matrices$a_spiked[features_sub, ] %>%
  t() %>% cor(method = "sp")
cor_rel <- fit_simulation$simulated_matrices$rel[features_sub, ] %>%
  t() %>% cor(method = "sp")
tb_3 <- list(cor_a, cor_rel) %>%
  purrr::map2_dfr(c("Abs", "Rel"),
                  function(cor_mat, class)
                    cor_mat %>%
                    as.data.frame() %>%
                    tibble::rownames_to_column("Feature1") %>%
                    tidyr::pivot_longer(cols = dplyr::one_of(features_sub),
                                        names_to = "Feature2",
                                        values_to = "Cor") %>%
                    dplyr::mutate(Feature1 = factor(Feature1, levels = features_sub),
                                  Feature2 = factor(Feature2, levels = features_sub)) %>%
                    dplyr::filter(class == "Abs" & as.integer(Feature1) > as.integer(Feature2) |
                                    class == "Rel" & as.integer(Feature1) < as.integer(Feature2))) %>%
  rbind(tibble::tibble(Feature1 = features_sub,
                       Feature2 = features_sub,
                       Cor = NA)) %>% 
  dplyr::mutate(case = "Spiked, Effect Size = 1")
df_spike <- 
  rbind(data.frame(
    metadata_datum = c(1, 1, 2, 2),
    feature_spiked = features_toSpike,
    associated_property = "abundance",
    effect_size = c(2, 2, 2, -2)),
    data.frame(
      metadata_datum = c(1, 1, 2, 2),
      feature_spiked = features_toSpike,
      associated_property = "prevalence",
      effect_size = c(2, 2, 2, -2)))
fit_simulation <- SparseDOSSA2::SparseDOSSA2(
  template = fit_EM,
  n_sample = 100000, 
  new_features = FALSE,
  spike_metadata = "both",
  feature_metadata_spike_df = df_spike,
  metadata_matrix = metadata
)
cor_a <- fit_simulation$simulated_matrices$a_spiked[features_sub, ] %>%
  t() %>% cor(method = "sp")
cor_rel <- fit_simulation$simulated_matrices$rel[features_sub, ] %>%
  t() %>% cor(method = "sp")
tb_4 <- list(cor_a, cor_rel) %>%
  purrr::map2_dfr(c("Abs", "Rel"),
                  function(cor_mat, class)
                    cor_mat %>%
                    as.data.frame() %>%
                    tibble::rownames_to_column("Feature1") %>%
                    tidyr::pivot_longer(cols = dplyr::one_of(features_sub),
                                        names_to = "Feature2",
                                        values_to = "Cor") %>%
                    dplyr::mutate(Feature1 = factor(Feature1, levels = features_sub),
                                  Feature2 = factor(Feature2, levels = features_sub)) %>%
                    dplyr::filter(class == "Abs" & as.integer(Feature1) > as.integer(Feature2) |
                                    class == "Rel" & as.integer(Feature1) < as.integer(Feature2))) %>%
  rbind(tibble::tibble(Feature1 = features_sub,
                       Feature2 = features_sub,
                       Cor = NA)) %>% 
  dplyr::mutate(case = "Spiked, Effect Size = 2")

tb_plot2 <- rbind(tb_1, tb_3, tb_4, tb_2) %>% 
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
      factor(levels = c("TP", "Null")))

tb_plot <- rbind(tb_plot1 %>% 
                   dplyr::mutate(dataset = "Stool"), 
                 tb_plot2 %>% 
                   dplyr::mutate(dataset = "Vaginal")) %>% 
  dplyr::mutate(Assoc. = dplyr::case_when(
    is.na(Cor) ~ " (pos)",
    Cor > 0 ~ " (pos)",
    TRUE ~ " (neg)"
  ) %>% 
    paste0(Assoc., .) %>% 
    factor(levels = c("TP (pos)", "TP (neg)",
                      "Null (pos)", "Null (neg)")))

p_cor <- tb_plot %>% 
  ggplot(aes(x = Feature1, y = Feature2, fill = Cor, color = Assoc.)) +
  geom_tile(width = 0.9, height = 0.9, size = 1) +
  theme_bw() +
  scale_fill_gradient2(breaks = seq(-1, 1, by = 0.5),
                       limits = c(-1, 1),
                       low = "blue", mid = "white", high = "red",
                       midpoint = 0, na.value = "grey") +
  scale_color_manual(values = c("TP (pos)" = "black",
                                "TP (neg)" = "black",
                                "Null (pos)" = "white",
                                "Null (neg)" = "white"),
                     labels = c("TP", "", "Null", "")) +
  # scale_size_manual(values = c("TP" = 1,
  #                               "Null" = 0)) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  facet_wrap(~dataset + case, nrow = 2, scales = "free") +
  xlab("Feature1") + ylab("Feature2") +
  ggtitle("Relative abundance Spearman (top left) vs.\nAbsolute abundance Spearman (bottom right)") +
  guides(color = guide_legend(override.aes = list(fill = c("red", "blue", "lightpink", "thistle2"))))
ggsave(filename = "Supp_Materials/SuppFigures/SuppFig4.pdf",
       p_cor, width = 12, height = 7)


# synergies
# stool
sim_stool <- SparseDOSSA2::SparseDOSSA2(template = "Stool",
                                        n_sample = 100000,
                                        new_features = FALSE,
                                        spike_metadata = "none")
mat_count <- sim_stool$simulated_data
mat_rel <- apply(mat_count, 2, function(x) x / sum(x))
cor_mat <- cor(t(mat_rel), method = "spearman")
features <- rownames(cor_mat)
phyla <- features %>% 
  stringr::str_replace(stringr::fixed("k__Bacteria|p__"), "") %>% 
  stringr::str_replace("\\|c\\_\\_.*", "")
features_firm <- features[phyla %in% c("Firmicutes")]
features_bact <- features[phyla %in% c("Bacteroidetes")]
features_firm_sel <- apply(mat_rel[features_firm, ], 1, mean) %>% 
  {order(-.)[seq(1, 5)]} %>% 
  {features_firm[.]}
features_bact_sel <- apply(mat_rel[features_bact, ], 1, mean) %>% 
  {order(-.)[seq(5, 1)]} %>% 
  {features_bact[.]}
features_firm_rename <- features_firm_sel %>% 
  stringr::str_replace("^.*\\|s\\_\\_", "") %>% 
  stringr::str_replace("\\_", " ") %>% 
  paste0(" (F)")
features_bact_rename <- features_bact_sel %>% 
  stringr::str_replace("^.*\\|s\\_\\_", "") %>% 
  stringr::str_replace("\\_", " ") %>% 
  paste0(" (B)")
df_cor <- as.data.frame(cor_mat[c(features_firm_sel, features_bact_sel), c(features_firm_sel, features_bact_sel)]) %>% 
  tibble::rownames_to_column("Feature1") %>% 
  tidyr::pivot_longer(cols = -Feature1, names_to = "Feature2", values_to = "Spearman_Cor") %>% 
  dplyr::filter(Feature1 != Feature2) %>% 
  dplyr::mutate(Feature1_rename = Feature1 %>% 
                  stringr::str_replace("^.*\\|s\\_\\_", "") %>% 
                  stringr::str_replace("\\_", " "),
                Feature2_rename = Feature2 %>% 
                  stringr::str_replace("^.*\\|s\\_\\_", "") %>% 
                  stringr::str_replace("\\_", " ")) %>% 
  dplyr::mutate(Feature1_rename = 
                  ifelse(Feature1 %in% features_firm_sel, 
                         paste0(Feature1_rename, " (F)"),
                         paste0(Feature1_rename, " (B)")),
                Feature2_rename = 
                  ifelse(Feature2 %in% features_firm_sel, 
                         paste0(Feature2_rename, " (F)"),
                         paste0(Feature2_rename, " (B)"))) %>% 
  dplyr::mutate(
    Feature1_rename = 
      factor(Feature1_rename, 
             levels = c(features_firm_rename, features_bact_rename)),
    Feature2_rename = 
      factor(Feature2_rename, 
             levels = c(features_firm_rename, features_bact_rename))) 

p1 <- df_cor %>% 
  ggplot(aes(x = Feature1_rename,
             y = Feature2_rename,
             fill = Spearman_Cor)) +
  geom_tile() +
  scale_fill_gradient2(
    limits = c(-0.2, 0.4),
    breaks = seq(-0.2, 0.4, by = 0.1),
    low = "blue", mid = "white", high = "red",
    midpoint = 0, na.value = "grey", name =
      "Cor") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust=1),
        axis.title = element_blank()) +
  ggtitle("Simulated Stool Spearman Correlation")

sim_vaginal <- SparseDOSSA2::SparseDOSSA2(template = "Vaginal", 
                                          n_sample = 100000, 
                                          new_features = FALSE, 
                                          spike_metadata = "none")
mat_count <- sim_vaginal$simulated_data
mat_rel <- apply(mat_count, 2, function(x) x / sum(x))
cor_mat <- cor(t(mat_rel), method = "spearman")
dimnames(cor_mat) <- list(rownames(mat_rel), rownames(mat_rel))
features <- rownames(cor_mat)
features_top <- apply(mat_rel, 1, mean) %>% 
  {order(-.)[1:8]} %>% 
  {rownames(mat_rel)[.]}
features_top_rename <- features_top %>% 
  stringr::str_replace("^.*\\|s\\_\\_", "") %>% 
  stringr::str_replace("\\_", " ")
df_cor <- as.data.frame(cor_mat[features_top, features_top]) %>% 
  tibble::rownames_to_column("Feature1") %>% 
  tidyr::pivot_longer(cols = -Feature1, names_to = "Feature2", values_to = "Spearman_Cor") %>% 
  dplyr::filter(Feature1 != Feature2) %>%
  dplyr::mutate(Feature1_rename = Feature1 %>% 
                  stringr::str_replace("^.*\\|s\\_\\_", "") %>% 
                  stringr::str_replace("\\_", " "),
                Feature2_rename = Feature2 %>% 
                  stringr::str_replace("^.*\\|s\\_\\_", "") %>% 
                  stringr::str_replace("\\_", " ")) %>% 
  dplyr::mutate(
    Feature1_rename = 
      factor(Feature1_rename, 
             levels = features_top_rename),
    Feature2_rename = 
      factor(Feature2_rename, 
             levels = features_top_rename))  

p2 <- df_cor %>% 
  ggplot(aes(x = Feature1_rename,
             y = Feature2_rename,
             fill = Spearman_Cor)) +
  geom_tile() +
  scale_fill_gradient2(breaks = seq(-0.3, 0.1, length.out = 5),
                       limits = c(-0.3, 0.1),
                       labels = round(seq(-0.3, 0.1, by = 0.1), digits = 1),
                       low = "blue", mid = "white", high = "red",
                       midpoint = 0, na.value = "grey", name =
                         "Cor") +
  theme_bw() +
  ggtitle("Simulated Vaginal Spearman Correlation") +
  theme(axis.text.x = element_text(angle = 30, hjust=1),
        axis.title = element_blank())

p <- cowplot::plot_grid(p1, p2, labels = c("A", "B"), nrow = 1,
                        align = "h")
ggsave(p, filename = "Supp_Materials/Response_Letter/fig_synergy.pdf",
       width = 15,
       height = 6)

mat_rel_clr <- apply(mat_rel, 2, function(x) {
  x <- x + min(setdiff(x, 0)) / 2
  log(x) - mean(log(x))
})
cor_mat <- sapply(
  seq(1, nrow(mat_rel_clr)),
  function(i) {
    sapply(
      seq(1, nrow(mat_rel_clr)),
      function(j) {
        var(mat_rel_clr[i, ] - mat_rel_clr[j, ]) /
          var(mat_rel_clr[i, ])
      })
  })
dimnames(cor_mat) <- list(rownames(mat_rel), rownames(mat_rel))
features <- rownames(cor_mat)
features_top <- apply(mat_rel, 1, mean) %>% 
  {order(-.)[1:8]} %>% 
  {rownames(mat_rel)[.]}
features_top_rename <- features_top %>% 
  stringr::str_replace("^.*\\|s\\_\\_", "") %>% 
  stringr::str_replace("\\_", " ")
df_cor <- as.data.frame(cor_mat[features_top, features_top]) %>% 
  tibble::rownames_to_column("Feature1") %>% 
  tidyr::pivot_longer(cols = -Feature1, names_to = "Feature2", values_to = "Spearman_Cor") %>% 
  dplyr::filter(Feature1 != Feature2) %>%
  dplyr::mutate(Feature1_rename = Feature1 %>% 
                  stringr::str_replace("^.*\\|s\\_\\_", "") %>% 
                  stringr::str_replace("\\_", " "),
                Feature2_rename = Feature2 %>% 
                  stringr::str_replace("^.*\\|s\\_\\_", "") %>% 
                  stringr::str_replace("\\_", " ")) %>% 
  dplyr::mutate(
    Feature1_rename = 
      factor(Feature1_rename, 
             levels = features_top_rename),
    Feature2_rename = 
      factor(Feature2_rename, 
             levels = features_top_rename))  

p21 <- df_cor %>% 
  ggplot(aes(x = Feature1_rename,
             y = Feature2_rename,
             fill = Spearman_Cor)) +
  geom_tile() +
  scale_fill_gradient(breaks = seq(0, 12, by = 2),
                       limits = c(0, 12),
                       low = "white", high = "red",
                       na.value = "grey", name =
                         "Cor") +
  theme_bw() +
  ggtitle("Simulated Vaginal Proportionality Coefficient") +
  theme(axis.text.x = element_text(angle = 30, hjust=1),
        axis.title = element_blank())


cor_mat <- 2 * cov(t(mat_rel_clr)) / 
  outer(apply(mat_rel_clr, 1, var), apply(mat_rel_clr, 1, var), "+")
dimnames(cor_mat) <- list(rownames(mat_rel), rownames(mat_rel))
features <- rownames(cor_mat)
features_top <- apply(mat_rel, 1, mean) %>% 
  {order(-.)[1:8]} %>% 
  {rownames(mat_rel)[.]}
features_top_rename <- features_top %>% 
  stringr::str_replace("^.*\\|s\\_\\_", "") %>% 
  stringr::str_replace("\\_", " ")
df_cor <- as.data.frame(cor_mat[features_top, features_top]) %>% 
  tibble::rownames_to_column("Feature1") %>% 
  tidyr::pivot_longer(cols = -Feature1, names_to = "Feature2", values_to = "Spearman_Cor") %>% 
  dplyr::filter(Feature1 != Feature2) %>%
  dplyr::mutate(Feature1_rename = Feature1 %>% 
                  stringr::str_replace("^.*\\|s\\_\\_", "") %>% 
                  stringr::str_replace("\\_", " "),
                Feature2_rename = Feature2 %>% 
                  stringr::str_replace("^.*\\|s\\_\\_", "") %>% 
                  stringr::str_replace("\\_", " ")) %>% 
  dplyr::mutate(
    Feature1_rename = 
      factor(Feature1_rename, 
             levels = features_top_rename),
    Feature2_rename = 
      factor(Feature2_rename, 
             levels = features_top_rename))  

p22 <- df_cor %>%
  ggplot(aes(x = Feature1_rename,
             y = Feature2_rename,
             fill = Spearman_Cor)) +
  geom_tile() +
  scale_fill_gradient2(breaks = seq(-0.1, 0.1, by = 0.05),
                      limits = c(-0.1, 0.1),
                      low = "blue", mid = "white", high = "red",
                      midpoint = 0, 
                      na.value = "grey", name =
                        "Cor") +
  theme_bw() +
  ggtitle("Simulated Vaginal Concordance Correlation Coefficient") +
  theme(axis.text.x = element_text(angle = 30, hjust=1),
        axis.title = element_blank())

p <- cowplot::plot_grid(p21, p22, nrow = 1,
                        align = "h")
ggsave(p, filename = "Supp_Materials/Response_Letter/fig_synergy_alt.jpeg",
       width = 15,
       height = 6)
