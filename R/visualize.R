# visualize
df_plot <- data.frame(read_depth = n_i, is_zero = y_i == 0)

# Biological vs. sequencing zeros
# Create zero (biological vs. sequencing) probabilities
df_plot_zeroDecomp <- df_plot %>% 
  dplyr::mutate(`Sequencing zero` = i_fit$hidden_param$w_i,
                `Biological zero` = 1 - `Sequencing zero`) 
# Formatting data
df_plot_zeroDecomp <- df_plot_zeroDecomp %>% 
  dplyr::filter(is_zero) %>% 
  tidyr::gather(key = "Type of zero",
                value = "Probability",
                `Sequencing zero`, `Biological zero`)
df_plot_zeroDecomp <- df_plot_zeroDecomp %>% 
  dplyr::filter(is_zero) %>% 
  dplyr::mutate(
    `Type of zero` = factor(`Type of zero`, 
                            levels = c("Sequencing zero", 
                                       "Biological zero"))) %>% 
  # to prevent duplicate read_depth from affecting plotting
  dplyr::group_by(read_depth, `Type of zero`) %>% 
  dplyr::mutate(read_depth_factor = 
                  read_depth + (0:(dplyr::n() - 1))*0.01) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(read_depth_factor) %>% 
  dplyr::mutate(read_depth_factor = 
                  factor(read_depth_factor,
                         levels = unique(read_depth_factor)))
p_zeroDecomp <- df_plot_zeroDecomp %>% 
  dplyr::filter(`Type of zero` == "Biological zero") %>% 
  ggplot(aes(x = read_depth,
             y = Probability)) +
  scale_x_continuous(trans = "log10") +
  geom_point() +
  xlab("Effective read depth") + ylab("P(0 is biological)")

# sample vs. posterior relative abundance estimates
df_plot_ra <- rbind(
  df_plot %>% 
    dplyr::mutate(`RA` = i_fit$hidden_param$mu_sample_i,
                  `sd` = sqrt(i_fit$hidden_param$sigma2_sample_i),
                  Estimation = "Sample"),
  df_plot %>% 
    dplyr::mutate(`RA` = i_fit$hidden_param$mu_posterior_i,
                  `sd` = sqrt(i_fit$hidden_param$sigma2_posterior_i),
                  Estimation = "Posterior")
)
df_plot_ra <- df_plot_ra %>% 
  dplyr::filter(!is_zero) %>% 
  dplyr::mutate(Estimation = factor(Estimation, 
                                    levels = c("Sample", "Posterior"))) %>% 
  # to prevent duplicate read_depth from affecting plotting
  dplyr::group_by(read_depth, Estimation) %>% 
  dplyr::mutate(read_depth_factor = 
                  read_depth + (0:(dplyr::n() - 1))*0.01) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(read_depth_factor) %>% 
  dplyr::mutate(read_depth_factor = 
                  factor(read_depth_factor,
                         levels = unique(read_depth_factor)))
width_reference <- (max(log10(df_plot_ra$read_depth)) - 
                      min(log10(df_plot_ra$read_depth))) / 50
p_ra <- df_plot_ra %>% 
  ggplot(aes(x = read_depth,
             y = RA,
             color = Estimation)) +
  geom_point(position = position_dodge(width = width_reference),
             size = 2) +
  geom_errorbar(aes(ymax = RA + sd,
                    ymin = RA - sd),
                width = width_reference,
                position = position_dodge(width = width_reference)) +
  scale_x_continuous(trans = "log10") +
  scale_color_manual(values = c("Sample" = "black",
                                "Posterior" = "red")) +
  geom_hline(yintercept = i_fit$theta["mu"],
             linetype = "dashed") +
  xlab("Effective read depth") + ylab("Logit(relative abundance)") +
  theme(legend.position = c(0, 0),
        legend.justification = c(0, 0),
        legend.background = element_blank())

label <- i_feature_simple %>% 
  paste0("\nPrevalence = ", round(mean(y_i != 0), digits = 4)) %>% 
  paste0("\nMean abundance = ", 
         format(exp(i_fit$theta_original["mu"]) / 
                  (1 + exp(i_fit$theta_original["mu"])),
                type = "e",
                digits = 4))
label_fake <- paste0("\n\n")
p_print <- cowplot::plot_grid(p_zeroDecomp + ggtitle(label),
                              p_ra + ggtitle(label_fake),
                              nrow = 1)
i_feature_simple %>% 
  paste0(dir_outputFig, ., ".pdf") %>% 
  ggsave(plot = p_print, width = 12, height = 7)
i_fit$feature <- i_feature_simple
return(i_fit)