rm(list = ls())
library(magrittr)
library(ggplot2)

# panel A
load("results/all_comparison/tb_summary_null.RData")
load("results/all_comparison/tb_summary_nonNull.RData")

tb_summary_null$effect <- 0

p_FDR <- rbind(tb_summary_null,
               tb_summary_nonNull[, colnames(tb_summary_null)]) %>% 
  dplyr::filter(Method %in% c("limmaVOOM", "ANCOM", "MaAsLin2")) %>% 
  ggplot(aes(x = factor(effect),
             y = FDR,
             color = Method)) +
  geom_point(position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = FDR - FDR_sd,
                    ymax = FDR + FDR_sd),
                width = 0.25,
                position = position_dodge(width = 0.25)) +
  theme_bw() +
  theme(legend.position = c(0.99, 0.99),
        legend.justification = c(1, 1),
        legend.background = element_blank()) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  xlab("Effect size") +
  ggtitle("Benchmarking analysis")
p_Recall <- tb_summary_nonNull %>% 
  dplyr::filter(Method %in% c("limmaVOOM", "ANCOM", "MaAsLin2")) %>% 
  ggplot(aes(x = factor(effect),
             y = Recall,
             color = Method)) +
  geom_point(position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = Recall - Recall_sd,
                    ymax = Recall + Recall_sd),
                width = 0.25,
                position = position_dodge(width = 0.25)) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Effect size") +
  ggtitle("")

### power analysis
load("results/Power_analysis/tb_summary.RData")
p_power <- tb_results_summary %>% 
  dplyr::mutate(effect = as.character(effect)) %>% 
  ggplot(aes(x = n, y = power, color = effect)) +
  geom_point(size = 2, aes(shape = effect)) + 
  geom_line(aes(linetype = effect)) +
  geom_errorbar(aes(ymin = power - sd_power / sqrt(500), 
                    ymax = power + sd_power / sqrt(500)),
                width = 5) +
  scale_color_manual(labels = paste0(c(0.5, 1, 2)),
                       name = "Log FC",
                       values = c("0.5" = "lightgrey",
                                  "1" = "darkgrey",
                                  "2" = "black")) +
  scale_shape_discrete(labels = paste0(c(0.5, 1, 2)),
                       name = "Log FC") +
  scale_linetype_discrete(labels = paste0(c(0.5, 1, 2)),
                       name = "Log FC") +
  scale_x_continuous(breaks = seq(200, 1000, by = 200)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.99, 0.01),
        legend.justification = c(1, 0),
        legend.background = element_blank()) +
  ggtitle("Power analysis on E. coli")

cowplot::plot_grid(p_FDR, p_Recall, p_power, nrow = 1,
                   labels = c("A", "", "B"),
                   rel_widths = c(1, 0.75, 1)) %>% 
  ggsave("figures/Fig4/Fig4.pdf", ., width = 10, height = 4)