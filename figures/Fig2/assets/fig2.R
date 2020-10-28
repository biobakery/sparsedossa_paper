rm(list = ls())
library(magrittr)
library(ggplot2)

load("results/Evaluate/results_R2.RData")
levels_dataset <- c("Stool",
                    "Vaginal",
                    "IBD")
levels_method <- c("DMM",
                   "metaSPARSim",
                   "SparseDOSSA 2")
tb_results_R2 <- 
  tb_results_R2 %>% 
  dplyr::mutate(
    dataset = dataset %>% 
      dplyr::recode(
        "vaginal" = "Vaginal",
        "stool" = "Stool"
      ) %>% 
      factor(levels = levels_dataset),
    method = method %>% 
      dplyr::recode(
        "SparseDOSSA2" = "SparseDOSSA 2"
      ) %>% 
      factor(levels = levels_method)
  )

# panel A

tb_A <- tb_results_R2 %>% 
  dplyr::filter(R == 1)
set.seed(1)
# ps_A <- c("stool",
#           "vaginal",
#           "IBD") %>% 
#   purrr::map2(
#     levels_dataset,
#     function(dataset, name_dataset) {
#       # original data
#       load(paste0("data/physeqs/", dataset, ".RData"))
#       original_data <- get(paste0("physeq_", dataset, "_bl")) %>% 
#         smar::otu_table2()
#       
#       # simulated data
#       name_data <- tb_A %>% 
#         dplyr::filter(method == "SparseDOSSA 2",
#                       dataset == name_dataset) %>% 
#         extract2("i_job") %>% 
#         paste0("results/Simulation/", ., ".RData")
#       load(name_data)
#       
#       dataset_combined <- cbind(original_data, sim_data) %>% 
#         apply(2, function(x) x / sum(x))
#       D <- vegan::vegdist(t(dataset_combined), method = "bray")
#       fit_pcoa <- ape::pcoa(D)
#       tb_plot <- 
#         tibble::tibble(sample = rep(c("Real-world",
#                                       "SparseDOSSA 2"),
#                                     each = ncol(original_data)) %>% 
#                          forcats::as_factor()) %>% 
#         dplyr::mutate(PCo1 = fit_pcoa$vectors[, 1],
#                       PCo2 = fit_pcoa$vectors[, 2])
#       perc_char <- fit_pcoa$values$Relative_eig[1:2] %>% 
#         {. * 100} %>% 
#         round(digits = 2) %>% 
#         paste0("(", ., "%)")
#       
#       p <- tb_plot %>% 
#         dplyr::slice_sample(prop = 1) %>% 
#         ggplot(aes(x = PCo1,
#                    y = PCo2,
#                    color = sample)) +
#         geom_point() +
#         scale_color_manual(values = c("Real-world" = "black",
#                                       "SparseDOSSA 2" = "red"),
#                            name = NULL) +
#         theme_bw() +
#         ggtitle(name_dataset) +
#         xlab(paste0("PCo1 ", perc_char[1])) +
#         ylab(paste0("PCo2 ", perc_char[2])) +
#         coord_fixed()
#       
#       if(name_dataset == "Stool")
#         p <- p + theme(legend.position = c(1, 1),
#                        legend.justification = c(1, 1),
#                        legend.box.background = element_rect(colour = "black"),
#                        legend.title=element_blank(),
#                        legend.margin=margin(c(0,0,0,0)))
#       else 
#         p <- p + theme(legend.position = "none")
#       
#       return(p)
#     }
#   )
# pA <- ps_A %>% 
#   cowplot::plot_grid(plotlist = ., 
#                      ncol = 1, 
#                      rel_heights = c(1, 1, 1), 
#                      align = "v")

pA <- c("stool",
        "vaginal",
        "IBD") %>%
  purrr::map2_dfr(
    levels_dataset,
    function(dataset, name_dataset) {
      # original data
      load(paste0("data/physeqs/", dataset, ".RData"))
      original_data <- get(paste0("physeq_", dataset, "_bl")) %>%
        smar::otu_table2()
      
      # simulated data
      name_data <- tb_A %>%
        dplyr::filter(method == "SparseDOSSA 2",
                      dataset == name_dataset) %>%
        extract2("i_job") %>%
        paste0("results/Simulation/", ., ".RData")
      load(name_data)
      
      sim_data1 <- sim_data
      
      name_data <- tb_A %>%
        dplyr::filter(method == "metaSPARSim",
                      dataset == name_dataset) %>%
        extract2("i_job") %>%
        paste0("results/Simulation/", ., ".RData")
      load(name_data)
      
      dataset_combined <- cbind(original_data, sim_data) %>%
        apply(2, function(x) x / sum(x))
      D <- vegan::vegdist(t(dataset_combined), method = "bray")
      fit_pcoa <- ape::pcoa(D)
      tb_plot <-
        tibble::tibble(sample = rep(c("Real-world",
                                      "SparseDOSSA 2"),
                                    each = ncol(original_data)) %>%
                         forcats::as_factor(),
                       dataset = name_dataset) %>%
        dplyr::mutate(PCo1 = fit_pcoa$vectors[, 1],
                      PCo2 = fit_pcoa$vectors[, 2])
      
    }) %>% 
  dplyr::mutate(dataset = factor(dataset, levels = levels_dataset)) %>% 
  ggplot(aes(x = PCo1, y = PCo2, color = sample)) +
  geom_point() +
  scale_color_manual(values = c("Real-world" = "black",
                                "SparseDOSSA 2" = "red")) +
  theme_bw() +
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.box.background = element_rect(colour = "black"),
        legend.title=element_blank(),
        legend.margin=margin(c(0,0,0,0)),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  facet_wrap(.~dataset, ncol = 1, scales = "free",
             strip.position = "left") +
  ggtitle("Ordination")


# panel B
# ps_B <- levels_dataset %>% 
#   purrr::map(function(name_dataset) {
#     p <- tb_results_R2 %>% 
#       dplyr::filter(dataset == name_dataset) %>% 
#       ggplot(aes(x = method, y = R2)) +
#       geom_boxplot() +
#       ggtitle("") +
#       theme_bw()
#     
#     if(name_dataset == "Stool")
#       p <- p + 
#         ylab("") +
#         xlab("") +
#         theme(axis.text.x = element_blank())
#     if(name_dataset == "Vaginal")
#       p <- p + 
#         ylab("R2 (real-world vs. simulated)") +
#         xlab("") +
#         theme(axis.text.x = element_blank())
#     if(name_dataset == "IBD")
#       p <- p + 
#         ylab("") +
#         xlab("Fitted model") +
#         theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))
#     return(p)
#   })
# 
# pB <- ps_B %>% 
#   cowplot::plot_grid(plotlist = ., 
#                      ncol = 1, 
#                      rel_heights = c(1, 1, 1), 
#                      align = "h")
pB <- tb_results_R2 %>% 
  ggplot(aes(x = method, y = R2)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(dataset~., scales = "free_y") +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 15, vjust = 1, hjust=1),
        axis.title.x = element_blank()) +
  ggtitle(expression(paste(R^{2}, " (real-world vs. simulated)"))) +
  ylab(expression(R^{2}))

cowplot::plot_grid(pA, pB, labels = c("A", "B"),
                   align = "h", axis = "tb")

# panel C
load("results/Evaluate/results_ks.RData")
pC <- tb_results_ks %>% 
  dplyr::mutate(
    dataset = dataset %>% 
      dplyr::recode("stool" = "Stool",
                    "vaginal" = "Vaginal") %>% 
      factor(levels = levels_dataset),
    method = method %>% 
      dplyr::recode("SparseDOSSA2" = "SparseDOSSA 2") %>% 
      factor(levels = levels_method)
  ) %>% 
  tidyr::pivot_wider(id_cols = c(dataset, feature),
                     names_from = method,
                     values_from = ks) %>% 
  tidyr::pivot_longer(cols = c(DMM, metaSPARSim),
                      names_to = "method",
                      values_to = "Other methods") %>% 
  ggplot(aes(x = `SparseDOSSA 2`, y = `Other methods`, color = method)) +
  geom_point() +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  facet_grid(dataset ~ .) +
  theme(legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.background = element_rect(colour = "black"),
        legend.title=element_blank(),
        legend.margin=margin(c(0,0,0,0)),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  ggtitle("Per-feature K-S statistic")

cowplot::plot_grid(pA, pB, pC, labels = c("A", "B", "C"),
                   align = "h", axis = "tb", nrow = 1,
                   rel_widths = c(0.9, 1.2, 1)) %>% 
  ggsave("figures/Fig2/Fig2.pdf",
       ., width = 8, height = 6)
 