rm(list = ls())
library(magrittr)
library(ggplot2)

# figure check time
load("results/Simulation/tb_job.RData")
tb_time <- tb_job %>%
  dplyr::group_by(i_job) %>% 
  dplyr::mutate(time = {
    load(paste0("results/Simulation/time/time_", i_job, ".RData"))
    as.numeric(time)
  },
  n = {
    load(paste0("data/physeqs/", dataset, ".RData"))
    x_samples <- 
      get(paste0("physeq_", dataset, "_bl")) %>%
      smar::otu_table2()
    round(ncol(x_samples) / 2)
  },
  p = {
    load(paste0("data/physeqs/", dataset, ".RData"))
    x_samples <- 
      get(paste0("physeq_", dataset, "_bl")) %>%
      smar::otu_table2()
    nrow(x_samples)
  }) %>% 
  dplyr::ungroup()
p_time <- tb_time %>% 
  dplyr::mutate(method = method %>%
           dplyr::recode("SparseDOSSA2" = "SparseDOSSA 2") %>% 
           factor(levels = c("DMM", "metaSPARSim", "SparseDOSSA 2")),
         dataset = dataset %>% 
           dplyr::recode("stool" = "Stool",
                         "vaginal" = "Vaginal") %>% 
           factor(levels = c("Stool", "Vaginal", "IBD"))) %>% 
  dplyr::arrange(dataset) %>% 
  dplyr::mutate(dataset_dim = paste0(dataset, 
                                     " (n=", n, ", p=", p, ")") %>% 
                  forcats::as_factor()) %>% 
  dplyr::filter(method %in% c("metaSPARSim", "SparseDOSSA 2")) %>% 
  ggplot(aes(x = method,
             y = time)) +
  geom_boxplot() +
  facet_wrap(~dataset_dim, scales = "free_y", nrow = 1) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  ylab("Time (seconds)") +
  ggtitle("Simulation time costs")


batchtools::loadRegistry(file.dir = "r_batchtools_reg/fit_SparseDOSSA2_full/", 
                         writeable = TRUE)
tb_time <- batchtools::getJobTable()
tb_time <- tb_time %>% 
  dplyr::select(job.id, time.running) %>% 
  dplyr::rename(i_job = job.id) %>% 
  dplyr::mutate(time.running = time.running * 30 %>% as.numeric())
load("/n/janson_lab/lab/sma/sparsedossa_paper/results/fitted/SparseDOSSA2_full/tb_job.RData")
tb_time <- tb_time %>% 
  dplyr::left_join(tb_job, by = "i_job") %>% 
  dplyr::group_by(i_job) %>% 
  dplyr::mutate(
    n = {
      load(paste0("data/physeqs/", dataset, ".RData"))
      x_samples <- 
        get(paste0("physeq_", dataset, "_bl")) %>%
        smar::otu_table2()
      ncol(x_samples)
    },
    p = {
      load(paste0("data/physeqs/", dataset, ".RData"))
      x_samples <- 
        get(paste0("physeq_", dataset, "_bl")) %>%
        smar::otu_table2()
      nrow(x_samples)
    },
    time_simulate = {
      load(paste0(
        "/n/janson_lab/lab/sma/sparsedossa_paper/results/fitted/SparseDOSSA2_full/", 
        "fit_", i_job, ".RData"))
      system.time(SparseDOSSA2::SparseDOSSA2(template = fit_EM, 
                                             n_sample = n,
                                             new_features = FALSE))[3]
    }
  ) %>% 
  dplyr::ungroup()

p_time2 <- tb_time %>%
  dplyr::mutate(time_running = as.numeric(time.running, units = "secs")) %>% 
  tidyr::pivot_longer(cols = c("time_running", "time_simulate"),
                      names_to = "category",
                      values_to = "time") %>% 
  dplyr::mutate(process = category %>% 
                  dplyr::recode("time_running" = "Fitting",
                                "time_simulate" = "Simulation")) %>% 
  dplyr::mutate(dataset = dataset %>% 
                  dplyr::recode("stool" = "Stool",
                                "vaginal" = "Vaginal",
                                "IBD" = "IBD") %>% 
                  factor(levels = c("Stool",
                                    "Vaginal",
                                    "IBD"))) %>% 
  dplyr::arrange(dataset) %>% 
  dplyr::mutate(dataset_dim = paste0(dataset, 
                                 " (n=", n, ", p=", p, ")") %>% 
                  forcats::as_factor()) %>% 
  ggplot(aes(x = process, y = log10(time))) +
  geom_boxplot() +
  facet_wrap(~dataset_dim, scales = "free") +
  theme_bw() +
  ylab("log10 time (seconds)") +
  ggtitle("SparseDOSSA 2 process time costs") +
  theme(axis.title.x = element_blank())

p <- cowplot::plot_grid(
  p_time,
  p_time2,
  labels = c("A", "B"),
  nrow = 2
)

ggsave(filename = "Supp_Materials/SuppFigures/SuppFig1.pdf", 
       p, width = 7.5, height = 6)

# figure check parameter distribution
set.seed(1)
l_p <- c("stool", "vaginal", "IBD") %>% 
  purrr::map2(c("Stool", "Vaginal", "IBD"),
                  function(i_dataset, i_dataset_name) {
    load(paste0("results/fitted/SparseDOSSA2_", i_dataset, ".RData"))
    param_original <- cbind("pi0" = fit_SparseDOSSA2$EM_fit$fit$pi0,
                            "mu" = fit_SparseDOSSA2$EM_fit$fit$mu,
                            "sigma" = fit_SparseDOSSA2$EM_fit$fit$sigma)
    param_new <- 
      SparseDOSSA2:::generate_featureParam(
        fit_SparseDOSSA2$F_fit,
        param_original = param_original, 
        n_feature = 1000)
    
    p <- rbind(data.frame(param_original, param = "Real-world"),
          data.frame(param_new, param = "SparseDOSSA 2")) %>% 
      dplyr::mutate(`logit(pi)` = log(pi0) - log(1 - pi0)) %>% 
      dplyr::select(-pi0) %>% 
      tidyr::pivot_longer(c("logit(pi)", "mu", "sigma"),
                          names_to = "Parameter",
                          values_to = "Value") %>% 
      dplyr::mutate(Parameter = factor(Parameter, levels = c("logit(pi)", "mu", "sigma"))) %>% 
      ggplot(aes(x = Value,
                 color = param)) +
      geom_density(alpha = 0.5,
                   fill = "grey") +
      scale_color_manual(values = c("Real-world"="black",
                                    "SparseDOSSA 2" = "red")) +
      # scale_fill_manual(values = c("Real-world"="black",
      #                               "SparseDOSSA 2" = "red")) +
      facet_wrap( ~ Parameter, scales = "free", nrow = 1) +
      ggtitle(i_dataset_name) +
      theme_bw() +
      theme(axis.title.x = element_blank())
    
    if(i_dataset == "stool")
      p <- 
      p + 
      theme(legend.position = c(0.01, 0.99),
                     legend.justification = c(0, 1),
                     legend.title = element_blank(),
                     legend.background = element_blank()) +
      ylab("Density")
    else
      p <- p + theme(legend.position = "none",
                     axis.title.y = element_blank())
    
    return(p)
  }) 
l_p %>% 
  cowplot::plot_grid(plotlist = ., 
                     nrow = 3, align = "v") %>% 
  ggsave(filename = "Supp_Materials/SuppFigures/SuppFig2.pdf",
         ., width = 8, height = 9)
  
df_plot %>% 
  ggplot(aes(x = Value,
             color = param)) +
  geom_density() +
  facet_grid(dataset ~ Parameter, scales = "free")


# figure 2
load("results/Evaluate/tb_job.RData")

levels_dataset <- c("stool" = "Stool",
                    "vaginal" = "Vaginal",
                    "IBD" = "IBD")
levels_method <- c("DMM" = "DM",
                   "MVN" = "Log MVN",
                   "metaSPARSim" = "metaSPARSim",
                   "SparseDOSSA2" = "SparseDOSSA 2")

# panel A
set.seed(1)
pA <- levels_dataset %>%
  purrr::imap_dfr(
    function(i_level_dataset, 
             i_dataset) {
      
      i_tb_A <- tb_job %>%
        dplyr::filter(method == "SparseDOSSA2",
                      dataset == i_dataset)
      
      # original data
      load(paste0("data/physeqs/", i_dataset, ".RData"))
      original_data <- get(paste0("physeq_", i_dataset, "_bl")) %>%
        smar::otu_table2()
      # original_data <- original_data[, setdiff(seq_len(ncol(original_data)),
      #                                          i_tb_A$training[[1]])]
      
      # simulated data
      sim_data_fill <- i_tb_A$i_job_sim %>% 
        purrr::map(function(i_job_sim) {
          load(paste0("results/Simulation/", i_job_sim, ".RData"))
          
          sim_data_fill <- matrix(0, nrow = nrow(original_data), ncol = ncol(sim_data))
          rownames(sim_data_fill) <- rownames(original_data)
          sim_data_fill[rownames(sim_data), ] <- sim_data
          
          return(sim_data_fill)
        }) %>% 
        purrr::reduce(cbind)
      
      sim_data_fill <- sim_data_fill[, sample.int(ncol(sim_data_fill), ncol(original_data))]

      dataset_combined <- cbind(original_data, sim_data_fill) %>%
        apply(2, function(x) x / sum(x))
      D <- vegan::vegdist(t(dataset_combined), method = "bray")
      fit_pcoa <- ape::pcoa(D)
      tb_plot <-
        tibble::tibble(sample = c(rep("Real-world", ncol(original_data)),
                                      rep("SparseDOSSA 2", ncol(sim_data_fill))) %>%
                         forcats::as_factor(),
                       dataset = i_level_dataset) %>%
        dplyr::mutate(PCo1 = fit_pcoa$vectors[, 1],
                      PCo2 = fit_pcoa$vectors[, 2])
      
      return(tb_plot)
    }) %>% 
  dplyr::mutate(dataset = factor(dataset, levels = levels_dataset)) %>% 
  ggplot(aes(x = PCo1, y = PCo2, color = sample)) +
  geom_point() +
  scale_color_manual(values = c("Real-world" = "black",
                                "SparseDOSSA 2" = "red")) +
  theme_bw() +
  theme(legend.position = c(0.01, 0.7),
        legend.justification = c(0, 0),
        legend.box.background = element_rect(colour = "black"),
        legend.title=element_blank(),
        legend.margin=margin(c(0,0,0,0)),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  facet_wrap(.~dataset, ncol = 1, scales = "free",
             strip.position = "left") +
  ggtitle("Bray-Curtis MDS")

# panel B
load("results/Evaluate/R2.RData")
dir_output_MVN <- "/n/janson_lab/lab/sma/sparsedossa_paper/results/fitted/MVN"
load(paste0(dir_output_MVN, "/tb_results_MVN.RData"))
tb_results_R2 <- 
  tb_job %>% 
  dplyr::mutate(R2 = R2) %>% 
  rbind(tb_results_MVN) %>% 
  dplyr::mutate(
    dataset = dataset %>% 
      dplyr::recode(
        "vaginal" = "Vaginal",
        "stool" = "Stool"
      ) %>% 
      factor(levels = levels_dataset),
    method = method %>% 
      dplyr::recode(
        "MVN" = "Log MVN",
        "DMM" = "DM",
        "SparseDOSSA2" = "SparseDOSSA 2"
      ) %>% 
      factor(levels = c("Original", levels_method)))

pB <- tb_results_R2 %>% 
  dplyr::filter(method != "Log MVN") %>% 
  ggplot(aes(x = method, y = R2)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(dataset~., scales = "free_y") +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 15, vjust = 1, hjust=1),
        axis.title.x = element_blank()) +
  ggtitle(expression(paste("PERMANOVA ", R^{2}, "\n vs. real-world"))) +
  ylab(expression(R^{2}))

pB_response <- tb_results_R2 %>% 
  ggplot(aes(x = method, y = R2)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~dataset, nrow = 1, scales = "free_y") +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 15, vjust = 1, hjust=1),
        axis.title.x = element_blank()) +
  ggtitle(expression(paste("PERMANOVA ", R^{2}, "\n vs. real-world"))) +
  ylab(expression(R^{2}))

ggsave(pB_response, filename = "Supp_Materials/Response_Letter/fig_MVN.jpeg",
       width = 10, height = 4)


tb_test <- tidyr::expand_grid(Dataset = c("Stool", "Vaginal", "IBD"),
                              `Method` = "SparseDOSSA 2",
                              `Existing method` = c("DM", "Log MVN", "metaSPARSim")) %>% 
  dplyr::mutate(`P value` = Dataset %>%
                  purrr::map2_dbl(
                    `Existing method`,
                    ~ {
                      R21 <- tb_results_R2 %>% 
                        dplyr::filter(dataset == .x,
                                      method == "SparseDOSSA 2") %>% 
                        extract2("R2")
                      
                      R22 <- tb_results_R2 %>% 
                        dplyr::filter(dataset == .x,
                                      method == .y) %>% 
                        extract2("R2")
                      
                      wilcox.test(R21, R22, paired = FALSE)$p.value
                    })) %>% 
  dplyr::filter(`Existing method` != "Log MVN")
readr::write_csv(tb_test, path = "Supp_Materials/SuppTables/SuppTable2.csv")

# panel C
features <- c("Prevotella copri" =
                "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella|s__Prevotella_copri",
              "Lactobacillus crispatus" =
                "k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Lactobacillus|s__Lactobacillus_crispatus",
              "Escherichia coli" = "k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia_coli"
)

tb_plot <- names(levels_dataset) %>%
  purrr::map2_dfr(
    names(features), 
    function(i_dataset, 
             i_feature) {
      
      # original data
      load(paste0("data/physeqs/", i_dataset, ".RData"))
      original_data <- get(paste0("physeq_", i_dataset, "_bl")) %>%
        smar::otu_table2()
      # original_data <- original_data[, setdiff(seq_len(ncol(original_data)),
      #                                          i_tb_A$training[[1]])]
      
      # simulated data
      i_tb <- tb_job %>%
        dplyr::filter(method != "Original") %>% 
        rbind(tb_results_MVN %>% dplyr::select(-R2)) %>% 
        dplyr::filter(dataset == i_dataset)
      
      
      sim_data_fill <- names(levels_method) %>% 
        purrr::map(function(i_method) {
          i_tb_C <- i_tb %>% 
            dplyr::filter(method == i_method)
          
          sim_data <- i_tb_C$i_job_sim %>% 
            purrr::map(function(i_job_sim) {
              if(i_tb_C$method[1] != "MVN") {
                load(paste0("results/Simulation/", i_job_sim, ".RData"))
              } else {
                load(paste0(dir_output_MVN, "/sim_", i_job_sim, ".RData"))
              }
              
              sim_data_fill <- matrix(0, nrow = nrow(original_data), ncol = ncol(sim_data))
              rownames(sim_data_fill) <- rownames(original_data)
              sim_data_fill[rownames(sim_data), ] <- sim_data
              
              return(sim_data_fill)
            }) %>% 
            purrr::reduce(cbind)
          return(sim_data)
        }) %>% 
        purrr::reduce(cbind)
      
      dataset_combined <- cbind(original_data, sim_data_fill) %>%
        apply(2, function(x) x / sum(x))
      
      tb_plot <- tibble::tibble(
        sample = c(rep("Real-world", ncol(original_data)),
                   rep(c("DM", "Log MVN", "metaSPARSim", "SparseDOSSA 2"), 
                       each = ncol(sim_data_fill) / 4)),
        Abundance = log10(dataset_combined[features[i_feature], ] + 1e-6),
        dataset = levels_dataset[i_dataset],
        feature = i_feature
      )
      
      return(tb_plot)
    }) %>% 
  dplyr::mutate(dataset = factor(dataset, levels = levels_dataset),
                feature = factor(feature, levels = names(features)),
                sample = factor(sample, levels = c("Real-world", levels_method)))

pC <- tb_plot %>% 
  dplyr::filter(sample != "Log MVN") %>% 
  ggplot(aes(x = Abundance, color = sample)) +
  stat_ecdf() +
  scale_color_manual(values = c("Real-world" = "black",
                                "DM" = "orange",
                                "Log MVN" = "darkgreen",
                                "metaSPARSim" = "blue",
                                "SparseDOSSA 2" = "red")) +
  theme_bw() +
  theme(legend.position = c(0.99, 0.69),
        legend.justification = c(1, 0),
        legend.box.background = element_rect(colour = "black"),
        legend.title=element_blank(),
        legend.margin=margin(c(0,0,0,0)),
        strip.background = element_blank()) +
  facet_grid(feature~.) +
  ggtitle("Per-feature distribution") +
  xlab("Log relative abundance") +
  ylab("Empirical CDF")

# panel D
tb_results_ks <- 
  tb_job %>% 
  dplyr::filter(method != "Original") %>% 
  rbind(tb_results_MVN %>% dplyr::select(-R2)) %>% 
  dplyr::group_by(dataset, method) %>% 
  dplyr::summarise(tb_results_ks = {
    load(paste0("data/physeqs/", dataset[1], ".RData"))
    
    original_data <- get(paste0("physeq_", dataset[1], "_bl")) %>% 
      smar::otu_table2()  %>% 
      apply(2, function(x) x / sum(x))
    
    sim_data <- i_job_sim %>% 
      purrr::map( ~ {
        if(method[1] != "MVN") {
          load(paste0("results/Simulation/", .x, ".RData")) 
        } else {
          load(paste0(dir_output_MVN, "/sim_", .x, ".RData")) 
        }
        sim_data_fill <- matrix(0, nrow = nrow(original_data), ncol = ncol(sim_data))
        rownames(sim_data_fill) <- rownames(original_data)
        sim_data_fill[rownames(sim_data), ] <- sim_data
        return(sim_data_fill)
      }) %>% 
      purrr::reduce(cbind) %>% 
      apply(2, function(x) x / sum(x))
    ks <- rownames(original_data) %>% 
      purrr::map_dbl(~ {
        ks.test(original_data[.x, ],
                sim_data[.x, ])$statistic
      })
    tb_results_ks <- tibble::tibble(
      dataset = dataset[1],
      method = method[1],
      feature = rownames(original_data),
      ks = ks
    )
    list(tb_results_ks)
  })
tb_results_ks <- tb_results_ks$tb_results_ks %>% 
  purrr::reduce(rbind)

pD <- tb_results_ks %>% 
  dplyr::mutate(
    dataset = dataset %>% 
      dplyr::recode("stool" = "Stool",
                    "vaginal" = "Vaginal") %>% 
      factor(levels = levels_dataset),
    method = method %>% 
      dplyr::recode(
        "MVN" = "Log MVN",
        "DMM" = "DM",
        "SparseDOSSA2" = "SparseDOSSA 2") %>% 
      factor(levels = levels_method)
  ) %>% 
  dplyr::filter(method != "Log MVN") %>% 
  tidyr::pivot_wider(id_cols = c(dataset, feature),
                     names_from = method,
                     values_from = ks) %>% 
  tidyr::pivot_longer(cols = -c(dataset, feature, `SparseDOSSA 2`),
                      names_to = "method",
                      values_to = "Other methods") %>% 
  dplyr::mutate(method = factor(method, levels = levels_method)) %>% 
  ggplot(aes(x = `SparseDOSSA 2`, y = `Other methods`, color = method)) +
  geom_point() +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  facet_grid(dataset~ .) +
  scale_color_manual(values = c("DM" = "orange", "Log MVN" = "darkgreen",
                                "metaSPARSim" = "blue")) +
  theme(legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.background = element_rect(colour = "black"),
        legend.title=element_blank(),
        legend.margin=margin(c(0,0,0,0)),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  ggtitle("Per-feature K-S statistic")

ggsave(pD, filename = "figures/Fig2/tmp.jpeg",
       width = 4,
       height = 12)

pD_response <- tb_results_ks %>% 
  dplyr::mutate(
    dataset = dataset %>% 
      dplyr::recode("stool" = "Stool",
                    "vaginal" = "Vaginal") %>% 
      factor(levels = levels_dataset),
    method = method %>% 
      dplyr::recode(
        "MVN" = "Log MVN",
        "DMM" = "DM",
        "SparseDOSSA2" = "SparseDOSSA 2") %>% 
      factor(levels = levels_method)
  ) %>% 
  tidyr::pivot_wider(id_cols = c(dataset, feature),
                     names_from = method,
                     values_from = ks) %>% 
  tidyr::pivot_longer(cols = -c(dataset, feature, `SparseDOSSA 2`),
                      names_to = "method",
                      values_to = "Other methods") %>% 
  dplyr::mutate(method = factor(method, levels = levels_method)) %>% 
  ggplot(aes(x = `SparseDOSSA 2`, y = `Other methods`, color = method)) +
  geom_point() +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  facet_grid(dataset~ .) +
  scale_color_manual(values = c("DM" = "orange", "Log MVN" = "darkgreen",
                                "metaSPARSim" = "blue")) +
  theme(legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.box.background = element_rect(colour = "black"),
        legend.title=element_blank(),
        legend.margin=margin(c(0,0,0,0)),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  ggtitle("Per-feature K-S statistic")

tb_results_prev <- 
  tb_job %>% 
  dplyr::filter(method != "Original") %>% 
  rbind(tb_results_MVN %>% dplyr::select(-R2)) %>% 
  dplyr::group_by(dataset, method) %>% 
  dplyr::summarise(tb_results_ks = {
    load(paste0("data/physeqs/", dataset[1], ".RData"))
    
    original_data <- get(paste0("physeq_", dataset[1], "_bl")) %>% 
      smar::otu_table2()  %>% 
      apply(2, function(x) x / sum(x))
    
    sim_data <- i_job_sim %>% 
      purrr::map( ~ {
        if(method[1] != "MVN") {
          load(paste0("results/Simulation/", .x, ".RData")) 
        } else {
          load(paste0(dir_output_MVN, "/sim_", .x, ".RData")) 
        }
        sim_data_fill <- matrix(0, nrow = nrow(original_data), ncol = ncol(sim_data))
        rownames(sim_data_fill) <- rownames(original_data)
        sim_data_fill[rownames(sim_data), ] <- sim_data
        return(sim_data_fill)
      }) %>% 
      purrr::reduce(cbind) %>% 
      apply(2, function(x) x / sum(x))
    prev_original <- apply(original_data > 0, 1, mean)
    prev_sim <- apply(sim_data > 0, 1, mean)
    tb_results_ks <- tibble::tibble(
      dataset = dataset[1],
      method = method[1],
      feature = rownames(original_data),
      prev_original = prev_original,
      prev_sim = prev_sim
    )
    list(tb_results_ks)
  })

tb_results_prev <- tb_results_prev$tb_results_ks %>% 
  purrr::reduce(rbind)

p_prev <- tb_results_prev %>% 
  ggplot(aes(x = prev_original,
             y = prev_sim,
             color = method)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  facet_grid(dataset~., scales = "free_y")

p_ks <- cowplot::plot_grid(
  pD_response,
  p_prev,
  nrow = 1
)

ggsave(p_ks, filename = "Supp_Materials/Response_Letter/fig_KS.jpeg",
       width = 8, height = 8)

pD_box <- tb_results_ks %>% 
  dplyr::filter(method != "MVN") %>% 
  dplyr::mutate(
    dataset = dataset %>% 
      dplyr::recode("stool" = "Stool",
                    "vaginal" = "Vaginal") %>% 
      factor(levels = levels_dataset),
    method = method %>% 
      dplyr::recode("DMM" = "DM",
                    "SparseDOSSA2" = "SparseDOSSA 2") %>% 
      factor(levels = levels_method)
  ) %>% 
  ggplot(aes(x = method, y = ks)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~dataset, scales = "free") +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) +
  ggtitle("Per-feature K-S statistic")
ggsave(pD_box, filename = "Supp_Materials/Response_Letter/fig_KS_boxplot.jpeg",
       width = 10, height = 4)

tb_test <- tidyr::expand_grid(Dataset = c("stool", "vaginal", "IBD"),
                              `Method` = "SparseDOSSA 2",
                              `Existing method` = c("DMM", "metaSPARSim")) %>% 
  dplyr::mutate(`P value` = Dataset %>%
                  purrr::map2_dbl(
                    `Existing method`,
                    ~ {
                      KS1 <- tb_results_ks %>% 
                        dplyr::filter(dataset == .x,
                                      method == "SparseDOSSA2") %>% 
                        extract2("ks")
                      
                      KS2 <- tb_results_ks %>% 
                        dplyr::filter(dataset == .x,
                                      method == .y) %>% 
                        extract2("ks")
                      
                      wilcox.test(KS1, KS2, paired = TRUE)$p.value
                    })) %>% 
  dplyr::mutate(Dataset = Dataset %>% 
                  dplyr::recode("stool" = "Stool",
                                "vaginal" = "Vaginal")) %>% 
  dplyr::mutate(`Existing method` = `Existing method` %>% 
                  dplyr::recode("DMM" = "DM"))
readr::write_csv(tb_test, path = "Supp_Materials/SuppTables/SuppTable3.csv")


cowplot::plot_grid(pA, pB, pC, pD, labels = c("A", "B", "C", "D"),
                   align = "h", axis = "tb", nrow = 1,
                   rel_widths = c(0.9, 1.2, 1, 1)) %>% 
  ggsave("figures/Fig2/Fig2.pdf",
       ., width = 11, height = 7)
 