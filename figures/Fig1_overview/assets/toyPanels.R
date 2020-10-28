mat_abundance <- matrix(c(0.03, 0.2, 0.77,
                          0.02, 0.3, 0.68,
                          0.04, 0, 0.96), 
                        nrow = 3, ncol = 3)
mat_count <- matrix(c(0.04*15000, 0.17*15000, 0.79*15000,
                      0, 0.33*3000, 0.67*3000,
                      0.02*20000, 0, 0.98*20000),
                    nrow = 3, ncol = 3)
rownames(mat_abundance) <- rownames(mat_count) <- paste0("Microbe ", 1:3)
colnames(mat_abundance) <- colnames(mat_count) <- paste0("Environment ", 1:3)

p_abundance <- mat_abundance %>% 
  as.data.frame(check.names = FALSE) %>% 
  tibble::rownames_to_column("Feature") %>% 
  tidyr::gather(key = "Environment", value = "Abundance", -Feature) %>% 
  ggplot(aes(x = "", y = Abundance, fill = Feature)) +
  geom_bar(stat = "identity") + 
  coord_polar("y", start = 0) +
  facet_grid(.~ Environment) +
  smar::rotate_xaxis(30) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")
ggsave("figures/fig_sequencedata/assets/toy_proportions.pdf",
       p_abundance,
       width = 4, height = 1.5)

df_count <- colnames(mat_count) %>% 
  purrr::map_dfr(function(i_environment){
    i_count <- round(mat_count[, i_environment] / 10)
    i_df <- data.frame(y = 1:i_count[3], Feature = "Microbe 3", stringsAsFactors = FALSE)
    if(i_count[2] > 0)
      i_df <- rbind(i_df, data.frame(y = (i_count[3] + 1):(i_count[2] + i_count[3]), 
                                     Feature = "Microbe 2", 
                                     stringsAsFactors = FALSE))
    if(i_count[1] > 0)
      i_df <- rbind(i_df, data.frame(y = (i_count[2] + i_count[3] + 1):sum(i_count), 
                                     Feature = "Microbe 1",
                                     stringsAsFactors = FALSE))
    i_df <- i_df %>% 
      dplyr::mutate(Environment = i_environment)
    return(i_df)
  }) %>% 
  dplyr::mutate(Environment = factor(Environment, levels = colnames(mat_count)),
                Feature = factor(Feature, levels = rownames(mat_count)))

png_sequencer <- "figures/fig_sequencedata/assets/table-graphic-system-hiseq.png"
p_count_no <- df_count %>% 
  ggplot(aes(x = Environment, y = as.factor(y))) +
  geom_tile(fill = "white", width = 0.8) +
  geom_segment(data = data.frame(x = c(0.6, 1.4, 1.6, 2.4, 2.6, 3.4),
                                 ystart = rep(0, 6), 
                                 yend = rep(apply(mat_count, 2, sum) / 10, each = 2)),
               aes(x = x, xend = x, y = ystart, yend = yend), size = 1) +
  geom_segment(data = data.frame(xstart = rep(c(0.5875, 1.5875, 2.5875), 2),
                                 xend = rep(c(1.4125, 2.4125, 3.4125), 2),
                                 y = c(apply(mat_count, 2, sum) / 10, rep(0.5, 3))),
               aes(x = xstart, xend = xend, y = y, yend = y), size = 1) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")
p_count_no <- cowplot::ggdraw(p_count_no) +
  cowplot::draw_image(png_sequencer, x = 0, y = 1, hjust = 0, vjust = 1, width = 0.25, height = 0.25)
ggsave("figures/fig_sequencedata/assets/toy_counts_no.pdf",
       p_count_no, 
       width = 4, height = 4)

random_ys <- sample.int(sum(mat_count[, 1]) / 10, size = 500)
colors <- smar::gg_color_hue(rownames(mat_count))

p_count_some1 <- df_count %>% 
  dplyr::mutate(Feature_plot = dplyr::case_when(
    Environment == "Environment 1" & y %in% random_ys[1] ~ as.character(Feature),
    TRUE ~ NA_character_
  ) %>% 
    factor(levels = rownames(mat_count))) %>% 
  dplyr::arrange(!is.na(Feature_plot)) %>% 
  ggplot(aes(x = Environment, y = as.factor(y))) +
  geom_tile(aes(fill = Feature_plot), width = 0.8, height = 2) +
  geom_segment(data = data.frame(x = c(0.6, 1.4, 1.6, 2.4, 2.6, 3.4),
                                 ystart = rep(0, 6), 
                                 yend = rep(apply(mat_count, 2, sum) / 10, each = 2)),
               aes(x = x, xend = x, y = ystart, yend = yend), size = 1) +
  geom_segment(data = data.frame(xstart = rep(c(0.5875, 1.5875, 2.5875), 2),
                                 xend = rep(c(1.4125, 2.4125, 3.4125), 2),
                                 y = c(apply(mat_count, 2, sum) / 10, rep(0.5, 3))),
               aes(x = xstart, xend = xend, y = y, yend = y), size = 1) +
  scale_fill_manual(values = colors, na.value = "white") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")
ggsave("figures/fig_sequencedata/assets/toy_counts_some1.pdf",
       p_count_some1, width = 4, height = 4)

p_count_some2 <- df_count %>% 
  dplyr::mutate(Feature_plot = dplyr::case_when(
    Environment == "Environment 1" & y %in% random_ys[1:10] ~ as.character(Feature),
    TRUE ~ NA_character_
  ) %>% 
    factor(levels = rownames(mat_count))) %>% 
  dplyr::arrange(!is.na(Feature_plot)) %>% 
  ggplot(aes(x = Environment, y = as.factor(y))) +
  geom_tile(aes(fill = Feature_plot), width = 0.8, height = 2) +
  geom_segment(data = data.frame(x = c(0.6, 1.4, 1.6, 2.4, 2.6, 3.4),
                                 ystart = rep(0, 6), 
                                 yend = rep(apply(mat_count, 2, sum) / 10, each = 2)),
               aes(x = x, xend = x, y = ystart, yend = yend), size = 1) +
  geom_segment(data = data.frame(xstart = rep(c(0.5875, 1.5875, 2.5875), 2),
                                 xend = rep(c(1.4125, 2.4125, 3.4125), 2),
                                 y = c(apply(mat_count, 2, sum) / 10, rep(0.5, 3))),
               aes(x = xstart, xend = xend, y = y, yend = y), size = 1) +
  scale_fill_manual(values = colors, na.value = "white") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")
ggsave("figures/fig_sequencedata/assets/toy_counts_some2.pdf",
       p_count_some2, width = 4, height = 4)

p_count_some3 <- df_count %>% 
  dplyr::mutate(Feature_plot = dplyr::case_when(
    Environment == "Environment 1" & y %in% random_ys[1:150] ~ as.character(Feature),
    TRUE ~ NA_character_
  ) %>% 
    factor(levels = rownames(mat_count))) %>% 
  dplyr::arrange(!is.na(Feature_plot)) %>% 
  ggplot(aes(x = Environment, y = as.factor(y))) +
  geom_tile(aes(fill = Feature_plot), width = 0.8, height = 2) +
  geom_segment(data = data.frame(x = c(0.6, 1.4, 1.6, 2.4, 2.6, 3.4),
                                 ystart = rep(0, 6), 
                                 yend = rep(apply(mat_count, 2, sum) / 10, each = 2)),
               aes(x = x, xend = x, y = ystart, yend = yend), size = 1) +
  geom_segment(data = data.frame(xstart = rep(c(0.5875, 1.5875, 2.5875), 2),
                                 xend = rep(c(1.4125, 2.4125, 3.4125), 2),
                                 y = c(apply(mat_count, 2, sum) / 10, rep(0.5, 3))),
               aes(x = xstart, xend = xend, y = y, yend = y), size = 1) +
  scale_fill_manual(values = colors, na.value = "white") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")
ggsave("figures/fig_sequencedata/assets/toy_counts_some3.pdf",
       p_count_some3, width = 4, height = 4)

p_count_some4 <- df_count %>% 
  dplyr::mutate(Feature_plot = dplyr::case_when(
    Environment == "Environment 1" ~ as.character(Feature),
    TRUE ~ NA_character_
  ) %>% 
    factor(levels = rownames(mat_count))) %>% 
  ggplot(aes(x = Environment, y = as.factor(y))) +
  geom_tile(aes(fill = Feature_plot), width = 0.8, height = 3) +
  geom_segment(data = data.frame(x = c(0.6, 1.4, 1.6, 2.4, 2.6, 3.4),
                                 ystart = rep(0, 6), 
                                 yend = rep(apply(mat_count, 2, sum) / 10, each = 2)),
               aes(x = x, xend = x, y = ystart, yend = yend), size = 1) +
  geom_segment(data = data.frame(xstart = rep(c(0.5875, 1.5875, 2.5875), 2),
                                 xend = rep(c(1.4125, 2.4125, 3.4125), 2),
                                 y = c(apply(mat_count, 2, sum) / 10, rep(0.5, 3))),
               aes(x = xstart, xend = xend, y = y, yend = y), size = 1) +
  scale_fill_manual(values = colors, na.value = "white") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")
ggsave("figures/fig_sequencedata/assets/toy_counts_some4.pdf",
       p_count_some4, width = 4, height = 4)

p_count_all <- df_count %>% 
  ggplot(aes(x = Environment, y = as.factor(y))) +
  geom_tile(aes(fill = Feature), width = 0.8, height = 3) +
  geom_segment(data = data.frame(x = c(0.6, 1.4, 1.6, 2.4, 2.6, 3.4),
                                 ystart = rep(0, 6), 
                                 yend = rep(apply(mat_count, 2, sum) / 10, each = 2)),
               aes(x = x, xend = x, y = ystart, yend = yend), size = 1) +
  geom_segment(data = data.frame(xstart = rep(c(0.5875, 1.5875, 2.5875), 2),
                                 xend = rep(c(1.4125, 2.4125, 3.4125), 2),
                                 y = c(apply(mat_count, 2, sum) / 10, rep(0.5, 3))),
               aes(x = xstart, xend = xend, y = y, yend = y), size = 1) +
  scale_fill_manual(values = colors, na.value = "white") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")
ggsave("figures/fig_sequencedata/assets/toy_counts_all.pdf",
       p_count_all, width = 4, height = 4)

# p_count <- mat_count %>% 
#   as.data.frame(check.names = FALSE) %>% 
#   tibble::rownames_to_column("Feature") %>% 
#   tidyr::gather(key = "Environment", value = "Count", -Feature) %>% 
#   ggplot(aes(x = Environment, y = Count, fill = Feature)) +
#   geom_bar(stat = "identity") +
#   theme(legend.position = "none",
#         axis.title.x = element_blank()) +
#   smar::rotate_xaxis(30)

# p_ra <- mat_count %>% 
#   apply(2, function(x) x / sum(x)) %>% 
#   as.data.frame(check.names = FALSE) %>% 
#   tibble::rownames_to_column("Feature") %>% 
#   tidyr::gather(key = "Environment", value = "Abundance", -Feature) %>% 
#   ggplot(aes(x = Environment, y = Abundance, fill = Feature)) +
#   geom_bar(stat = "identity") +
#   theme(legend.position = "none",
#         axis.title.x = element_blank()) +
#   smar::rotate_xaxis(30)
# 
# cowplot::plot_grid(p_abundance, p_count, p_ra,
#                    ncol = 1, align = "hv") %>% 
#   ggsave("figures/fig_sequencedata/assets/toypanels.pdf", ., width = 4, height = 12)

p_countNumbers <- mat_count %>% 
  as.data.frame(check.names = FALSE) %>% 
  tibble::rownames_to_column("Feature") %>% 
  tidyr::gather(key = "Environment", value = "Count", -Feature) %>% 
  dplyr::mutate(Feature = factor(Feature, levels = c("Microbe 3", "Microbe 2", "Microbe 1"))) %>% 
  ggplot(aes(x = Environment, y = Feature)) +
  geom_tile(fill = "white", color = "black", size = 2) +
  geom_text(aes(label = Count), size = 8) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank())
ggsave("figures/fig_sequencedata/assets/table.pdf", p_countNumbers,
       width = 4, height = 3)
