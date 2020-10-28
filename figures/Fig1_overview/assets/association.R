mat_null <- rbinom(n = 4, size = 50, prob = 0.5) / 50
mat_null <- rbind(mat_null, 1 - mat_null)
dimnames(mat_null) <- list(paste0("Feature", 1:2), paste0("Sample", 1:4))
covar <- c(-1, 0, 0, 1)
covar_norm <- covar/sd(covar)

mat_1 <- mat_null
bug <- mat_1[1, ]
bug_mean <- mean(bug)
bug_sd <- sd(bug)
bug_stand <- (bug - bug_mean) / bug_sd
covar_norm <- covar_norm * sd(covar)
bug_mod <- (bug + (covar_norm + bug_mean)) / 2
mat_1[1, ] <- bug_mod
mat_1_renorm <- apply(mat_1, 2, MMUPHin:::tss)

 (mat_null %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Feature") %>% 
  tidyr::gather(key = "Sample", value = "Abundance", -Feature) %>% 
  ggplot(aes(x = Sample, y = Abundance, fill = Feature)) +
  geom_bar(stat = "identity")) %>% 
  ggsave("figures/fig_model/rel1.pdf", ., width = 5, height = 4)

(mat_1 %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("Feature") %>% 
    tidyr::gather(key = "Sample", value = "Abundance", -Feature) %>% 
    ggplot(aes(x = Sample, y = Abundance, fill = Feature)) +
    geom_bar(stat = "identity")) %>% 
  ggsave("figures/fig_model/rel2.pdf", ., width = 5, height = 4)

(mat_1_renorm %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("Feature") %>% 
    tidyr::gather(key = "Sample", value = "Abundance", -Feature) %>% 
    ggplot(aes(x = Sample, y = Abundance, fill = Feature)) +
    geom_bar(stat = "identity")) %>% 
  ggsave("figures/fig_model/rel3.pdf", ., width = 5, height = 4)
