DA <- function(data, metadata,
               method, features_TP, i) {
  if(method == "Maaslin2") {
    tmp <- capture.output(
      fit_DA <- Maaslin2::Maaslin2(
        input_data = data,
        input_metadata = metadata,
        normalization = "TSS", 
        fixed_effects = "Exposure",
        plot_heatmap = FALSE,
        plot_scatter = FALSE,
        output = paste0("results/DA/Maaslin2_tmp/", i)))
    n_TP <- sum(fit_DA$results$pval < 0.05 & 
                  fit_DA$results$feature %in% features_TP,
                na.rm = TRUE)
    n_FP <- sum(fit_DA$results$pval < 0.05 & 
                  !(fit_DA$results$feature %in% features_TP),
                na.rm = TRUE)
  }
  if(method == "DESeq2") {
    tmp <- capture.output(
      fit_DA <- DESeq2::DESeqDataSetFromMatrix(
      countData = data + 1,
      colData = metadata,
      design= ~ Exposure))
    tmp <- capture.output(
      fit_DA <- DESeq2::DESeq(fit_DA))
    results <- DESeq2::results(fit_DA)
    n_TP <- sum(results$pvalue < 0.05 & 
                  rownames(results) %in% features_TP,
                na.rm = TRUE)
    n_FP <- sum(results$pvalue < 0.05 & 
                  !(rownames(results) %in% features_TP),
                na.rm = TRUE)
  }
  if(method == "metagenomeSeq") {
    tmp <- capture.output(
      Mobj <- metagenomeSeq::newMRexperiment(
        counts = data,
        phenoData = Biobase::AnnotatedDataFrame(metadata),
        featureData = Biobase::AnnotatedDataFrame(data.frame(rownames(data),
                                                             row.names = rownames(data)))
      ))
    tmp <- capture.output(
      Mobj_norm <- metagenomeSeq::cumNorm(Mobj))
    mod <- model.matrix(~ 1 + Exposure, data = metadata)
    tmp <- capture.output(fit_DA <- metagenomeSeq::fitFeatureModel(Mobj_norm,
                                                                   mod))
    results <- metagenomeSeq::MRcoefs(fit_DA, number = nrow(data))
    n_TP <- sum(results$pvalues < 0.05 & 
                  rownames(results) %in% features_TP,
                na.rm = TRUE)
    n_FP <- sum(results$pvalues < 0.05 & 
                  !(rownames(results) %in% features_TP),
                na.rm = TRUE)
  }
  return(c(n_TP, n_FP)) 
}
