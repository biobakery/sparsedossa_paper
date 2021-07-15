library(DESeq2)

###########################
# Fit DESeq2 To A Dataset #
###########################

fit.DESeq2<-function(features, metadata, transformation = "NONE", MultTestCorrection = "BH"){
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default DESeq2 model. Use NONE.')
  
  
  # Random Effect Adjustment
  # if(!length(ID)==length(unique(ID))){
  #   stop('edgeR random effect model is currently not implemented.')
  # }
  
  formula <- as.formula(paste('~', paste(colnames(metadata), collapse = "+"), sep=''))
  
  
  # Fit Model
  x <- DESeqDataSetFromMatrix(countData = t(features), colData = metadata, design = formula)
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(x), 1, gm_mean)
  x = estimateSizeFactors(x, geoMeans = geoMeans)
  fit <- DESeq(x)
  if(dim(metadata)[2]>1){
    coef.vector<-rename.features(coef(fit)[,-1], 'coef')
    pvalMatrix<-get_pval_DESeq2(fit)
    pvalue.vector<-rename.features(pvalMatrix[,-1], 'pval')
    DD<-cbind.data.frame(coef.vector, pvalue.vector)
    DD<-DD[, !duplicated(colnames(DD))]
    DD$qval<-as.numeric(p.adjust(DD$pval, method = MultTestCorrection, n = nrow(DD)))
    DD<-DD[order(DD$qval, decreasing=FALSE),]
    DD<-dplyr::select(DD, c('feature', 'metadata'), everything())
    rownames(DD)<-NULL;
  }
  else{
    coef<-coef(fit)[,-1]
    pval<-results(fit,name=resultsNames(fit)[2])$pvalue
    DD<-cbind.data.frame(coef,pval)
    DD$qval<-as.numeric(p.adjust(DD$pval, method = MultTestCorrection, n = nrow(DD)))
    DD<-DD[order(DD$qval, decreasing=FALSE),]
    DD$feature<-rownames(DD)
    DD$metadata<- names(metadata)
    DD<-dplyr::select(DD, c('feature', 'metadata'), everything())
    rownames(DD)<-NULL;
  }
  return(DD)
}

# Get P-values from DESeq2 Fit
get_pval_DESeq2<-function(fit){
  List <- list()
  for(i in 1:length(resultsNames(fit))){
    List[[i]] <- results(fit,name=resultsNames(fit)[i])$pvalue
  }
  Matrix = do.call(cbind, List)
  rownames(Matrix)<-names(fit)
  colnames(Matrix)<-resultsNames(fit)
  return(Matrix)
}

###########################
# ANCOM#
###########################

devtools::source_url("https://github.com/FrederickHuangLin/ANCOM/blob/master/scripts/ancom_v2.1.R?raw=TRUE")

###########################
# Fit ANCOM To A Dataset #
###########################

fit.ANCOM = function(features, metadata, transformation = "NONE", MultTestCorrection = "BH") {
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default ANCOM model. Use NONE.')
  
  #############################################################################################
  # ANCOM standard pipeline for DA (adopted from https://github.com/FrederickHuangLin/ANCOM) #
  #############################################################################################
  
  ###############################
  # Step 1: ANCOM preprocessing #
  ###############################
  
  features_ancom<-as.data.frame(t(features))
  metadata_ancom<-tibble::rownames_to_column(metadata, 'ID')
  preprocess.ancom = feature_table_pre_process(feature_table = features_ancom, 
                                               meta_data = metadata_ancom, 
                                               sample_var = "ID", 
                                               group_var = NULL,
                                               out_cut = 0.05,
                                               zero_cut = 0.90,
                                               lib_cut = 1000,
                                               neg_lb = FALSE)
  feature_table = preprocess.ancom$feature_table 
  meta_data = preprocess.ancom$meta_data 
  struc_zero = preprocess.ancom$structure_zeros 
  
  
  #################
  # Step 2: ANCOM #
  #################
  
  # if(!length(ID)==length(unique(ID))){
  #   res <- ANCOM(feature_table = feature_table, 
  #                meta_data = meta_data, 
  #                struc_zero = struc_zero,
  #                main = colnames(metadata),
  #                alpha = 0.05,
  #                adj_formula = NULL,
  #                rand_formula = "~ 1 | ID",
  #                p_adj_method = MultTestCorrection) 
  # } else{
    res <- ANCOM(feature_table = feature_table, 
                 meta_data = meta_data, 
                 struc_zero = struc_zero,
                 main = colnames(metadata),
                 alpha = 0.05,
                 adj_formula = NULL,
                 rand_formula = NULL,
                 p_adj_method = MultTestCorrection)
  # }
  
  #########################################################
  # Step 3: Extract results and enforce meaningful format #
  #########################################################
  
  df <- data.frame(coef = res$out)
  df$pval<-1 # Fake p-values
  df$feature = df$coef.taxa_id
  df$metadata = names(metadata)
  df$pval[df$coef.detected_0.7] <- 0 # Fake p-values
  df$qval<-df$pval # Fake q-values
  df<-df[order(df$qval, decreasing=FALSE),]
  df<-dplyr::select(df, c('feature', 'metadata', 'pval', 'qval'))
  rownames(df)<-NULL;
  return(df)
}

library(limma)
fit.limmaVOOM = function(features, metadata, transformation = "NONE", MultTestCorrection = "BH") {
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default limmaVOOM model. Use NONE.')
  
  # Fit limmaVOOM
  x<-t(as.matrix(features)+1) # Convert to matrix, round up to nearest integer, and transpose
  design <- model.matrix(~., data=metadata)
  y <- voom(x,design,plot=FALSE)
  
  # Fit limma 
  # Random Effect Adjustment
  
  fit <- limma::lmFit(y,design)
  
  # Empirical Bayes Adjustment
  fit <- limma::eBayes(fit)
  
  
  if(dim(metadata)[2]>1){
    coef.vector<-rename.features(fit$coefficients[,-1], 'coef')
    pvalue.vector<-rename.features(fit$p.value[,-1], 'pval')
    DD<-cbind.data.frame(coef.vector, pvalue.vector)
    DD<-DD[, !duplicated(colnames(DD))]
    DD$qval<-as.numeric(p.adjust(DD$pval, method = MultTestCorrection, n = nrow(DD)))
    DD<-DD[order(DD$qval, decreasing=FALSE),]
    DD<-dplyr::select(DD, c('feature', 'metadata'), everything())
    rownames(DD)<-NULL;
  } 
  else{
    coef<-fit$coefficients[,-1]
    pval<-fit$p.value[,-1]
    DD<-cbind.data.frame(coef,pval)
    DD$qval<-as.numeric(p.adjust(DD$pval, method = MultTestCorrection, n = nrow(DD)))
    DD<-DD[order(DD$qval, decreasing=FALSE),]
    DD$feature<-rownames(DD)
    DD$metadata<- names(metadata)
    DD<-dplyr::select(DD, c('feature', 'metadata'), everything())
    rownames(DD)<-NULL;
  }
  return(DD)
}


##########################
# Fit edgeR To A Dataset #
##########################
library(edgeR)

fit.edgeR = function(features, metadata, transformation = "NONE", MultTestCorrection = "BH") {
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default edgeR model. Use NONE.')
  
  d <- DGEList(counts=t(features))
  d <- edgeR::calcNormFactors(d, method='TMM')
  
  # Random Effect Adjustment
  design <- model.matrix(~., data=metadata)
  d <- estimateGLMCommonDisp(d,design)
  d <- estimateGLMTrendedDisp(d,design)
  d <- estimateGLMTagwiseDisp(d,design)
  fit <- glmFit(d,design)
  if(dim(metadata)[2]>1){
    coef.vector<-rename.features(fit$coefficients[,-1], 'coef')
    pvalMatrix<-get_pval_edgeR(fit)
    pvalue.vector<-rename.features(pvalMatrix[,-1], 'pval')
    DD<-cbind.data.frame(coef.vector, pvalue.vector)
    DD<-DD[, !duplicated(colnames(DD))]
    DD$qval<-as.numeric(p.adjust(DD$pval, method = MultTestCorrection, n = nrow(DD)))
    DD<-DD[order(DD$qval, decreasing=FALSE),]
    DD<-dplyr::select(DD, c('feature', 'metadata'), everything())
    rownames(DD)<-NULL;
  }
  else{
    fit<-glmLRT(fit, 2)
    coef<-fit$coefficients[,-1]
    pval<-fit$table$PValue
    DD<-cbind.data.frame(coef,pval)
    DD$qval<-as.numeric(p.adjust(DD$pval, method = MultTestCorrection, n = nrow(DD)))
    DD<-DD[order(DD$qval, decreasing=FALSE),]
    DD$feature<-rownames(DD)
    DD$metadata<- names(metadata)
    DD<-dplyr::select(DD, c('feature', 'metadata'), everything())
    rownames(DD)<-NULL;
  }
  return(DD)
}


###############################################################
# Summarization funciton
###############################################################
eval_res_list = function(resi, alpha=0.05) {
  
  # Replace Missing Q-values to the Highest Possible Value 1.0
  resi[is.na(resi[, "qval"]), "qval"] <- 1
  
  # Evaluate Detection Performance
  
  time= mean(resi[,"time"], na.rm=TRUE)
  wh.pred = (resi[, "qval"] < alpha)
  wh.pos = which(wh.pred)
  wh.neg = which(!wh.pred)
  wh.TP = grep("[[:print:]]+\\_TP$", resi[, "pairwiseAssociation"])
  FPs = sum(!wh.pos %in% wh.TP)
  TPs = sum(wh.pos %in% wh.TP)
  TNs = sum(!wh.neg %in% wh.TP)
  FNs = sum(wh.neg %in% wh.TP)
  
  
  # Sensitivity: True Positives Divided by All Positives (Sum of True
  # Positives and False Negatives)
  Sensitivity = TPs/(TPs + FNs)
  
  # Specificity: True Negatives Divided by All Negatives (Sum of True
  # Negatives and False Positives)
  Specificity = TNs/(TNs + FPs)
  
  # False Discovery Rate: False Positives Divided by All Detected Positives
  FDR = if ((TPs + FPs) == 0) 
    0 else FPs/(TPs + FPs)
  # If no true positives, return NA's for irrelevant measures
  
  # FScore
  FScore<-2*TPs/(2*TPs + FNs + FPs) 
  
  # Matthew's Correlation Coefficient
  numerator <- (TPs * TNs - FPs * FNs)
  denominator <- sqrt((TPs + FPs)*(TPs + FNs)*(TNs + FPs)*(TNs + FNs))
  if(denominator == 0) denominator <- 1
  MCC <- numerator/denominator
  
  # AUC, pAUC (FPR < 0.20)
  wh.truth = (1:nrow(resi) %in% wh.TP)
  if (all(wh.truth=='TRUE')) {
    AUC=1
    pAUC=1
    fAUC=1
  } else if (all(wh.truth=='FALSE')) {
    AUC=0
    pAUC=0
  } else {
    pred <- prediction(as.numeric(wh.pred), factor(wh.truth, levels=c("TRUE", "FALSE")))
    AUC = performance(pred, "auc")@y.values[[1]]
    pAUC = performance(pred, "auc", fpr.stop=0.20)@y.values[[1]][[1]]
  }
  
  # Departure from Uniformity Under the Null
  AucAocVals <- AucAocFun(resi[!wh.truth, "pval"], plotIt = FALSE, pch = 20, type = "l")
  totalArea = AucAocVals["conservArea"] + AucAocVals["liberalArea"]
  names(totalArea) = "totalArea"  
  
  # Return
  return(c(Sensitivity = Sensitivity,
           Specificity=Specificity,
           FDR=FDR,
           FScore=FScore,
           MCC=MCC,
           AUC=AUC,
           pAUC=pAUC,
           AucAocVals["conservArea"],
           AucAocVals["liberalArea"],
           totalArea,
           time=time))
}