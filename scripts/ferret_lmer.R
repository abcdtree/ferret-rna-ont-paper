library(lme4)
library(readr)
library(rjson)
library(dplyr)
library(ggplot2)


gene_lm_test <- function(df, gene, plot=FALSE){
  df_sub <- df[df$gene_id==gene,]
  
  #in vitro mock vs treated
  df_sub_vitro <- df_sub[df_sub$sample=="vitro",]
  m_count <- df_sub_vitro %>% count(barcode)
  m_count_2 <- df_sub_vitro %>% count(stype)
  if (nrow(m_count[m_count$n <=10, ]) != 0 | nrow(m_count_2) < 2) {
    coefficients_vitro <- 0
    pv_vitro <- 1
  }
  else{
    lm_vitro <- tryCatch(
      {lm(polyA_length ~ stype, data = df_sub_vitro)},
      error=function(e){
        print(gene)
        NA
      }
      )
    if (anyNA(lm_vitro)){
      print(gene)
    }
    coefficients_vitro <- summary(lm_vitro)$coefficients[,1][[2]]
    pv_vitro <- summary(lm_vitro)$coefficients[,4][[2]]
    df_sub_vitro$mtype <- ifelse(df_sub_vitro$stype=="Mock", 0, 1)
    if(plot == TRUE){
      #print("here")
      k <- ggplot(df_sub_vitro, aes(x = mtype, y = polyA_length)) +
        geom_point() +
        stat_smooth(method = "lm")
      print(k)
    }
  }
  #lm_vitro <- lm(polyA_length ~ stype, data = df_sub_vitro)
  #coefficients_vitro <- summary(lm_vitro)$coefficients[,1][[2]]
  #in vivo mock vs 24h
  df_sub_vivo <- df_sub[df_sub$sample=="vivo",]
  df_sub_vivo <- df_sub_vivo[df_sub_vivo$stype=="Mock"| df_sub_vivo$stype=="24H",]
  m_count <- df_sub_vivo %>% count(barcode)
  m_count_2 <- df_sub_vivo %>% count(stype)
  if (nrow(m_count[m_count$n <=10, ]) != 0 | nrow(m_count_2) < 2) {
    coefficients_vivo <- 0
    pv_vivo <- 1
  }
  else{
    #lm_vivo <- lm(polyA_length ~ stype, data = df_sub_vivo)
    lm_vivo <- tryCatch(
      {lm(polyA_length ~ stype, data = df_sub_vivo)},
      error=function(e){
        print(gene)
        NA
      }
    )
    if (anyNA(lm_vivo)){
      print(gene)
    }
    coefficients_vivo <- summary(lm_vivo)$coefficients[,1][[2]]
    pv_vivo <- summary(lm_vivo)$coefficients[,4][[2]]
    df_sub_vivo$mtype <- ifelse(df_sub_vivo$stype=="Mock", 0, 1)
    if(plot == TRUE){
      k <- ggplot(df_sub_vivo, aes(x = mtype, y = polyA_length)) +
        geom_point() +
        stat_smooth(method = "lm")
      print(k)
    }
  }
  #lm_vivo <- lm(polyA_length ~ stype, data = df_sub_vivo)
  #coefficients_vivo <- summary(lm_vivo)$coefficients[,1][[2]]
  return(list(gene, coefficients_vitro, pv_vitro, coefficients_vivo, pv_vivo))
  
}

poly_a_d <- read.csv("ployA_with_gene_assigned_all_type.csv", header=TRUE)
#n_distinct(poly_a_d$gene_id)

gene_list <- unique(poly_a_d$gene_id)
raw_list = list()
count <- 0
for (gene in gene_list){
  t <- gene_lm_test(poly_a_d, gene)
  raw_list <- c(raw_list, t )
  count <- count +1
#  if (count >10){
#    break
# }
}
#raw_list
# #unlist(raw_list)
raw_matrix <- matrix(unlist(raw_list), nrow=5)
raw_df <- as.data.frame(t(raw_matrix))
colnames(raw_df) <- c("Gene","vitro_coefficient","vitro_pvalue","vivo_coefficient","vivo_pvalue")
#raw_df
write.csv(raw_df, "polyA_gene_test_2.csv")
