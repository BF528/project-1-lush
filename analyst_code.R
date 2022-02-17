library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(abind)

# ----------------------- Functions  ------------------------

read_expression_table <- function(filename) {
  
  file <- read_csv(filename)
  
  cnames <- pull(file,1)
  
  file <- t(file)
  
  file <- file[-c(1),]
  
  colnames(file) <- cnames
  
  rownames(file) <- gsub('_.*',"",rownames(file))
  
  return (as_tibble(file,rownames=NA))
}

chi_sq <- function(data,sigma0,conf.level = 0.99){
  
  df = length(data) - 1
  chilower = qchisq((1 - conf.level)/2, df)
  chiupper = qchisq((1 - conf.level)/2, df, lower.tail = FALSE)
  v = var(data)
  testchi = df*v/(sigma0)
  testchi>chilower & testchi<chiupper
}

median_var <- function(x){
  median(
    apply(
      x, 
      MARGIN = 2, 
      var
    )
  )
}

eval_filters <- function(expr_mat) {
  
  expr_mat <- as.matrix(expr_mat)
  
  numeric_table <- matrix(as.numeric(expr_mat),
                          ncol = ncol(expr_mat))
  colnames(numeric_table) <- colnames(expr_mat)
  
  coefficient_of_variation_threshold <- apply(
    numeric_table,
    MARGIN = 2,
    function(x){sd(x)/mean(x)*100})>0.186
  
  coef_genes <- names(
    coefficient_of_variation_threshold[which(coefficient_of_variation_threshold==TRUE)])
  
  expression_threshold <- apply(
    numeric_table,
    MARGIN = 2,
    function(y){
      length(
        which(
          unlist(
            lapply(
              y, function(x){
                x>log2(15)}))==TRUE))/length(
                  unlist(
                    lapply(
                      y, 
                      function(x){
                        x>log2(15)})))*100}) > 20
  
  expr_genes <- names(expression_threshold[which(expression_threshold==TRUE)])
  
  sigma0 <- median_var(numeric_table)
  
  chi_threshold<-apply(
    numeric_table, 
    MARGIN = 2, 
    function(data){
      chi_sq(data,sigma0)
    })
  
  chi_genes <- names(chi_threshold[which(chi_threshold==TRUE)])
  
  passed_genes <- intersect(coef_genes,chi_genes)
  
  passed_genes <- intersect(passed_genes,expr_genes)
  
  return (passed_genes)
}

filter_matrix <- function(expr_mat) {
  
  passed_genes <- eval_filters(expr_mat)
  
  filtered_matrix <- expr_mat %>%
    
    select(all_of(passed_genes)) 
  
  
  filtered_matrix <- as.matrix(filtered_matrix)
  
  numeric_table <- matrix(as.numeric(filtered_matrix),
                          ncol = ncol(filtered_matrix))
  
  colnames(numeric_table) <- colnames(filtered_matrix)
  
  rownames(numeric_table) <- rownames(expr_mat)
  
  return (numeric_table)
}

# ----------------------- Script  ------------------------

filename <- ('projectnb/bf528/users/lush_2022/project_1/expression_data.csv')

expr_mat <- read_expression_table(filename)

filtered_matrix <- filter_matrix(expr_mat)

print(paste(length(colnames(filtered_matrix)),'genes passed filtering'))

write.csv(filtered_matrix,file = 'filtered_matrix.csv')

clusters <- hclust(dist(filtered_matrix))

plot(clusters)

clusterCut <- cutree(clusters, h=90)

print(paste(table(clusterCut),'samples in cluster'))

metadata <- read_csv('project/bf528/project_1/doc/proj_metadata.csv')

metadata <- metadata %>%
  
  select('geo_accession','cit-coloncancermolecularsubtype')

subtypes <- c()
for (i in rownames(filtered_matrix)) {
  idx <- which(metadata$geo_accession==i)
  subtypes <- c(subtypes,metadata$`cit-coloncancermolecularsubtype`[idx])
}

subtypes <- replace(subtypes,which(subtypes=="C3"),'red')

subtypes <- replace(subtypes,which(subtypes=="C4"),'blue')

heatmap(t(filtered_matrix), ColSideColors=subtypes)

cluster_idxs <- cbind(names(clusterCut),clusterCut)

group1 <- which(cluster_idxs[,2]==1)

group2 <- which(cluster_idxs[,2]==2)

t.test_results <- data.frame(filtered_matrix) %>%
  summarise(
    probe=colnames(filtered_matrix),
    t.statistic=apply(
      filtered_matrix, 2, function(x) { 
        t.test(x[group1], x[group2])$statistic }),
    p_val=apply(
      filtered_matrix, 2, function(x) {
        t.test(x[group1], x[group2])$p.value }),
    adjp_val=p.adjust(
      apply(
        filtered_matrix, 2, function(x) {
          t.test(x[group1], x[group2])$p.value }),
      method='fdr'))

write.csv(t.test_results,file = 't.test_table.csv')

print(paste('Genes with p adj < 0.05:', length(pull(filter(t.test_results,adjp_val<0.05),probe))))

sig_genes_res <- filter(t.test_results[order(t.test_results$adjp_val),],adjp_val<0.05)

sig_genes_mat <- as.matrix(sig_genes_res)

top_ten_percent <- sig_genes_mat[c(1:(length(rownames(sig_genes_res))*0.1)),]

cords <- c()
cluster_genes <- top_ten_percent[,1]
for (i in cluster_genes) {
  cords <- c(cords,which(colnames(filtered_matrix)==i))
}

cluster_genes

cluster_mat <- filtered_matrix[,cords]

heatmap(t(cluster_mat), ColSideColors=subtypes)