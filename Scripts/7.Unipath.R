library("UniPath", lib.loc= "/home/cell9/R/x86_64-pc-linux-gnu-library/4.1")

data("human_null_model") 

Pval_human <-binorm(human_null_data)

symbols <- read.csv("/home/haseena/endoderm/velocity/symbols.csv", row.names = 1) 

combine_human = combine(symbols, human_null_data, rownames (human_null_data), Pval_human, thr=5)

expres_count_matrix_data <- read.csv("/home/haseena/endoderm/velocity/endo_fpkm_matrix.csv", row.names=1)

Pval_count_matrix_data <- binorm(expres_count_matrix_data)

combine_exp_count_matrix <-combine(symbols, expres_count_matrix_data, rownames(expres_count_matrix_data), Pval_count_matrix_data, thr=5)

Exp_Uni_count_matrix <- adjust(combine_exp_count_matrix, combine_human)  

score = Exp_Uni_count_matrix[["adjpvalog"]]

write.csv(score,"0h_unipath.csv")
