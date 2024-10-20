library(ggplot2)

pathway_score_all_cells <- read.csv(paste0("/data1/haseena_IP/pan_new/fastq/all_tsv_files/breast_blood_unipath_score_adjpvalogs.csv"), row.names = 1)    
#colnames(pathway_score_all_cells) <- gsub("_abundance.tsv")
colnames(pathway_score_all_cells) <- gsub("_abundance.tsv", "", colnames(pathway_score_all_cells))

selected_cellid <- gsub("converted_lofreq_", "", selected_cellid)
selected_cellid <- gsub(".hg38_multianno.txt", "", selected_cellid)
selected_cellid


rna_2_dna <- function(patient_id, id_rna) {
list_id_dna <- list()
for (i in 1:length(id_rna)) {
pos_id = which(cell_id_rna[, patient_id] %in% id_rna[i])
list_id_dna <- append(list_id_dna, cell_id_dna[pos_id, patient_id])
}
return(unlist(list_id_dna))
}


dna_2_rna <- function(patient_id, id_dna) {
list_id_rna <- list()
for (i in 1:length(id_dna)) {
pos_id = which(cell_id_dna[, patient_id] %in% id_dna[i])
list_id_rna <- append(list_id_rna, cell_id_rna[pos_id, patient_id])
}
return(unlist(list_id_rna))
}


setwd("/home/haseena")
cancer_hall_mark_data <- read.csv("/data1/haseena_IP/pan_new/fastq/all_tsv_files/GMT_Breast_cancer.csv")

wilcoxson_function <- function(data, group){
clases = unique(group)
colnames(data) <- group
group1 = data[,which(colnames(data) %in% clases[1])]
group2 = data[,which(colnames(data) %in% clases[2])]
wil_mat = matrix(1, nrow(data) , 1) ;
FC_mat = wil_mat

pos = which(group == "CTC") ;
pos1 = which(group == "non-CTC") ;

for (k in 1:length(data[,1])) {
output = wilcox.test( data[k, pos] , data[k, pos1] ) 
wil_mat[k,1] = output$p.value
FC = median(data[k,pos]) - median(data[k, pos1])
#FC_mat[k,1] = sign(FC) * -log2(abs(FC))
FC_mat[k,1] = FC
}

final_mat = as.data.frame(cbind(wil_mat, FC_mat))
colnames(final_mat) <- c('wilx_pvalue', 'fc')
rownames(final_mat) <- rownames(data)
final_mat[is.na(final_mat)] <- 0
return(final_mat)
}


plot_volcano <- function(feature, feature_non, patient_id, fc_threshold) {
name = patient_id
plot_name <- paste("/home/haseena/violin_plot_enrich_", name, ".pdf", sep = "")
title_name <- paste('Unipath results: ', toupper(name), sep = "")

#rownames(feature) <- cancer_hall_mark_data[,1]
#rownames(feature_non) <- cancer_hall_mark_data[,1]

data1 <- as.data.frame(cbind(feature, feature_non))
data <- apply(data1, 2, function(x) as.numeric(x))
rownames(data) <- rownames(data1) 

group <- c(rep('CTC', ncol(feature)), rep('non-CTC',  ncol(feature_non)))

group = as.matrix(group)
data = as.matrix(data)

wilx =  wilcoxson_function(data, group)
wilx$wilx_pvalue <- -log10(wilx$wilx_pvalue)

wilx$diffexpressed <- "NO"
wilx$diffexpressed[wilx$fc > 0 & wilx$wilx_pvalue > 1] <- "UP"
wilx$diffexpressed[wilx$fc < 0 & wilx$wilx_pvalue > 1] <- "DOWN"

labels_pre = rownames(wilx)
labels_pre[wilx$fc == 0 ] = ""
labels_pre[wilx$fc > -fc_threshold  & wilx$fc < fc_threshold] = ""

wilx = as.data.frame(wilx)
p <- ggplot(data=as.data.frame(wilx), aes(x=fc, y=wilx_pvalue, col=diffexpressed)) +
geom_point() +
theme_minimal() + 
#geom_text(label=labels_pre, size = 2) +
ggrepel::geom_text_repel(
    aes(label=labels_pre), size=2, box.padding = unit(0.5, "lines")
  )

p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
#p2 <- p + geom_vline(xintercept=c(-threshold_fc, threshold_fc), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red") 

#p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)

pdf(plot_name, width = 10, height = 10)
plot(p3)
dev.off()
system(paste('code', plot_name))

}

filter_cells_on_pathway_score <- function(pathway_scores) {

zero_freq_mat = cbind(colnames(pathway_scores), NA) 
rownames(zero_freq_mat) = zero_freq_mat[,1]
unlist_scores = unlist(as.data.frame(pathway_scores))
unlist_scores = unlist_scores[which(unlist_scores != 0)]

for (col_i in 1:ncol(pathway_scores)) {
  # print(col_i)
zero_freq_mat[colnames(pathway_scores)[col_i], 2] = length(which(pathway_scores[,col_i] >= quantile(unlist_scores, 0.75)))
#zero_freq_mat[colnames(pathway_scores)[col_i], 2] = length(which(pathway_scores[,col_i] >= 0.9))

}
rownames(zero_freq_mat) = NULL 
zero_freq_mat = as.data.frame(as.matrix(zero_freq_mat))
zero_freq_mat[,2] <- as.numeric(zero_freq_mat[,2])
higher_count_cell_id = zero_freq_mat[which(zero_freq_mat$V2 >= quantile(zero_freq_mat$V2, 0.8)),"V1"]

return(pathway_scores[,higher_count_cell_id])

}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#/data1/haseena_IP/pan_new/ALL_annovar_output/ALL_RANK.csv
#/data1/haseena_IP/pan_new/new_pan_ctc/filter_ctc_bam/annovar/CTC_RANK.csv
#/data1/haseena_IP/pan_new/new_pan_wbc/filtered_wbc_bam/converted/annovar/WBC_RANK.csv
#/data1/haseena_IP/pan_new/new_pan_xeno/filter_xeno_bam/converted/annovar/XENO_RANK.csv
selected_g50 <- read.csv("/data1/haseena_IP/pan_new/ALL_annovar_output/ALL_RANK.csv")
#selected_g51 <- read.csv("/data1/haseena_IP/pan_new/new_pan_ctc/filter_ctc_bam/annovar/CTC_RANK.csv")
#selected_g52 <- read.csv("/data1/haseena_IP/pan_new/new_pan_wbc/filtered_wbc_bam/converted/annovar/WBC_RANK.csv")
#selected_g53 <- read.csv("/data1/haseena_IP/pan_new/new_pan_xeno/filter_xeno_bam/converted/annovar/XENO_RANK.csv")

selected_g50_cell_id_pre <- gsub("converted_lofreq_", "", selected_g50$cell_id)
selected_g50$cell_id <- gsub(".hg38_multianno.txt", "", selected_g50_cell_id_pre)
 
 
#selected_g51_cell_id_pre <- gsub("converted_lofreq_", "", selected_g51$cell_id)
#selected_g51$cell_id <- gsub(".hg38_multianno.txt", "", selected_g51_cell_id_pre)


#selected_g52_cell_id_pre <- gsub("converted_lofreq_", "", selected_g52$cell_id)
#selected_g52$cell_id <- gsub(".hg38_multianno.txt", "", selected_g52_cell_id_pre)


#selected_g53_cell_id_pre <- gsub("converted_lofreq_", "", selected_g53$cell_id)
#selected_g53$cell_id <- gsub(".hg38_multianno.txt", "", selected_g53_cell_id_pre)

# blood_marks <- c("HALLMARK_COAGULATION", "BLOOD_DN", "BLOOD_UP", "RED_BLOOD_CELL", "WHITE BLOOD CELL:BLOOD", "LYMPHOCYTES_RNA_SEQ_UP", 
# "LYMPHOCYTES_T", "CD4_T_CELLS", "B_CELLS", "CD8_T_CELLS", "MONOCYTE_CELLS", "NK_CELLS", "MEMORY_T_CELL", "MEMORY_B_CELL", "NAIVE_B_CELL", "NAIVE_T_CELL",
# "NEUTROPHIL", "EOSINOPHIL"
# "ZHENG_CORD_BLOOD_C10_MULTILYMPHOID_PROGENITOR",
# "ZHENG_CORD_BLOOD_C1_PUTATIVE_MEGAKARYOCYTE_PROGENITOR", "ZHENG_CORD_BLOOD_C2_PUTATIVE_BASOPHIL_EOSINOPHIL_MAST_CELL_PROGENITOR", 
# "ZHENG_CORD_BLOOD_C3_MEGAKARYOCYTE_ERYTHROID_PROGENITOR", "ZHENG_CORD_BLOOD_C4_PUTATIVE_EARLY_ERYTHROID_COMMITMENT", 
# "ZHENG_CORD_BLOOD_C5_SIMILAR_TO_HSC_C6_PUTATIVE_ALTERED_METABOLIC_STATE", "ZHENG_CORD_BLOOD_C6_HSC_MULTIPOTENT_PROGENITOR",
# "ZHENG_CORD_BLOOD_C7_PUTATIVE_LYMPHOID_PRIMED_MULTIPOTENT_PROGENITOR_1", "ZHENG_CORD_BLOOD_C8_PUTATIVE_LYMPHOID_PRIMED_MULTIPOTENT_PROGENITOR_2",
# "ZHENG_CORD_BLOOD_C9_GRANULOCYTE_MACROPHAGE_PROGENITOR")

blood_marks <- c("CD4_T_CELLS", "CD8_T_CELLS", "B_CELLS", "NK_CELLS", "MEMORY_B_CELL", "MEMORY_T_CELL", "NEUTROPHIL", "EOSINOPHIL",
"NAIVE_B_CELL", "NAIVE_T_CELL")

breast_pathways <- rownames(pathway_score_pre)[grep("BREAST", rownames(pathway_score_pre), ignore.case = TRUE)]

select_positives_negatives_for_ml <- function(selected_g, patient_id, num_cells_percent, fc_threshold) {

pathway_score_pre <- pathway_score_all_cells[,which(colnames(pathway_score_all_cells) %in% selected_g$cell_id)]  # select the pathway scores for the cohort
pathway_scores = preprocessCore::normalize.quantiles(as.matrix(pathway_score_pre)) # qunatile normalize the cells
colnames(pathway_scores) <- colnames(pathway_score_pre)
rownames(pathway_scores) <- rownames(pathway_score_pre)

cell_id_to_check <- selected_g$cell_id[1:(nrow(selected_g) * num_cells_percent)] # selected top ranking cells
breast_pathway <- pathway_scores[breast_pathways, cell_id_to_check]  # check breastse related enrichment in 
row_sums <- rowSums(breast_pathway)
row_sums <- row_sums[order(row_sums, decreasing = TRUE)]

# pathway_scores = pathway_score_pre
breast_enriched <- breast_pathway[order(rowSums(breast_pathway), decreasing = TRUE),]
breast_enriched_rowsums <- cbind(breast_enriched, row_sums)
breast_enriched_rowsums = breast_enriched_rowsums[which(breast_enriched_rowsums[,"row_sums"] > 0),]

pathway_scores_filter = filter_cells_on_pathway_score(pathway_scores)
message("No. of cells filtered out: ", ncol(pathway_scores) - ncol(pathway_scores_filter), " out of ", ncol(pathway_scores) )
cor_pathway <- cor(t(pathway_scores_filter))

anti_cor_pathways <- list()

for (i in rownames(breast_enriched_rowsums)) {
print(i)
anti_cor_pathways <- append(anti_cor_pathways, names(which(cor_pathway[i,] < -0.1)) )
}

anti_cor_pathway <- as.data.frame(table(unlist(anti_cor_pathways)))

#anti_cor_pathway_non_pan <- anti_cor_pathway[-grep("breast", anti_cor_pathway$Var1),]
anti_cor_pathway_non_pan <- anti_cor_pathway[which(anti_cor_pathway$Var1 %in% blood_marks),]

anti_cor_pathway_non_pan <- anti_cor_pathway_non_pan[order(anti_cor_pathway_non_pan$Freq, decreasing = TRUE), ]
message("selecting negatives having enrichment of following anti-correlated pathways: ")
print(as.character(anti_cor_pathway_non_pan$Var1))

pos_neg_check <- which(colnames(pathway_scores_filter) %in% cell_id_to_check)
if (length(pos_neg_check) != 0) {pathway_scores_filter <- pathway_scores_filter[,-which(colnames(pathway_scores_filter) %in% cell_id_to_check)]}

enrich_count_mat <- as.data.frame(matrix(NA, ncol(pathway_scores_filter), 2)); colnames(enrich_count_mat) <- c("cell_id", "num_of_enrich")

for (col in 1:ncol(pathway_scores_filter)) {
enrich_count_mat[col, "cell_id"] <- colnames(pathway_scores_filter)[col] 
enrich_count_mat[col, "num_of_enrich"] <- length(which(pathway_scores_filter[anti_cor_pathway_non_pan$Var1, col] > 0.1)) # check by varing
}

#enrichment filter
enrich_count_mat = enrich_count_mat[order(enrich_count_mat$num_of_enrich, decreasing = TRUE),]
rna_id_ctc_non_pre = enrich_count_mat[which(enrich_count_mat$num_of_enrich >= quantile(enrich_count_mat$num_of_enrich, 0.7)), "cell_id"]
message("No. of non-CTCs after anti-corr pathway enrich : ", length(rna_id_ctc_non_pre))


#mutation filter
dna_id_ctc_non_mut_count <- selected_g[which(selected_g$cell_id %in% rna_id_ctc_non_pre),c("cell_id", "no_somatic_mutations")]
dna_id_ctc_non_mut_count = dna_id_ctc_non_mut_count[order(dna_id_ctc_non_mut_count$no_somatic_mutations),]
rna_id_ctc_non = dna_id_ctc_non_mut_count[which(dna_id_ctc_non_mut_count$no_somatic_mutations <= quantile(dna_id_ctc_non_mut_count$no_somatic_mutations, 0.1)), "cell_id"] # selection of negatives based on the number of somatic mutaions
message("No. of non-CTCs after mutation filter: ", length(rna_id_ctc_non))


df_ids <- gdata::cbindX(as.data.frame( selected_g$cell_id[1:round(nrow(selected_g) * num_cells_percent)]), as.data.frame(rna_id_ctc_non))

df_ids[is.na(df_ids)] <- ""
colnames(df_ids) <- c("CTC", "non-CTC")



if (length(rna_id_ctc_non) > 10) {set.seed(1); rna_id_ctc_non = sample(rna_id_ctc_non, 10)} else {rna_id_ctc_non = rna_id_ctc_non}

message("num of non-CTCs: ", length(rna_id_ctc_non))
non_ctc_score <- pathway_score_all_cells[, rna_id_ctc_non]
non_ctc_score [non_ctc_score < 0] = 0
non_ctc_score <- as.matrix(non_ctc_score)
ctc_score <- pathway_score_all_cells[, df_ids$CTC[which(df_ids$CTC != "") ]]

message("CTC number: ", length(which(df_ids[,1] != "")))
message("non-CTC number: ", length(which(df_ids[,2] != "")))

plot_volcano(ctc_score, non_ctc_score, patient_id, fc_threshold)

message("no error")


return(df_ids)
}

# arguments are selected_g, patient_id, num_cells_percent, fc_threshold 
train_ml_g50 <- select_positives_negatives_for_ml(selected_g50, "g50", 0.05, 0.01) # pass
#train_ml_g51 <- select_positives_negatives_for_ml(selected_g51, "g51", 0.08, 0.01) # pass
#train_ml_g52 <- select_positives_negatives_for_ml(selected_g52, "g52", 0.08, 0.01) # pass
#train_ml_g53 <- select_positives_negatives_for_ml(selected_g53, "g53", 0.08, 0.01) # pass


pathway_g50 <- read.csv(paste0("data/g50_unipath_score_pl1_norm.csv"), row.names = 1)    

ctc_id_rna_g50 = ctc_id_rna_g50_pre[which(ctc_id_rna_g50_pre != "")]
ctc_non_id_rna_g50 = ctc_non_id_rna_g50_pre[which(ctc_non_id_rna_g50_pre != "")]
ctc_score_g50 <- pathway_g50[, ctc_id_rna_g50]
ctc_non_score_g50 <- pathway_g50[, ctc_non_id_rna_g50]

#pathway_g51 <- read.csv(paste0("data/g51_unipath_score_pl1_norm.csv"), row.names = 1)    

#ctc_id_rna_g51_pre <- dna_2_rna("G51", train_ml_g51$CTC)
#ctc_id_rna_g51 = ctc_id_rna_g51_pre[which(ctc_id_rna_g51_pre != "")]
#ctc_non_id_rna_g51_pre <- dna_2_rna("G51", train_ml_g51$non_CTC)
#ctc_non_id_rna_g51 = ctc_non_id_rna_g51_pre[which(ctc_non_id_rna_g51_pre != "")]
#ctc_score_g51 <- pathway_g51[, ctc_id_rna_g51]
#ctc_non_score_g51 <- pathway_g51[, ctc_non_id_rna_g51]

#pathway_g52 <- read.csv(paste0("data/g52_unipath_score_pl1_norm.csv"), row.names = 1)    

#ctc_id_rna_g52_pre <- dna_2_rna("G52", train_ml_g52$CTC)
#ctc_id_rna_g52 = ctc_id_rna_g52_pre[which(ctc_id_rna_g52_pre != "")]
#ctc_non_id_rna_g52_pre <- dna_2_rna("G52", train_ml_g52$non_CTC)
#ctc_non_id_rna_g52 = ctc_non_id_rna_g52_pre[which(ctc_non_id_rna_g52_pre != "")]
#ctc_score_g52 <- pathway_g52[, ctc_id_rna_g52]
#ctc_non_score_g52 <- pathway_g52[, ctc_non_id_rna_g52]

#pathway_g53 <- read.csv(paste0("data/g53_unipath_score_pl1_norm.csv"), row.names = 1)    

#ctc_id_rna_g53_pre <- dna_2_rna("G53", train_ml_g53$CTC)
#ctc_id_rna_g53 = ctc_id_rna_g53_pre[which(ctc_id_rna_g53_pre != "")]
#ctc_non_id_rna_g53_pre <- dna_2_rna("G53", train_ml_g53$non_CTC)
#ctc_non_id_rna_g53 = ctc_non_id_rna_g53_pre[which(ctc_non_id_rna_g53_pre != "")]
#ctc_score_g53 <- pathway_g53[, ctc_id_rna_g53]
#ctc_non_score_g53 <- pathway_g53[, ctc_non_id_rna_g53]

#pathway_g54 <- read.csv(paste0("data/g54_unipath_score_pl1_norm.csv"), row.names = 1)    
#ctc_id_rna_g54_pre <- dna_2_rna("G54", train_ml_g54$CTC)
#ctc_id_rna_g54 = ctc_id_rna_g54_pre[which(ctc_id_rna_g54_pre != "")]
#ctc_non_id_rna_g54_pre <- dna_2_rna("G54", train_ml_g54$non_CTC)
#ctc_non_id_rna_g54 = ctc_non_id_rna_g54_pre[which(ctc_non_id_rna_g54_pre != "")]
#ctc_score_g54 <- pathway_g54[, ctc_id_rna_g54]
#ctc_non_score_g54 <- pathway_g54[, ctc_non_id_rna_g54]

# pathway_g56 <- read.csv(paste0("data/g56_unipath_score_pl1_norm.csv"), row.names = 1)    
# ctc_id_rna_g56_pre <- dna_2_rna("G56", train_ml_g56$CTC)
# ctc_id_rna_g56 = ctc_id_rna_g56_pre[which(ctc_id_rna_g56_pre != "")]
# ctc_non_id_rna_g56_pre <- dna_2_rna("G56", train_ml_g56$non_CTC)
# ctc_non_id_rna_g56 = ctc_non_id_rna_g56_pre[which(ctc_non_id_rna_g56_pre != "")]
# ctc_score_g56 <- pathway_g56[, ctc_id_rna_g56]
# ctc_non_score_g56 <- pathway_g56[, ctc_non_id_rna_g56]

#pathway_g57 <- read.csv(paste0("data/g57_unipath_score_pl1_norm.csv"), row.names = 1)    
#ctc_id_rna_g57_pre <- dna_2_rna("G57", train_ml_g57$CTC)
#ctc_id_rna_g57 = ctc_id_rna_g57_pre[which(ctc_id_rna_g57_pre != "")]
#ctc_non_id_rna_g57_pre <- dna_2_rna("G57", train_ml_g57$non_CTC)
#ctc_non_id_rna_g57 = ctc_non_id_rna_g57_pre[which(ctc_non_id_rna_g57_pre != "")]
#ctc_score_g57 <- pathway_g57[, ctc_id_rna_g57]
#ctc_non_score_g57 <- pathway_g57[, ctc_non_id_rna_g57]

#pathway_g59 <- read.csv(paste0("data/g59_unipath_score_pl1_norm.csv"), row.names = 1)    
#ctc_id_rna_g59_pre <- dna_2_rna("G59", train_ml_g59$CTC)
#ctc_id_rna_g59 = ctc_id_rna_g59_pre[which(ctc_id_rna_g59_pre != "")]
#ctc_non_id_rna_g59_pre <- dna_2_rna("G59", train_ml_g59$non_CTC)
#ctc_non_id_rna_g59 = ctc_non_id_rna_g59_pre[which(ctc_non_id_rna_g59_pre != "")]
#ctc_score_g59 <- pathway_g59[, ctc_id_rna_g59]
#ctc_non_score_g59 <- pathway_g59[, ctc_non_id_rna_g59]

#pathway_g60 <- read.csv(paste0("data/g60_unipath_score_pl1_norm.csv"), row.names = 1)    
#ctc_id_rna_g60_pre <- dna_2_rna("G60", train_ml_g60$CTC)
#ctc_id_rna_g60 = ctc_id_rna_g60_pre[which(ctc_id_rna_g60_pre != "")]
#ctc_non_id_rna_g60_pre <- dna_2_rna("G60", train_ml_g60$non_CTC)
#ctc_non_id_rna_g60 = ctc_non_id_rna_g60_pre[which(ctc_non_id_rna_g60_pre != "")]
#ctc_score_g60 <- pathway_g60[, ctc_id_rna_g60]
#ctc_non_score_g60 <- pathway_g60[, ctc_non_id_rna_g60]

#pathway_hth5 <- read.csv(paste0("data/hth5_unipath_score_pl1_norm.csv"), row.names = 1)    
#ctc_id_rna_hth5_pre <- dna_2_rna("HTH5", train_ml_hth5$CTC)
#ctc_id_rna_hth5 = ctc_id_rna_hth5_pre[which(ctc_id_rna_hth5_pre != "")]
#ctc_non_id_rna_hth5_pre <- dna_2_rna("HTH5", train_ml_hth5$non_CTC)
#ctc_non_id_rna_hth5 = ctc_non_id_rna_hth5_pre[which(ctc_non_id_rna_hth5_pre != "")]
#ctc_score_hth5 <- pathway_hth5[, ctc_id_rna_hth5]
#ctc_non_score_hth5 <- pathway_hth5[, ctc_non_id_rna_hth5]


ctc_pathway_all <- cbind(ctc_score_g50)
#ctc_pathway_all <- cbind(ctc_score_g50, ctc_score_g51, ctc_score_g52, ctc_score_g54,
#ctc_score_g57, ctc_score_g59, ctc_score_g60, ctc_score_hth5)
ctc_non_pathway_all <- cbind(ctc_non_score_g50)
#ctc_non_pathway_all <- cbind(ctc_non_score_g50, ctc_non_score_g51, ctc_non_score_g52, ctc_non_score_g54,
#ctc_non_score_g57, ctc_non_score_g59, ctc_non_score_g60, ctc_non_score_hth5)

plot_volcano(ctc_pathway_all, ctc_non_pathway_all, "all_patients", 0.01)


func_change_id_df <- function(patient_id, train_ml_g) {
df_if_ctc_non_g <- gdata::cbindX(as.data.frame(dna_2_rna(patient_id, train_ml_g$CTC)), as.data.frame(dna_2_rna(patient_id, train_ml_g$non_CTC)))
df_if_ctc_non_g[is.na(df_if_ctc_non_g)] <- ""
colnames(df_if_ctc_non_g) <- c("rna_id_ctc", "rna_id_ctc_non")
df_if_ctc_non_g = df_if_ctc_non_g[!apply(df_if_ctc_non_g == "", 1, all),]
return(df_if_ctc_non_g)
}

df_if_ctc_non_g50 = func_change_id_df("G50", train_ml_g50)
#df_if_ctc_non_g51 = func_change_id_df("G51", train_ml_g51)
#df_if_ctc_non_g52 = func_change_id_df("G52", train_ml_g52)
#df_if_ctc_non_g53 = func_change_id_df("G53", train_ml_g53)

write.csv(df_if_ctc_non_g50, "data/df_id_ctc_and_non_g50.csv")
#write.csv(df_if_ctc_non_g51, "data/df_id_ctc_and_non_g51.csv")
#write.csv(df_if_ctc_non_g52, "data/df_id_ctc_and_non_g52.csv")
#write.csv(df_if_ctc_non_g53, "data/df_id_ctc_and_non_g53.csv")