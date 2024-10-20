library(ggplot2)
library(preprocessCore)
library(ggrepel)
library(gdata)

# Define ranges for fine-tuning
quantile_range <- seq(0.6, 0.9, by = 0.05)
higher_count_cell_id_range <- seq(0.6, 0.9, by = 0.05)
anti_cor_value_range <- seq(-0.05, -0.13, by = -0.02)
enrich_threshold_range <- seq(0.05, 0.13, by = 0.02)
mutation_quantile_range <- seq(0.1, 0.3, by = 0.05)

# Load data
pathway_score_all_cells <- read.csv("/data1/haseena_IP/pan_new/fastq/all_tsv_files/breast_blood_unipath_score_adjpvalogs.csv", row.names = 1)
colnames(pathway_score_all_cells) <- gsub("_abundance.tsv", "", colnames(pathway_score_all_cells))

# Load selected cell data
selected_g50 <- read.csv("/data1/haseena_IP/pan_new/ALL_annovar_output/ALL_RANK.csv")
selected_g51 <- read.csv("/data1/haseena_IP/pan_new/new_pan_ctc/filter_ctc_bam/annovar/CTC_RANK.csv")
selected_g52 <- read.csv("/data1/haseena_IP/pan_new/new_pan_wbc/filtered_wbc_bam/converted/annovar/WBC_RANK.csv")
selected_g53 <- read.csv("/data1/haseena_IP/pan_new/new_pan_xeno/filter_xeno_bam/converted/annovar/XENO_RANK.csv")

# Clean cell IDs
selected_g50_cell_id_pre <- gsub("converted_lofreq_", "", selected_g50$cell_id)
selected_g50$cell_id <- gsub(".hg38_multianno.txt", "", selected_g50_cell_id_pre)

selected_g51_cell_id_pre <- gsub("converted_lofreq_", "", selected_g51$cell_id)
selected_g51$cell_id <- gsub(".hg38_multianno.txt", "", selected_g51_cell_id_pre)

selected_g52_cell_id_pre <- gsub("converted_lofreq_", "", selected_g52$cell_id)
selected_g52$cell_id <- gsub(".hg38_multianno.txt", "", selected_g52_cell_id_pre)

selected_g53_cell_id_pre <- gsub("converted_lofreq_", "", selected_g53$cell_id)
selected_g53$cell_id <- gsub(".hg38_multianno.txt", "", selected_g53_cell_id_pre)

# Load cancer hallmark data
setwd("/home/haseena")
cancer_hall_mark_data <- read.csv("/data1/haseena_IP/pan_new/fastq/all_tsv_files/GMT_Breast_cancer.csv")

# Blood markers and breast pathways
blood_marks <- c("CD4_T_CELLS", "CD8_T_CELLS", "B_CELLS", "NK_CELLS", "MEMORY_B_CELL", "MEMORY_T_CELL", "NEUTROPHIL", "EOSINOPHIL", "NAIVE_B_CELL", "NAIVE_T_CELL")

# Extract breast pathways
breast_pathways <- rownames(pathway_score_all_cells)[grep("BREAST", rownames(pathway_score_all_cells), ignore.case = TRUE)]

# Function for Wilcoxon rank-sum test
wilcoxson_function <- function(data, group) {
    clases = unique(group)
    colnames(data) <- group
    group1 = data[, which(colnames(data) %in% clases[1])]
    group2 = data[, which(colnames(data) %in% clases[2])]
    wil_mat = matrix(1, nrow(data), 1)
    FC_mat = wil_mat

    pos = which(group == "CTC")
    pos1 = which(group == "non-CTC")

    for (k in 1:length(data[, 1])) {
        # Convert to numeric and check for NAs
        data1 <- as.numeric(data[k, pos])
        data2 <- as.numeric(data[k, pos1])

        # Identify and print NA values
        if (any(is.na(data1)) || any(is.na(data2))) {
            cat("NA values found in row:", k, "\n")
            cat("Data1 (CTC):", data1, "\n")
            cat("Data2 (non-CTC):", data2, "\n")
        }

        # Only perform the test if there are no NAs
        if (!any(is.na(data1)) && !any(is.na(data2))) {
            output = wilcox.test(data1, data2)
            wil_mat[k, 1] = output$p.value
            FC = median(data1) - median(data2)
            FC_mat[k, 1] = FC
        } else {
            wil_mat[k, 1] = NA
            FC_mat[k, 1] = NA
        }
    }

    final_mat = as.data.frame(cbind(wil_mat, FC_mat))
    colnames(final_mat) <- c('wilx_pvalue', 'fc')
    rownames(final_mat) <- rownames(data)
    final_mat[is.na(final_mat)] <- 0
    return(final_mat)
}

# Function to plot volcano
plot_volcano <- function(feature, feature_non, patient_id, fc_threshold) {
    name = patient_id
    plot_name <- paste("/home/haseena/volcano/violin_plot_enrich_", name, ".pdf", sep = "")
    title_name <- paste('Unipath results: ', toupper(name), sep = "")

    data1 <- as.data.frame(cbind(feature, feature_non))
    data <- apply(data1, 2, function(x) as.numeric(x))
    rownames(data) <- rownames(data1) 

    group <- c(rep('CTC', ncol(feature)), rep('non-CTC', ncol(feature_non)))

    group = as.matrix(group)
    data = as.matrix(data)

    wilx = wilcoxson_function(data, group)
    wilx$wilx_pvalue <- -log10(wilx$wilx_pvalue)

    wilx$diffexpressed <- "NO"
    wilx$diffexpressed[wilx$fc > 0 & wilx$wilx_pvalue > 1] <- "UP"
    wilx$diffexpressed[wilx$fc < 0 & wilx$wilx_pvalue > 1] <- "DOWN"

    labels_pre = rownames(wilx)
    labels_pre[wilx$fc == 0] = ""
    labels_pre[wilx$fc > -fc_threshold & wilx$fc < fc_threshold] = ""

    wilx = as.data.frame(wilx)
    p <- ggplot(data=as.data.frame(wilx), aes(x=fc, y=wilx_pvalue, col=diffexpressed)) +
        geom_point() +
        theme_minimal() + 
        ggrepel::geom_text_repel(aes(label=labels_pre), size=2, box.padding = unit(0.5, "lines"))

    p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red") 

    mycolors <- c("blue", "red", "black")
    names(mycolors) <- c("DOWN", "UP", "NO")
    p3 <- p2 + scale_colour_manual(values = mycolors)

    pdf(plot_name, width = 10, height = 10)
    plot(p3)
    dev.off()
    system(paste('code', plot_name))
}

# Function to filter cells based on pathway scores
filter_cells_on_pathway_score <- function(pathway_scores, quantile_thresh, count_thresh) {
    zero_freq_mat = cbind(colnames(pathway_scores), NA) 
    rownames(zero_freq_mat) = zero_freq_mat[, 1]
    unlist_scores = unlist(as.data.frame(pathway_scores))
    unlist_scores = unlist_scores[which(unlist_scores != 0)]
    
    for (col_i in 1:ncol(pathway_scores)) {
        zero_freq_mat[colnames(pathway_scores)[col_i], 2] = length(which(pathway_scores[, col_i] >= quantile(unlist_scores, quantile_thresh)))
    }
    rownames(zero_freq_mat) = NULL 
    zero_freq_mat = as.data.frame(as.matrix(zero_freq_mat))
    zero_freq_mat[, 2] <- as.numeric(zero_freq_mat[, 2])
    higher_count_cell_id = zero_freq_mat[which(zero_freq_mat$V2 >= quantile(zero_freq_mat$V2, count_thresh)), "V1"]

    return(pathway_scores[, higher_count_cell_id])
}

# Function to select positives and negatives for ML
select_positives_negatives_for_ml <- function(selected_g, patient_id, num_cells_percent, fc_threshold, quantile_thresh, count_thresh, anti_cor_value, enrich_threshold, mutation_quantile) {
    pathway_score_pre <- pathway_score_all_cells[, which(colnames(pathway_score_all_cells) %in% selected_g$cell_id)]
    pathway_scores = preprocessCore::normalize.quantiles(as.matrix(pathway_score_pre))
    colnames(pathway_scores) <- colnames(pathway_score_pre)
    rownames(pathway_scores) <- rownames(pathway_score_pre)

    cell_id_to_check <- selected_g$cell_id[1:(nrow(selected_g) * num_cells_percent)]
    breast_pathway <- pathway_scores[breast_pathways, cell_id_to_check]
    row_sums <- rowSums(breast_pathway)
    row_sums <- row_sums[order(row_sums, decreasing = TRUE)]

    breast_enriched <- breast_pathway[order(rowSums(breast_pathway), decreasing = TRUE),]
    breast_enriched_rowsums <- cbind(breast_enriched, row_sums) 
    breast_enriched_rowsums = breast_enriched_rowsums[which(breast_enriched_rowsums[, "row_sums"] > 0),]

    pathway_scores_filter = filter_cells_on_pathway_score(pathway_scores, quantile_thresh, count_thresh)

    # Message for number of cells filtered out
    cat("No. of cells filtered out: ", ncol(pathway_scores) - ncol(pathway_scores_filter), " out of ", ncol(pathway_scores), "\n")
    
    cor_pathway <- cor(t(pathway_scores_filter))

    anti_cor_pathways <- list()

    for (i in rownames(breast_enriched_rowsums)) {
        anti_cor_pathways <- append(anti_cor_pathways, names(which(cor_pathway[i, ] < anti_cor_value)))
    }

    anti_cor_pathway <- as.data.frame(table(unlist(anti_cor_pathways)))
    anti_cor_pathway_non_pan <- anti_cor_pathway[which(anti_cor_pathway$Var1 %in% blood_marks),]

    anti_cor_pathway_non_pan <- anti_cor_pathway_non_pan[order(anti_cor_pathway_non_pan$Freq, decreasing = TRUE),]
    
    # Message for selected negatives
    cat("selecting negatives having enrichment of following anti-correlated pathways: \n")
    cat(as.character(anti_cor_pathway_non_pan$Var1), "\n")

    pos_neg_check <- which(colnames(pathway_scores_filter) %in% cell_id_to_check)
    if (length(pos_neg_check) != 0) {
        pathway_scores_filter <- pathway_scores_filter[, -which(colnames(pathway_scores_filter) %in% cell_id_to_check)]
    }

    enrich_count_mat <- as.data.frame(matrix(NA, ncol(pathway_scores_filter), 2))
    colnames(enrich_count_mat) <- c("cell_id", "num_of_enrich")

    for (col in 1:ncol(pathway_scores_filter)) {
        enrich_count_mat[col, "cell_id"] <- colnames(pathway_scores_filter)[col] 
        enrich_count_mat[col, "num_of_enrich"] <- length(which(pathway_scores_filter[anti_cor_pathway_non_pan$Var1, col] > enrich_threshold))
    }

    enrich_count_mat = enrich_count_mat[order(enrich_count_mat$num_of_enrich, decreasing = TRUE),]
    rna_id_ctc_non_pre = enrich_count_mat[which(enrich_count_mat$num_of_enrich >= quantile(enrich_count_mat$num_of_enrich, 0.9)), "cell_id"]

    # Message for non-CTCs after anti-correlation filtering
    cat("No. of non-CTCs after anti-corr pathway enrich: ", length(rna_id_ctc_non_pre), "\n")

    dna_id_ctc_non_mut_count <- selected_g[which(selected_g$cell_id %in% rna_id_ctc_non_pre), c("cell_id", "no_somatic_mutations")]
    dna_id_ctc_non_mut_count = dna_id_ctc_non_mut_count[order(dna_id_ctc_non_mut_count$no_somatic_mutations),]
    rna_id_ctc_non = dna_id_ctc_non_mut_count[which(dna_id_ctc_non_mut_count$no_somatic_mutations <= quantile(dna_id_ctc_non_mut_count$no_somatic_mutations, mutation_quantile)), "cell_id"]
    
    # Message for non-CTCs after mutation filtering
    cat("No. of non-CTCs after mutation filter: ", length(rna_id_ctc_non), "\n")

    df_ids <- gdata::cbindX(as.data.frame(selected_g$cell_id[1:round(nrow(selected_g) * num_cells_percent)]), as.data.frame(rna_id_ctc_non))
    df_ids[is.na(df_ids)] <- ""
    colnames(df_ids) <- c("CTC", "non-CTC")

    if (length(rna_id_ctc_non) > 10) {
        set.seed(1)
        rna_id_ctc_non = sample(rna_id_ctc_non, 10)
    } else {
        rna_id_ctc_non = rna_id_ctc_non
    }

    # Message for number of non-CTCs
    cat("num of non-CTCs: ", length(rna_id_ctc_non), "\n")
    
    non_ctc_score <- pathway_score_all_cells[, rna_id_ctc_non]
    non_ctc_score[non_ctc_score < 0] = 0
    non_ctc_score <- as.matrix(non_ctc_score)
    ctc_score <- pathway_score_all_cells[, df_ids$CTC[which(df_ids$CTC != "")]]

    # Message for number of CTCs and non-CTCs in final dataset
    cat("CTC number: ", length(which(df_ids[, 1] != "")), "\n")
    cat("non-CTC number: ", length(which(df_ids[, 2] != "")), "\n")

    # Ensure that the scores are numeric
    ctc_score <- apply(ctc_score, 2, as.numeric)
    non_ctc_score <- apply(non_ctc_score, 2, as.numeric)

    # Return scores for further evaluation
    return(list(ctc_score = ctc_score, non_ctc_score = non_ctc_score, df_ids = df_ids))
}

# Function to score the parameter combination
# Function to score the parameter combination
score_parameters <- function(ctc_score, non_ctc_score, fc_threshold) {
    wilx = wilcoxson_function(as.data.frame(cbind(ctc_score, non_ctc_score)), 
                              c(rep('CTC', ncol(ctc_score)), rep('non-CTC', ncol(non_ctc_score))))
    wilx$wilx_pvalue <- -log10(wilx$wilx_pvalue)

    wilx$diffexpressed <- "NO"
    wilx$diffexpressed[wilx$fc > 0 & wilx$wilx_pvalue > 1] <- "UP"
    wilx$diffexpressed[wilx$fc < 0 & wilx$wilx_pvalue > 1] <- "DOWN"

    # Print the wilx dataframe
    print(wilx)

    # Score based on the number of significant "UP" pathways
    score = sum(wilx$diffexpressed == "UP")
    return(score)
}


# File to store outputs
output_file <- "/home/haseena/volcano_output/wilx_output_file.txt"  # Change this to your desired output file path

# Redirecting output to a file
sink(output_file)

# Fine-tuning loop to explore different parameter combinations
best_score <- -Inf
best_parameters <- list()
best_result <- list()

for (quantile_thresh in quantile_range) {
    for (count_thresh in higher_count_cell_id_range) {
        for (anti_cor_value in anti_cor_value_range) {
            for (enrich_threshold in enrich_threshold_range) {
                for (mutation_quantile in mutation_quantile_range) {
                    
                    # Message for the current parameter combination
                    cat("Testing parameters: quantile_thresh = ", quantile_thresh, 
                        ", count_thresh = ", count_thresh,
                        ", anti_cor_value = ", anti_cor_value, 
                        ", enrich_threshold = ", enrich_threshold,
                        ", mutation_quantile = ", mutation_quantile, "\n")
                    
                    # Run the analysis for each parameter combination
                    result_g50 <- select_positives_negatives_for_ml(selected_g50, "g50", 0.08, 0.01,
                                                                      quantile_thresh, count_thresh, anti_cor_value,
                                                                      enrich_threshold, mutation_quantile)
                    result_g51 <- select_positives_negatives_for_ml(selected_g51, "g51", 0.08, 0.01,
                                                                      quantile_thresh, count_thresh, anti_cor_value,
                                                                      enrich_threshold, mutation_quantile)
                    result_g52 <- select_positives_negatives_for_ml(selected_g52, "g52", 0.08, 0.01,
                                                                      quantile_thresh, count_thresh, anti_cor_value,
                                                                      enrich_threshold, mutation_quantile)
                    result_g53 <- select_positives_negatives_for_ml(selected_g53, "g53", 0.08, 0.01,
                                                                      quantile_thresh, count_thresh, anti_cor_value,
                                                                      enrich_threshold, mutation_quantile)
                    
                    # Combine results to evaluate scoring
                    ctc_combined <- cbind(result_g50$ctc_score, result_g51$ctc_score, result_g52$ctc_score, result_g53$ctc_score)
                    non_ctc_combined <- cbind(result_g50$non_ctc_score, result_g51$non_ctc_score, result_g52$non_ctc_score, result_g53$non_ctc_score)

                    # Score the current parameter combination
                    score <- score_parameters(ctc_combined, non_ctc_combined, 0.01)
                    
                    # Message for the current score
                    cat("Score: ", score, "\n\n")
                    
                    # Check if this is the best score so far
                    if (score > best_score) {
                        best_score <- score
                        best_parameters <- list(quantile_thresh = quantile_thresh,
                                                count_thresh = count_thresh,
                                                anti_cor_value = anti_cor_value,
                                                enrich_threshold = enrich_threshold,
                                                mutation_quantile = mutation_quantile)
                        best_result <- list(ctc_score = ctc_combined, non_ctc_score = non_ctc_combined)
                    }
                }
            }
        }
    }
}

# Stop redirecting output
sink()

# Print the best parameters to the console
cat("Best parameters found:\n")
print(best_parameters)

# Plot the best parameter combination
plot_volcano(best_result$ctc_score, best_result$non_ctc_score, "best_combination", 0.01)
