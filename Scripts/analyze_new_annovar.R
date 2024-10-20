

setwd("/data1/haseena_IP/tsv_files/sandisk/sandisk_output3/BAM/annovar/")

file_list <- list.files(pattern = "*.hg38_multianno.txt") # path to the indivdual annovar output 
genes_selected <- read.csv("/data1/haseena_IP/tsv_files/sandisk/sandisk_output3/BAM/all_work/genes_to_select_pancreas_cancer.csv") #  genes which you taken


data_merge_pat1 <-  data.table::fread("/data1/haseena_IP/tsv_files/sandisk/sandisk_output3/BAM/filtered_bam/vcf/annovar/final_spatial_converted.hg38_multianno.txt", header = TRUE, sep = "\t") # merged annovar output normal

#data_merge_pat1 <- data_merge_pat1[which(data_merge_pat1$Gene.refGene %in% genes_selected[,1]),]

mat_data_pat <- matrix(NA, 0, 10)
colnames(mat_data_pat) <- c("Chr", "Start", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", "AAChange.refGene", "Otherinfo11", "cell_id")

for (i in 1:length(file_list)) {
  file_path <- paste0("/data1/haseena_IP/tsv_files/sandisk/sandisk_output3/BAM/annovar/", file_list[i])
  
  if (file.info(file_path)$size == 0) { # Check if the file is empty
    next
  }
  
  data_pre <- read.table(file_path, header = TRUE, sep = "\t")
  data_pre1 <- data_pre[which(data_pre$Gene.refGene %in% genes_selected$x),]
  if (nrow(data_pre1) == 0) {next}
  
  data <- data_pre1[, c("Chr", "Start", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", "AAChange.refGene", "Otherinfo11")] 
  name_pre <- gsub("lofreq_", "", file_list[i])
  data$cell_id <- gsub("_.hg38_multianno.txt", "", name_pre)
  mat_data_pat <- rbind(mat_data_pat, data)
}

mat_data_pat$Chr_Start <- paste0(mat_data_pat$Chr, "_", mat_data_pat$Start)



build_mut_freq_per_cell <- function(mat_data) {

new_mat <- matrix(NA, length(unique(paste(mat_data$Chr, mat_data$Start, sep="_"))),  length(unique(mat_data$cell_id))+2)
new_mat[,1] <- unique(paste(mat_data$Chr, mat_data$Start, sep="_"))
colnames(new_mat) <- c("Chr_Start", "number_of_cells",  unique(mat_data$cell_id))

for (i in 1:nrow(new_mat)) {
    pos_list <- which(mat_data$Chr_Start %in% new_mat[i,1])
    for (pos in pos_list) {
          cell_info = paste(mat_data$Ref, mat_data$Alt, mat_data$Gene.refGene, mat_data$ExonicFunc.refGene, mat_data$RS.ID, sep = ",")[pos]
    cell_col_name <- mat_data$cell_id[pos]
    cell_info <- as.data.frame(cell_info); colnames(cell_info) <- cell_col_name 
    new_mat[i, mat_data$cell_id[pos]] <-  as.character(cell_info[,1])
      }

    new_mat[i, "number_of_cells"] = length(pos_list)

}
new_mat <- as.data.frame(new_mat)
}


check_if_matching_alterations <- function (profile_af_pre, mat_data) {

af_value = as.numeric(profile_af_pre$AF) * length(unique(mat_data$cell_id))

if (length(which(af_value > 10)) != 0) {profile_af <- profile_af_pre[-which(af_value > 10),]} else {
  profile_af <- profile_af_pre}

profile_af_present <- as.data.frame(matrix(NA, 0, ncol(profile_af))) 

for (i in 1:nrow(mat_data)) {

pos_pre <- which(mat_data$Start[i] %in% profile_af$Start  )

if (length(pos_pre) == 0) {#message("locus not present");
 next}

profile_af_i <- profile_af[which(profile_af$Start %in% mat_data$Start[i]),]
# message("dimension ", nrow(profile_af_i))

if (paste0(mat_data$Ref[i], "-", mat_data$Alt[i]) %in% paste0(profile_af_i$Ref, "-", profile_af_i$Alt) ) {
  pos_alt <- which(paste0(profile_af_i$Ref, "-", profile_af_i$Alt) %in% paste0(mat_data$Ref[i], "-", mat_data$Alt[i]))
  # message(paste0(mat_data$Ref[i], "-", mat_data$Alt[i]), "   ", paste0(profile_af_i$Ref, "-", profile_af_i$Alt))
  # print(i)
  profile_af_present <- rbind(profile_af_present, profile_af_i[pos_alt,])}

}
return(profile_af_present)
}


select_ctcs_by_mutations <- function (mat_data, selected_mutations) {
mat_data$considered_somatic_mutation <- 0

mat_data$considered_somatic_mutation[which(mat_data$Start %in% selected_mutations)] <- 1
selected_ctcs <- mat_data[which(mat_data$considered_somatic_mutation == 1),]
no_som_mut_per_cell <- as.data.frame(table(selected_ctcs$cell_id))
no_som_mut_per_cell <- no_som_mut_per_cell[order(no_som_mut_per_cell$Freq, decreasing = T), ]
colnames(no_som_mut_per_cell) <- c("cell_id", "no_somatic_mutations")


cells_alt_mut_info <- as.data.frame(matrix(NA, length(unique(selected_ctcs$cell_id)), 7  * no_som_mut_per_cell$no_somatic_mutations[1] ))
col_names_mut <- make_colnames_for_cell_alt_mut_1(rep("something", no_som_mut_per_cell$no_somatic_mutations[1]))
colnames(cells_alt_mut_info) <- col_names_mut
rownames(cells_alt_mut_info) <- no_som_mut_per_cell$cell_id

for ( i in 1:length(unique(selected_ctcs$cell_id))) {

pos_cells <- which(selected_ctcs$cell_id %in% unique(selected_ctcs$cell_id)[i])

  if (length(pos_cells) == 1 ) { 
    row_to_append <- selected_ctcs[pos_cells, -which(colnames(selected_ctcs[pos_cells,]) %in% c("End", "considered_somatic_mutation", "cell_id", "Chr_Start", "GeneDetail.refGene", "AAChange.refGene" ))]
    rownames(row_to_append) <- NULL
    col_names_ii <- make_colnames_for_cell_alt_mut_1(1)
    rownames(row_to_append) <- unique(selected_ctcs$cell_id)[i]
    colnames(row_to_append) <- col_names_ii

    cells_alt_mut_info[unique(selected_ctcs$cell_id)[i], 1:length(row_to_append)] <-  row_to_append } 


  if (length(pos_cells) > 1 )  {

    row_i_df <-  data.frame(NA, 1, 0)

  for (row in pos_cells) { 
    row_to_add_pre <- selected_ctcs[row, ]
    row_to_add <- row_to_add_pre[ , -which(colnames(row_to_add_pre) %in% c("End", "considered_somatic_mutation", "cell_id", "Chr_Start", "GeneDetail.refGene", "AAChange.refGene" ))]
    rownames(row_to_add) <- NULL
    colnames(row_to_add) <- NULL
    rownames(row_to_add) <- unique(selected_ctcs$cell_id)[i]
    row_i_df <- cbind(row_i_df, row_to_add)
  }

    row_i_df_final <- row_i_df[,-which(colnames(row_i_df) %in% c("NA.", "X1", "X0"))]

    col_names_list <- make_colnames_for_cell_alt_mut_1(pos_cells)
    colnames(row_i_df_final) <- col_names_list 
    cells_alt_mut_info[unique(selected_ctcs$cell_id)[i], 1:length(row_i_df_final)] <-  row_i_df_final
}

}

no_som_mut_per_cell_final <- cbind(no_som_mut_per_cell, cells_alt_mut_info)
no_som_mut_per_cell_final[is.na(no_som_mut_per_cell_final)] <- ""
return(no_som_mut_per_cell_final)

}

make_colnames_for_cell_alt_mut_1 <- function(pos_cells) {
    col_names_list <- list()
    for (p in 1:length(pos_cells)) { ## this  is for making the columnames for the row_i_df_final
      col_names <-  c("Chr", "Start", "Ref", "Alt", "Func.refGene", "Gene.refGene",  "Otherinfo11")
      col_names_suffix <- unlist(lapply(col_names, paste0, paste0("_mut", p)))
      col_names_list <- append(col_names_list, col_names_suffix)
    }
    col_names_list <- unlist(col_names_list)
  return(col_names_list)
}


overlap_lofreq_merge_vs_lofreq_ctc <- function(lofreq_merge_af, mut_data, AF_threshold, cell_freq_threshold) {

profile_af_present <- check_if_matching_alterations(lofreq_merge_af, mut_data)
profile_af_present_f1 <- profile_af_present[which(profile_af_present$AF < AF_threshold),] # Rule 1: Must have AF of less than 60%
mut_per_cell <- build_mut_freq_per_cell(mut_data)

skipped_locus <- list()
selected_mutations <- list()


for (mut in 1:nrow(profile_af_present_f1)) { 
    locus <- profile_af_present_f1[mut, "Start"] 
    # print(locus)
    locus_row <- as.numeric(mut_per_cell[which(gsub(".*_", "", mut_per_cell$Chr_Start) %in%  factor(locus)), "number_of_cells"])
    # message("row_id: ", mut)
    # message("locus_row: ", locus_row)
    if (length(locus_row) == 0 ) { # Rule 1: Must be less than threshold number of cells
      skipped_locus <- append(skipped_locus, locus); next} 

    if (locus_row <= cell_freq_threshold * (ncol(mut_per_cell) - 2 ))  {
      # message("selected locus: ", locus)
      selected_mutations <- append(selected_mutations, locus)
        # message(locus_row)
    } else {  skipped_locus <- append(skipped_locus, locus) } 

}

    if (mut == nrow(profile_af_present_f1)) {
        message(length(skipped_locus), " locus skipped on quntile filter") 
        message(length(unique(selected_mutations)), " loci selected out of ", nrow(mut_per_cell), " loci") 
        }

selected_cells_pre  <- select_ctcs_by_mutations(mut_data, selected_mutations)
cell_ids_not_ranked <- unique(mut_data$cell_id)[-which(unique(mut_data$cell_id) %in% selected_cells_pre$cell_id)]
cell_ids_not_ranked <- as.data.frame(cell_ids_not_ranked)
colnames(cell_ids_not_ranked) <- "cell_id"

selected_cells <- plyr::rbind.fill(selected_cells_pre, cell_ids_not_ranked)
selected_cells[which(as.character(selected_cells$cell_id) %in% cell_ids_not_ranked$cell_id), "no_somatic_mutations"] <- 0
selected_cells[is.na(selected_cells)] <- ""
message("No. of cells ranked: ", nrow(selected_cells_pre), " out of ", length(unique(mut_data$cell_id)), " cells" )
print(selected_cells_pre$no_somatic_mutations)
return(selected_cells)


}



















split_new_pat1_all <- splitstackshape::cSplit(data_merge_pat1, "Otherinfo11", sep = ";", direction = "wide")
lofreq_merge_af_pat1_all <- cbind(data_merge_pat1, gsub("AF=", "", split_new_pat1_all$Otherinfo11_01))
colnames(lofreq_merge_af_pat1_all) <- c(colnames(data_merge_pat1), "AF")
cell_rank_pat1_all <-  overlap_lofreq_merge_vs_lofreq_ctc(lofreq_merge_af_pat1_all, mat_data_pat, 0.6, 0.2)






file_path <- "spatial_1.csv"

# Write the data to a CSV file
write.csv(cell_rank_pat1_all, file = file_path, row.names = FALSE)

