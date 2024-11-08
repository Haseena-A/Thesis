# Read the gene expression matrix and transpose it
expression_matrix <- t(read.csv("/home/haseena/0h/all_0h_stochastic.csv", row.names=1, header=TRUE, stringsAsFactors=FALSE))



# Define the binorm function
binorm <- function(x) {
  xdim = dim(x)
  xpval = matrix(1, xdim[1], xdim[2], dimnames=list(toupper(rownames(x)), colnames(x)))
  for (i in 1:xdim[2]) {
    pos = which(abs(x[, i]) > 0.00000000001)
    pnz = length(pos) / xdim[1]
    meannz = 0
    sdnz = sd(x[pos, i])
    xpval[pos, i] = (1 - pnorm(x[pos, i], meannz, sdnz))
  }
  return(xpval)
}

# Define population variance and standard deviation functions
pop.var <- function(x) var(x) * (length(x) - 1) / length(x)
pop.sd <- function(x) sqrt(pop.var(x))

# Define the transformData function
transformData <- function(data_vector) {
  dvm = mean(data_vector)
  dvsd = pop.sd(data_vector)
  s = (data_vector - dvm) / dvsd
  distr = ecdf(s)
  sapply(s, function(a) -2 * log(distr(a)))
}

# Define the calculateCovariances function
calculateCovariances <- function(data_matrix) {
  transformed_data_matrix = apply(data_matrix, MARGIN=1, FUN=transformData)
  covar_matrix = cov(transformed_data_matrix)
  covar_matrix
}

# Define the combinePValues function
combinePValues <- function(covar_matrix, p_values, extra_info = FALSE) {
  N = ncol(covar_matrix) # number of samples
  df_fisher = 2.0 * N
  Expected  = 2.0 * N
  cov_sum <- (2 * sum(covar_matrix[lower.tri(covar_matrix, diag=FALSE)]))
  Var = 4.0 * N + cov_sum
  c = Var / (2.0 * Expected)
  df_brown = (2.0 * Expected^2) / Var
  if (df_brown > df_fisher) {
    df_brown = df_fisher
    c = 1.0
  }
  x = 2.0 * sum(-log(p_values))

  p_brown = pchisq(df = df_brown, q = x / c, lower.tail = FALSE)
  p_fisher = pchisq(df = df_fisher, q = x, lower.tail = FALSE)

  if (extra_info) {
    return(list(P_test = p_brown, P_Fisher = p_fisher, Scale_Factor_C = c, DF = df_brown))
  } else {
    return(p_brown)
  }
}

# Define the combine function
combine <- function(gene_file, expression_matrix, gnames, Pval1, thr = 2) {
  combined_pvals <- matrix(0, nrow(gene_file), ncol(expression_matrix), dimnames=list(rownames(gene_file), colnames(expression_matrix)))
  gene_file = data.frame(apply(gene_file, 2, toupper))
  rownames(expression_matrix) = toupper(rownames(expression_matrix))
  gnames = toupper(gnames)

  for (i in 1:nrow(gene_file)) {
    pathway <- gene_file[i, ]
    genes = pathway[2:ncol(pathway)]

    pose = which(genes == "")
    if (dim(as.matrix(pose))[1] > 0) {
      genes = genes[-pose]
    }
    genes = as.matrix(genes)
    pos = which(genes %in% gnames)
    pos = as.matrix(pos)

    if (nrow(pos) > 5) {
      exprlocal = abs(expression_matrix[genes[pos], ])
      plocal <- Pval1[genes[pos], ]
      plocal = t(plocal)
      covar_matrix = cov(plocal)

      for (j in 1:nrow(plocal)) {
        gpos = which(exprlocal[, j] > 0)
        if (nrow(as.matrix(gpos)) > thr) {
          lcovar = covar_matrix[gpos, gpos]
          temp = combinePValues(lcovar, p_values = plocal[j, gpos], extra_info = TRUE)
          combined_pvals[i, j] = -log2(temp$P_test)
        }
      }
    }
  }
  return(combined_pvals)
}

# Example usage
gene_file <- read.csv("/home/haseena/hesc_lineage.csv", header = TRUE, stringsAsFactors = FALSE)
Pval1 <- binorm(expression_matrix)

# Combine P-values
combined_pvals <- combine(gene_file, expression_matrix, gnames = rownames(expression_matrix), Pval1, thr = 2)

# Write the combined P-values to a CSV file
write.csv(combined_pvals, "/home/haseena/combined_HOUR0_velocity_stochaastic.csv", row.names = TRUE)

