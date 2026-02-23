# ------------------------------------------------------ #
# ------------------- Testing Function ----------------- #
# ------------------------------------------------------ #

source("helpers.R")
library(tidyverse)

A_star_search <- function(data,
                   A0 = character(0),
                   method = c("kw", "permute"),
                   alpha = 0.05,
                   n_perm = 1000) {
  
  method <- match.arg(method)
  A <- A0
  reps <- 0
  A_history <- list()
  
  if (method == "kw") {
    taxa <- rownames(data[[1]])
    uvw <- compute_uvw(data)
  }
  
  if (method == "permute") {
    corr_matrices <- lapply(data, function(mat) cor(t(mat)))
    taxa <- rownames(corr_matrices[[1]])
    permed_corr_matrices <- permute_data(data, n_perm)
  }
  
  repeat {
    
    # -----------------------------
    # METHOD-SPECIFIC TEST STAT
    # -----------------------------
    
    if (method == "kw") {
      p_values <- compute_kw_pvalues(data, uvw, A, taxa)
    }
    
    if (method == "permute") {
      p_values <- compute_permuted_pvalues(
        corr_matrices,
        permed_corr_matrices,
        A,
        taxa,
        n_perm
      )
    }
    
    # -----------------------------
    # B-H FDR correction
    # -----------------------------
    
    p.adj <- p.adjust(p_values, method = "BH")
    sig <- names(p_values)[!is.na(p.adj) & p.adj < alpha]
    
    if (length(A) == 1) {
      sig <- union(sig, A)
    }
    
    if (setequal(A, sig)) break
    
    # -----------------------------
    # Handling Detection Cases
    # -----------------------------
    
    sig_str <- paste(sort(sig), collapse = ",")
    history_strings <- sapply(A_history, paste, collapse = ",")
    
    if (sig_str %in% history_strings) {
      message("Exiting early: detected a repeating pattern in A.")
      break
    } else {
      A_history[[length(A_history) + 1]] <- sort(sig)
    }
    
    print(A)
    A <- sig
    reps <- reps + 1
    print(reps)
  }
  
  return(A)
}


# -------------------------------------------------------- #
# ----------------- K-W test computation ----------------- #
# -------------------------------------------------------- #

compute_kw_pvalues <- function(data, uvw, A, taxa) {
  
  n1 <- length(data[[1]])
  n2 <- length(data[[2]])
  n3 <- length(data[[3]])
  
  pvals <- numeric(length(taxa))
  names(pvals) <- taxa
  
  for (j in taxa) {
    
    M <- compute_m_vecs(uvw, A, j)
    
    m1 <- M$m1
    m2 <- M$m2
    m3 <- M$m3
    
    u_tild <- uvw[[1]][j,] * m1
    v_tild <- uvw[[2]][j,] * m2
    w_tild <- uvw[[3]][j,] * m3
    
    x <- c(u_tild, v_tild, w_tild)
    g <- c(rep(1, n1), rep(2, n2), rep(3, n3))
    
    pvals[j] <- kruskal.test(x = x, g = g)$p.value
  }
  
  return(pvals)
}

# ----------------------------------------------------- #
# --------------- p-vals from permuation -------------- #
# ----------------------------------------------------- #

compute_permuted_pvalues <- function(corr_matrices,
                                   permed_corr_mats_list,
                                   A,
                                   taxa,
                                   n_perm) {
  
  delta_values <- numeric(length(taxa))
  names(delta_values) <- taxa
  
  for (tax in taxa) {
    if (tax %in% A) {
      delta_values[tax] <- compute_delta(corr_matrices, tax, setdiff(A, tax))
    } else {
      delta_values[tax] <- compute_delta(corr_matrices, tax, A)
    }
  }
  
  delta_perm <- matrix(NA, nrow = length(taxa), ncol = n_perm)
  rownames(delta_perm) <- taxa
  
  for (perm in seq_len(n_perm)) {
    delta_perm[, perm] <- sapply(taxa, function(j) {
      compute_delta(permed_corr_mats_list[[perm]], j, A)
    })
  }
  
  p_values <- sapply(seq_along(delta_values), function(i) {
    mean(delta_perm[i, ] >= delta_values[i])
  })
  
  names(p_values) <- taxa
  return(p_values)
}



