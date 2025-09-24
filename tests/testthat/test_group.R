context("Group-level utilities")

test_that("mppi_group_ebayes detects signal", {
  set.seed(200)
  r <- 4L
  n_subj <- 6L
  # create small symmetric matrices with a common signal on edge (1,2)
  M_list <- vector("list", n_subj)
  for (s in seq_len(n_subj)) {
    base <- matrix(rnorm(r * r, sd = 0.1), r, r)
    base <- (base + t(base)) / 2
    base[1,2] <- base[2,1] <- base[1,2] + 0.5
    M_list[[s]] <- base
  }
  res <- mppi_group_ebayes(M_list, use_limma = FALSE)
  idx_mat <- matrix(FALSE, res$V, res$V)
  idx_mat[upper.tri(idx_mat)] <- res$p < 0.1
  expect_true(idx_mat[1,2])
})
