test_that("pwf_stack retains directed edges", {
  W1 <- matrix(c(0, 2, -3,
                 -2, 0, 4,
                  3, -4, 0), nrow = 3, byrow = TRUE)
  W2 <- -W1
  stk <- pwf_stack(list(subjA = W1, subjB = W2), labels = c("A", "B", "C"),
                   triangle = "both")
  expect_s3_class(stk, "mppi_stack")
  expect_true(stk$directed)
  expect_equal(ncol(stk$matrix), 6)
  expect_true(all(grepl("->", colnames(stk$matrix))))
  mask_vals <- as.numeric(W1[stk$mask])
  expect_equal(unname(stk$matrix[1, ]), unname(mask_vals))
})

test_that("mppi_group_glm handles directed stacks", {
  W1 <- matrix(c(0, 1, -2,
                 -1, 0, 3,
                  2, -3, 0), nrow = 3, byrow = TRUE)
  W2 <- matrix(c(0, 2, -1,
                 -2, 0, 4,
                  1, -4, 0), nrow = 3, byrow = TRUE)
  stk <- pwf_stack(list(W1, W2), labels = c("L", "M", "R"), triangle = "upper")
  design <- matrix(1, nrow = 2, ncol = 1)
  colnames(design) <- "Intercept"
  res <- mppi_group_glm(stk, design, method = "lm")
  expect_true(all(res$directed))
  expect_equal(res$triangle[1], "upper")
  # manual check for first edge (L->M corresponds to entry 1->2)
  edge_idx <- which(res$edge == "L->M")
  expect_equal(res$estimate[edge_idx], mean(c(W1[1,2], W2[1,2])))
})

test_that("mppi_group_ebayes supports directed option", {
  set.seed(5)
  mats <- replicate(4, {
    M <- matrix(0, 3, 3)
    M[upper.tri(M)] <- rnorm(3)
    M <- M - t(M)
    M
  }, simplify = FALSE)
  out_dir <- mppi_group_ebayes(mats, use_limma = FALSE, directed = TRUE,
                               triangle = "upper")
  expect_length(out_dir$t, 3)
  out_both <- mppi_group_ebayes(mats, use_limma = FALSE, directed = TRUE,
                                triangle = "both")
  expect_length(out_both$t, 6)
})
