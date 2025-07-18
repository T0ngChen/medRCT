# simulate data
K <- 4
n <- 100
data <- data.table::data.table(
  X = as.factor(sample(0:1, n, replace = TRUE)),
  M1 = rnorm(n),
  M2 = rbinom(n, 1, 0.5),
  M3 = rnorm(n),
  M4 = rnorm(n),
  C = rbinom(n, 1, 0.5)
)
dat2 <- data.table::copy(data)
k <- 1
first = 2
lnzero = 1
interactions_XC <- "X:C"
exposure_level <- c(0, 1)
fam_type = family_type(data, c("M1", "M2", "M3", "M4"))


test_that("test joint_dist", {
  result <- joint_dist(
    k = k,
    K = K,
    data = data,
    dat2 = dat2,
    fam_type = fam_type,
    interactions_XC = interactions_XC,
    exposure_level = exposure_level,
    separation_method = "discard",
    n = n
  )

  expect_true(all(c("M1", "M2", "M3") %in% colnames(result)))
  expect_true(all(
    paste0(
      "m",
      k,
      "_",
      exposure_level,
      "_",
      paste0(rep("m", K), collapse = "")
    ) %in%
      colnames(result)
  ))
  expect_equal(as.numeric(as.character(unique(result$X))), 1)
})


test_that("test marg_dist", {
  result <- marg_dist(
    k = k,
    first = 1,
    K = K,
    data = data,
    dat2 = dat2,
    fam_type = fam_type,
    interactions_XC = interactions_XC,
    separation_method = "discard",
    n = n
  )

  expect_true(paste0("m", k, "_0_", strrep("m", K)) %in% colnames(result))
  expect_equal(as.numeric(as.character(unique(result$X))), 0)
})


test_that("test joint_X_nonzero", {
  k = 3
  MM = 2
  index = setdiff(first:K, MM)
  result <- joint_X_nonzero(
    MM = MM,
    k = k,
    first = first,
    K = K,
    data = data,
    dat2 = dat2,
    fam_type = fam_type,
    interactions_XC = interactions_XC,
    lnzero = lnzero,
    n = n,
    index = index,
    separation_method = "discard"
  )

  for (a in lnzero) {
    var_name <- paste0(
      "m",
      k,
      "_",
      a,
      "_",
      paste0(
        c(
          rep(paste0(a), min(k - 1, MM - 1)),
          "m",
          rep(paste0(a), max(k - 1 - MM, 0)),
          rep("m", K - 1 - min(k - 1, MM - 1) - max(k - 1 - MM, 0))
        ),
        collapse = ""
      )
    )
    expect_true(var_name %in% colnames(result))
  }
  expect_equal(length(result[[var_name]]), n)

  # Check that dat2 is returned early when MM != 1 and k == index[1]
  MM <- 2
  k <- index[1]
  result_early <- joint_X_nonzero(
    MM = MM,
    k = k,
    first = first,
    K = K,
    data = data,
    dat2 = dat2,
    fam_type = fam_type,
    interactions_XC = interactions_XC,
    lnzero = lnzero,
    n = n,
    index = index
  )
  expect_equal(result_early, dat2)

  expect_equal(as.numeric(as.character(unique(result$X))), 1)
})


test_that("test con_exposed", {
  MM = 1
  k = 2
  result <- con_exposed(
    MM = MM,
    k = k,
    K = K,
    data = data,
    dat2 = dat2,
    fam_type = fam_type,
    interactions_XC = interactions_XC,
    lnzero = lnzero,
    n = n,
    separation_method = "brglm"
  )

  a = lnzero
  var_name <- paste0(
    "m",
    k,
    "_",
    a,
    "_",
    paste0(
      c(
        rep(paste0(a), MM - 1),
        0,
        rep(paste0(a), max(k - 1 - MM, 0)),
        rep("m", K - MM - max(k - 1 - MM, 0))
      ),
      collapse = ""
    )
  )
  expect_true(var_name %in% colnames(result))
  expect_equal(as.numeric(as.character(unique(result$X))), 1)
})


test_that("test joint_unexposed", {
  first = 2
  for (k in first:K) {
    dat2 <- marg_dist(
      k = k,
      first = first,
      K = K,
      data = data,
      dat2 = dat2,
      fam_type = fam_type,
      interactions_XC = interactions_XC,
      separation_method = "discard",
      n = n
    )
  }
  k = 3
  result <- joint_unexposed(
    k = k,
    first = first,
    K = K,
    data = data,
    dat2 = dat2,
    fam_type = fam_type,
    interactions_XC = interactions_XC,
    n = n
  )

  l <- first:(k - 1)
  var_name <- med_outcome_all(l = k, first = first, a = 0, K = K)
  expect_true(var_name %in% colnames(result))
  expect_equal(unique(as.numeric(as.character(result$X))), 0)
})
