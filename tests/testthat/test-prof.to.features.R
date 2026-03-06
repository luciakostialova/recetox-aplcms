patrick::with_parameters_test_that(
  "test prof.to.features_bigauss",
  {
    testdata <- file.path("..", "testdata")
    input_path <- file.path(testdata, "filtered", filename)
    extracted_features <- arrow::read_parquet(input_path)

    actual <- prof.to.features(
      profile = extracted_features,
      bandwidth = 0.25,
      min_bandwidth = 0.25,
      max_bandwidth = 1,
      sd_cut = sd_cut,
      sigma_ratio_lim = sigma_ratio_lim,
      shape_model = shape_model,
      peak_estim_method = "moment",
      moment_power = 1,
      component_eliminate = 0.01,
      BIC_factor = 2,
      do.plot = do.plot
    )

    expected_path <- file.path(testdata, "features", expected_filename)
    expected <- arrow::read_parquet(expected_path)



    expect_equal(dplyr::arrange(actual, dplyr::pick(area)), dplyr::arrange(expected, dplyr::pick(area)), tolerance = 0.01) # sort by area to avoid switched values
  },
  patrick::cases(
    mbr_test0 = list(
      filename = c("mbr_test0.parquet"),
      expected_filename = "mbr_test0_features.parquet",
      sd_cut = c(0.1, 100),
      sigma_ratio_lim = c(0.1, 10),
      shape_model = "bi-Gaussian",
      do.plot = FALSE
    )
    #,
    # thermo_raw_profile = list(
    #   filename = c("thermo_raw_profile.parquet"),
    #   expected_filename = "thermo_raw_profile_features.parquet",
    #   sd_cut = c(0.1, 1),
    #   sigma_ratio_lim = NA,
    #   shape_model = "Gaussian",
    #   do.plot = FALSE
    # ),
    # RCX_06_shortened_gaussian = list(
    #   filename = c("RCX_06_shortened.parquet"),
    #   expected_filename = "RCX_06_shortened_gaussian_features.parquet",
    #   sd_cut = c(0.01, 500),
    #   sigma_ratio_lim = c(0.01, 100),
    #   shape_model = "Gaussian",
    #   do.plot = FALSE
    # ),
    # RCX_06_shortened_v2 = list(
    #   filename = c("RCX_06_shortened.parquet"),
    #   expected_filename = "RCX_06_shortened_features.parquet",
    #   sd_cut = c(0.01, 500),
    #   min_bandwidth = 0.5,
    #   max_bandwidth = 1.5,
    #   sigma_ratio_lim = c(0.01, 100),
    #   shape_model = "bi-Gaussian",
    #   do.plot = FALSE
    # )
   )
  )


test_that("Unit testing for normix function", {

  # Test 1: small Gaussian mixture
  x <- seq(0, 6, length.out = 61)
  y <- 50 * dnorm(x, 2, 0.3) + 30 * dnorm(x, 4, 0.4)
  rt_profile <- data.frame('rt' = x, 'intensity' = y)

  pks <- c(2, 4)
  vlys <- c(-Inf, 3, Inf)
  
  actual <- normix(rt_profile, pks, vlys, ignore=0, max.iter=50, aver_diff=mean(diff(x)))
  expected <- cbind(
    miu = c(2, 4),
    sigma = c(0.3, 0.4),
    scale = c(50, 30)    
  )
  expect_equal(actual, expected, tolerance = 1)

  # Test 2: large Gaussian mixture
  n <- 2000
  x <- seq(0, 190, length.out = n)
  peaks_def <- list(
    list(mu=10, sigma=3, scale=30),
    list(mu=34, sigma=2.5, scale=40),
    list(mu=50, sigma=2.5, scale=45),
    list(mu=77, sigma=2.2, scale=39),
    list(mu=101, sigma=3, scale=25),
    list(mu=124, sigma=2.5, scale=33),
    list(mu=144, sigma=2.7, scale=41),
    list(mu=169, sigma=3.3, scale=54)
  )
  y <- rep(0, n)
  for (peak in peaks_def) {
    y <- y + peak$scale * dnorm(x, peak$mu, peak$sigma)
  }
  rt_profile <- data.frame('rt' = x, 'intensity' = y)
  pks <- sapply(peaks_def, function(peak) peak$mu)
  vlys <- c(-Inf, 22, 42, 64, 87, 113, 134, 155, Inf)
  actual <- normix(rt_profile, pks, vlys, ignore=0, max.iter=50, aver_diff=mean(diff(x)))
  expected <- cbind(
    miu = sapply(peaks_def, function(peak) peak$mu),
    sigma = sapply(peaks_def, function(peak) peak$sigma),
    scale = sapply(peaks_def, function(peak) peak$scale)
  )
  expect_equal(actual, expected, tolerance = 1) 
})

test_that("Unit testing for normix function, overlapping peaks", {
  # Overlaps
  n <- 1000
  x <- seq(0, 100, length.out = n)

  peaks_def <- list(
    list(mu = 25, sigma = 4, scale = 40),
    list(mu = 36, sigma = 5, scale = 45), 
    list(mu = 60, sigma = 4, scale = 40),
    list(mu = 75, sigma = 7, scale = 40)   
  )

  y <- rep(0, n)
  for (peak in peaks_def) {
    y <- y + peak$scale * dnorm(x, peak$mu, peak$sigma)
  }

  rt_profile <- data.frame('rt' = x, 'intensity' = y)
  turns <- find.turn.point(y)
  pks <- rt_profile[, 1][turns$pks]
  vlys <- c(-Inf, rt_profile[, 1][turns$vlys], Inf)
  
  actual <- normix(rt_profile, pks, vlys, ignore=0, max.iter=50, aver_diff=mean(diff(x)))
  expected <- cbind(
    miu = sapply(peaks_def, function(peak) peak$mu),
    sigma = sapply(peaks_def, function(peak) peak$sigma),
    scale = sapply(peaks_def, function(peak) peak$scale)
  )

  expect_equal(actual, expected, tolerance = 1)  # works, scale estimates are a bit off
   
})