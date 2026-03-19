#' Validate that provided model and method inputs match expected, exit execution otherwise.
#' @param shape_model The mathematical model for the shape of a peak. There are two choices - "bi-Gaussian" and "Gaussian".
#'  When the peaks are asymmetric, the bi-Gaussian is better. The default is "bi-Gaussian".
#' @param peak_estim_method The estimation method for the bi-Gaussian peak model. Two possible values: moment and EM.
#' @export
validate_model_method_input <- function(shape_model, peak_estim_method) {
  if (!shape_model %in% c("Gaussian", "bi-Gaussian")) {
    stop("shape_model argument must be 'Gaussian' or 'bi-Gaussian'")
  }
  if (!peak_estim_method %in% c("moment", "EM")) {
    stop("peak_estim_method argument must be 'moment' or 'EM'")
  }
}


#' Initialize minimum and maximum bandwidth values if none given. Ensure that minimum bandwidth is lower that maximum, else set minimum to 1/4 of maximum value.
#' @param min_bandwidth The minimum bandwidth to use in the kernel smoother.
#' @param max_bandwidth The maximum bandwidth to use in the kernel smoother.
#' @param profile Profile table with shape number-of-features*4. The table contains following columns:
#' \itemize{
#'   \item mz - float - mass-to-charge ratio of feature
#'   \item rt - float - retention time of features
#'   \item intensity - float - intensity of features
#'   \item group_number - integer - group number assigned to each feature based on their rt similarity
#' }
#' @return Returns a list object with the following objects in it:
#' \itemize{
#'   \item min_bandwidth - float - Minimum bandwidth.
#'   \item max_bandwidth - float - Maximum bandwidth
#' }
#' @export
preprocess_bandwidth <- function(min_bandwidth, max_bandwidth, profile) {
  if (is.na(min_bandwidth)) {
    min_bandwidth <- diff(range(profile[, 'rt'], na.rm = TRUE)) / 60
  }
  if (is.na(max_bandwidth)) {
    max_bandwidth <- diff(range(profile[, 'rt'], na.rm = TRUE)) / 15
  }
  if (min_bandwidth >= max_bandwidth) {
    min_bandwidth <- max_bandwidth / 4
  }
  return(data.frame("min_bandwidth" = min_bandwidth, "max_bandwidth" = max_bandwidth))
}

#' Convert matrix to dataframe and rename columns.
#' @param profile Profile table with shape number-of-features*4. The table contains following columns:
#' \itemize{
#'   \item float - mass-to-charge ratio of feature
#'   \item float - retention time of features
#'   \item float - intensity of features
#'   \item integer - group number assigned to each feature based on their rt similarity
#' }
#' @return  Returns a dataframe with shape number-of-features*4. The columns are as follows:
#' \itemize{
#'   \item mz - float - mass-to-charge ratio of feature
#'   \item rt - float - retention time of features
#'   \item intensity - float - intensity of features
#'   \item group_number - integer - group number assigned to each feature based on their rt similarity
#' }
#' @export
preprocess_profile <- function(profile) {
  keys <- c("mz", "rt", "intensity", "group_number")
  colnames(profile) <- keys
  return(data.frame(profile))
}

#' Compute parameters of chromatographic peak shape if peaks are considered to be gaussian
#' @param rt_profile A matrix with two columns: "base.curve" (rt) and "intensity".
#' @param bw Bandwidth vector to use in the kernel smoother.
#' @param component_eliminate When a component accounts for a proportion of intensities less than this value, the component will be ignored.
#' @param BIC_factor The factor that is multiplied on the number of parameters to modify the BIC criterion. If larger than 1,
#'  models with more peaks are penalized more.
#' @param aver_diff Average retention time difference across RTs of all features.
#' @return Returns a single-row vector or a table object with the following items/columns:
#' \itemize{
#'   \item miu - float - mean value of the gaussian curve
#'   \item sigma - float - standard deviation of the gaussian curve
#'   \item sigma - float - standard deviation of the gaussian curve
#'   \item scale - float - estimated total signal strength (total area of the estimated normal curve)
#'}
#' @export
compute_gaussian_peak_shape <- function(rt_profile, bw, component_eliminate, BIC_factor, aver_diff) {
  rt_peak_shape <- normix.bic(rt_profile, bw = bw, eliminate = component_eliminate, BIC_factor = BIC_factor, aver_diff = aver_diff)$param # a matrix is returned, even if has one row
  if (any(is.na(rt_peak_shape))) return(NULL)
  rt_peak_shape <- rt_peak_shape[, c(1, 2, 2, 3)]  
  return(rt_peak_shape)
}

#' Compute 'a' parameter of bi-gaussian chromatographic peak shape - the m param
#' @description
#' This function solves peak summit using the x, t, peak summit from the previous step, 
#' and sigma.1, and sigma.2 (original authors' comment).
#' @param x A vector of numerical values (intensities).
#' @param t A vector of numerical values (rt).
#' @param a A vector of peak summits.
#' @param sigma.1 left standard deviation of the gaussian curve
#' @param sigma.2 right standard deviation of the gaussian curve
#' @export
solve_a <- function(x, t, a, sigma.1, sigma.2) {
  # This function is a part of bigauss.esti.EM and is not covered by any of test-cases
  w <- x * (as.numeric(t < a) / sigma.1 + as.numeric(t >= a) / sigma.2)
  return(sum(t * w) / sum(w))
}

#' Compute 'u' and 'v' parameters of bi-gaussian chromatographic peak shape
#' @description
#' This function prepares the parameters required for latter computation.
#' u, v, and sum of x (original authors' comment).
#' @param x A vector of numerical values (intensities).
#' @param t A vector of numerical values (rt).
#' @param a A vector of peak summits.
#' @export
prep_uv <- function(x, t, a) {
  # This function is a part of bigauss.esti.EM and is not covered by any of test-cases
  temp <- (t - a)^2 * x
  u <- sum(temp * as.numeric(t < a))
  v <- sum(temp * as.numeric(t >= a))
  return(list(
    u = u,
    v = v,
    x.sum = sum(x)
  ))
}

#' Solve sigma
#' @description
#' Calculates the square estimated of left and right standard deviations.
#' @param x A vector of numerical values (intensities).
#' @param t A vector of numerical values (rt).
#' @param a A vector of peak summits.
#' @return A vector of:
#' \itemize{
#'   \item standard deviation at the left side of the gaussian curve
#'   \item standard deviation at the right side of the gaussian curve
#' }
#' @export
solve_sigma <- function(x, t, a) {
  # This function is a part of bigauss.esti.EM and is not covered by any of test-cases
  tt <- prep_uv(x, t, a)
  sigma.1 <- tt$u / tt$x.sum * ((tt$v / tt$u)^(1 / 3) + 1)
  sigma.2 <- tt$v / tt$x.sum * ((tt$u / tt$v)^(1 / 3) + 1)
  return(list(
    sigma.1 = sigma.1,
    sigma.2 = sigma.2
  ))
}

#' Compute bi-gaussian parameters using Expectation-Maximization algorithm
#' @description
#' Computes bi-gaussian parameters using an iterative method - Expectation-Maximization algorithm.
#' The algorithm iteratively updates the parameters until convergence or until the maximum number of iterations is reached.
#' @param t A vector of numerical values (rt).
#' @param x A vector of numerical values (intensities).
#' @param max.iter Maximum number of iterations.
#' @param epsilon Threshold for continuing the iteration
#' @param sigma_ratio_lim A vector of two. It enforces the belief of the range of the ratio between the left-standard deviation 
#' and the right-standard deviation of the bi-Gaussian function.
#' @return A vector with length 4. The items are as follows going from first to last:
#' \itemize{
#'   \item mean of gaussian curve
#'   \item standard deviation at the left side of the gaussian curve
#'   \item standard deviation at the right side of the gaussian curve
#'   \item estimated total signal strength (total area of the estimated normal curve)
#' }
#' @export
bigauss.esti.EM <- function(t, x, max.iter = 50, epsilon = 0.005, do.plot = FALSE, sigma_ratio_lim = c(0.3, 1)) {
  # This function is not covered by any test case
  sel <- which(x > 1e-10)
  if (length(sel) == 0) {
    return(c(median(t), 1, 1, 0))
  }
  if (length(sel) == 1) {
    return(c(t[sel], 1, 1, 0))
  }
  t <- t[sel]
  x <- x[sel]

  ## epsilon is the threshold for continuing the iteration. change in
  ## a smaller than epsilon will terminate the iteration.
  ## epsilon <- min(diff(sort(t)))/2

  ## using the median value of t as the initial value of a (peak summit)
  a.old <- t[which(x == max(x))[1]]
  a.new <- a.old
  change <- 10 * epsilon

  ## n.iter is the number of iteration covered so far.
  n.iter <- 0

  while ((change > epsilon) & (n.iter < max.iter)) {
    a.old <- a.new
    n.iter <- n.iter + 1
    sigma <- solve_sigma(x, t, a.old)
    if (n.iter == 1) {
        sigma[is.na(sigma)] <- as.numeric(sigma[which(!is.na(sigma))])[1] / 10
    }
    a.new <- solve_a(x, t, a.old, sigma$sigma.1, sigma$sigma.2)
    change <- abs(a.old - a.new)
  }
  d <- x
  sigma$sigma.2 <- sqrt(sigma$sigma.2)
  sigma$sigma.1 <- sqrt(sigma$sigma.1)

  d[t < a.new] <- dnorm(t[t < a.new], mean = a.new, sd = sigma$sigma.1) * sigma$sigma.1
  d[t >= a.new] <- dnorm(t[t >= a.new], mean = a.new, sd = sigma$sigma.2) * sigma$sigma.2
  scale <- exp(sum(d[d > 1e-3]^2 * log(x[d > 1e-3] / d[d > 1e-3])) / sum(d[d > 1e-3]^2))
  return(c(a.new, sigma$sigma.1, sigma$sigma.2, scale))
}

#' Reverse cumulative sum
#' @description
#' Computes vector of cumulative sums on reversed input. Returns cumulative sum vector going from the sum of all elements to one.
#' @param x float - vector of numerical values
#' @return Returns a vector
#' @export
rev_cum_sum <- function(x) {
  x <- rev(x)
  return(rev(cumsum(x)))
}

#' Computes initial bound of set of values.
#' @param x Cumulative intensity values.
#' @param left_sigma_ratio_lim Left-standard deviation of the bi-Gaussian function.
#' @return Returns end bound.
#' @export
compute_start_bound <- function(x, left_sigma_ratio_lim) {
  start_bound <- 1  
  len_x <- length(x)
  idx <- which(x >= left_sigma_ratio_lim / (left_sigma_ratio_lim + 1) * x[len_x])
  if (length(idx) > 0) {
    start_bound <- max(1, min(idx))
  }
  return (start_bound)
}

#' Computes final bound of set of values.
#' @param x Cumulative intensity values.
#' @param right_sigma_ratio_lim Right-standard deviation of the bi-Gaussian function.
#' @return Returns end bound.
#' @export
compute_end_bound <- function(x, right_sigma_ratio_lim) {
  len_x <- length(x)
  end_bound <- len_x - 1
  idx <- which(x <= right_sigma_ratio_lim / (right_sigma_ratio_lim + 1) * x[len_x])
  if (length(idx) > 0) {
    end_bound <- min(len_x - 1, max(idx))
  }
  return (end_bound)
}

#' Computes initial and final bounds of set of values.
#' @param x Cumulative intensity values.
#' @param sigma_ratio_lim A vector of two. It enforces the belief of the range of the ratio between the left-standard deviation.
#' and the right-standard deviation of the bi-Gaussian function used to fit the data.
#' @return Returns a list with bounds with following items:
#' \itemize{
#'   \item start - start bound
#'   \item end - end bound
#'}
#' @export
compute_bounds <- function(x, sigma_ratio_lim) {
  start <- compute_start_bound(x, sigma_ratio_lim[1])
  end <- compute_end_bound(x, sigma_ratio_lim[2])
  return(list(start = start, end = end))
}

#' Custom difference bounded to 4-fold minimum difference
#' @description
#' Compute difference between neighbouring elements of a vector and optionally apply a 
#' mask such that the maximum difference is no higher than 4-fold minimum difference.
#' @param x - float - a vector of numerical values.
#' @param apply_mask - boolean - whether to apply threshold mask to the output vector.
#' @return Returns vector of numeric differences between neighbouring values.
#' @export
compute_dx <- function(x, apply_mask=TRUE) {
  l <- length(x)
  diff_x <- diff(x)
  if (l == 2) {
      dx <- rep(diff_x, 2)
  } else {
    dx <- c(
      x[2] - x[1],
      diff(x, lag = 2) / 2,
      x[l] - x[l - 1]
    )
  }
  if (apply_mask) {
    diff_threshold <- min(diff_x) * 4
    dx <- pmin(dx, diff_threshold)
  }
  return (dx)
}

#' Generate chromatographic profile for a feature.
#' @description
#' Find base.curve RTs that lay within RT range of the whole feature table and append intensities to these RTs.
#' @param profile Profile table with shape number-of-features*4 in dataframe.The table contains following columns:
#' \itemize{
#'   \item mz - float - mass-to-charge ratio of feature
#'   \item rt - float - retention time of features
#'   \item intensity - float - intensity of features
#'   \item group_number - integer - group number assigned to each feature based on their rt similarity
#' }
#' @param base.curve Matrix that contains rts of feature in the same rt cluster.
#' @return dataframe with two columns
#' @export
compute_chromatographic_profile <- function(profile, base.curve) {
  rt_range <- range(profile[, "rt"])
  rt_profile <- base.curve %>%
    dplyr::filter(dplyr::between(rt, min(rt_range), max(rt_range))) %>% dplyr::mutate(intensity = 0)
  rt_profile[rt_profile[, "rt"] %in% profile[, "rt"], 'intensity'] <- profile[, "intensity"]
  return (rt_profile)
}

#' Estimates total signal strength (total area of the estimated normal curve).
#' @param y - float - a vector of intensities.
#' @param d - float - a vector of \emph{y} values in a gaussian curve.
#' @return scale - float - a vector of scaled intensity values.
#' @export
compute_scale <- function(y, d) {
  dy_ratio <- d^2 * log(y / d)
  dy_ratio[is.na(dy_ratio)] <- 0
  dy_ratio[is.infinite(dy_ratio)] <- 0

  scale <- exp(sum(dy_ratio) / sum(d^2))
  return(scale)
}

#' Estimate the parameters of Bi-Gaussian curve by Method of Moments
#' @param x Vector of RTs that lay in the same RT cluster.
#' @param y Intensities that belong to x.
#' @param moment_power The parameter for data transformation when fitting the bi-Gaussian or Gaussian mixture model in an EIC.
#' @param sigma_ratio_lim A vector of two. It enforces the belief of the range of the ratio between the left-standard deviation
#'  and the right-standard deviation of the bi-Gaussian function used to fit the data.
#' @return A vector with length 4. The items are as follows going from first to last:
#' \itemize{
#'   \item mean of gaussian curve
#'   \item standard deviation at the left side of the gaussian curve
#'   \item standard deviation at the right side of the gaussian curve
#'   \item estimated total signal strength (total area of the estimated normal curve)
#'}
#' @export
bigauss.esti <- function(x, y, moment_power = 1, do.plot = FALSE, sigma_ratio_lim = c(0.3, 3)) {
  # even producing a dataframe with x and y as columns without actually using it causes the test to run forever
  sel <- which(y > 1e-10)
  if (length(sel) < 2) return (c(median(x), 1, 1, 0))

  x <- x[sel]
  y <- y[sel]

  y.0 <- y
  max.y.0 <- max(y.0, na.rm = TRUE)
  y <- (y / max.y.0)^moment_power

  dx <- compute_dx(x)

  y.cum <- cumsum(y * dx)
  x.y.cum <- cumsum(y * x * dx)
  xsqr.y.cum <- cumsum(y * x^2 * dx)

  y.cum.rev <- rev_cum_sum(y * dx)
  x.y.cum.rev <- rev_cum_sum(x * y * dx)
  xsqr.y.cum.rev <- rev_cum_sum(y * x^2 * dx)

  bounds <- compute_bounds(y.cum, sigma_ratio_lim)
  end <- bounds$end
  start <- bounds$start

  if (end <= start) {
    m <- min(mean(x[start:end]), x[max(which(y.cum.rev > 0))])
  } else {
    m.candi <- x[start:end] + diff(x[start:(end + 1)]) / 2
    rec <- matrix(numeric(0), nrow = 0, ncol = 3)

    s1 <- sqrt((xsqr.y.cum[start:end] + m.candi^2 * y.cum[start:end] - 2 * m.candi * x.y.cum[start:end]) / y.cum[start:end])
    s2 <- sqrt((xsqr.y.cum.rev[start:end + 1] + m.candi^2 * y.cum.rev[start:end + 1] - 2 * m.candi * x.y.cum.rev[start:end + 1]) / y.cum.rev[start:end + 1])
    rec <- rbind(rec, cbind(s1, s2, y.cum[start:end] / y.cum.rev[start:end + 1]))

    d <- log(rec[,1] / rec[,2]) - log(rec[,3])
    if (min(d, na.rm = TRUE) * max(d, na.rm = TRUE) < 0) {
      sel <- c(which(d == max(d[d < 0]))[1], which(d == min(d[d >= 0])))
      m <- (sum(abs(d[sel]) * m.candi[sel])) / (sum(abs(d[sel])))
    } else {
      d <- abs(d)
      m <- m.candi[which(d == min(d, na.rm = TRUE))[1]]
    }
  }

  sel1 <- which(x < m)
  sel2 <- which(x >= m)
  s1 <- sqrt(sum((x[sel1] - m)^2 * y[sel1] * dx[sel1]) / sum(y[sel1] * dx[sel1]))
  s2 <- sqrt(sum((x[sel2] - m)^2 * y[sel2] * dx[sel2]) / sum(y[sel2] * dx[sel2]))

  s1 <- s1 * sqrt(moment_power)
  s2 <- s2 * sqrt(moment_power)

  d1 <- dnorm(x[sel1], sd = s1, mean = m)
  d2 <- dnorm(x[sel2], sd = s2, mean = m)
  d <- c(d1 * s1, d2 * s2) # notice this "density" does not integrate to 1. Rather it integrates to (s1+s2)/2
  y <- y.0

  scale <- compute_scale(y, d)

  if (do.plot) {
    plot(x, y)
    abline(v = m)
    lines(x[y > 0], d * scale, col = "red")
  }

  to.return <- c(m, s1, s2, scale)
  if (sum(is.na(to.return)) > 0) {
    m <- sum(x * y) / sum(y)
    s1 <- s2 <- sum(y * (x - m)^2) / sum(y)
    scale <- sum(y) / s1
    to.return <- c(m, s1, s2, scale)
  }
  
  return(to.return)
}

#' Calculates the three initial bi-gaussian parameters (sd1, sd2, and scaling factor)
#' @param rt_profile A matrix with two columns: "base.curve" (rt) and "intensity".
#' @param pks A vector of sorted RT-peak values at which the kernel estimate was computed.
#' @param vlys A vector of sorted RT-valley values at which the kernel estimate was computed.
#' @param dx Difference between neighbouring RT values with step 2.
#' @return A matrix. The items are as follows going from first to last:
#' \itemize{
#'   \item s1: standard deviation at the left side of the gaussian curve
#'   \item s2: standard deviation at the right side of the gaussian curve
#'   \item delta: estimated total signal strength (total area of the estimated normal curve)
#' }
#' probably a matrix is better though
#' @export
compute_initiation_params <- function(rt_profile, pks, vlys, dx) {
  miu <- s1 <- s2 <- delta <- pks
  for (i in 1:length(miu)) {
    ind.1 <- which(rt_profile[, "rt"] >= max(vlys[vlys < miu[i]]) & rt_profile[, "rt"] < miu[i])
    s1[i] <- sqrt(sum((rt_profile[ind.1, "rt"] - miu[i])^2 * rt_profile[ind.1, "intensity"] * dx[ind.1]) / sum(rt_profile[ind.1, "intensity"] * dx[ind.1]))

    ind.2 <- which(rt_profile[, "rt"] >= miu[i] & rt_profile[, "rt"] < min(vlys[vlys > miu[i]]))
    s2[i] <- sqrt(sum((rt_profile[ind.2, "rt"] - miu[i])^2 * rt_profile[ind.2, "intensity"] * dx[ind.2]) / sum(rt_profile[ind.2, "intensity"] * dx[ind.2]))

    delta[i] <- (sum(rt_profile[ind.1, "intensity"] * dx[ind.1]) + sum(rt_profile[ind.2, "intensity"] * dx[ind.2])) / 
    ((sum(dnorm(rt_profile[ind.1, "rt"], mean = miu[i], sd = s1[i])) * s1[i] / 2) + (sum(dnorm(rt_profile[ind.2, "rt"], mean = miu[i], sd = s2[i])) * s2[i] / 2))
  }
  return (cbind(s1=s1, s2=s2, delta=delta))
}

#' Computes the expectation step of the EM method.
#' @param m A vector of sorted RT-peak values at which the kernel estimate was computed.
#' @param rt_profile A matrix with two columns: "base.curve" (rt) and "intensity".
#' @param delta Parameter computed by the initiation step.
#' @param s1 Parameter computed by the initiation step.
#' @param s2 Parameter computed by the initiation step.
#' @export
compute_e_step <- function(rt_profile, miu, s1, s2, delta) {
  fit <- matrix(numeric(0), ncol = length(miu), nrow = length(rt_profile[, "rt"])) # this is the matrix of fitted values
  cuts <- c(-Inf, miu, Inf)
  for (j in 2:length(cuts)) {
    ind <- which(dplyr::between(rt_profile[, "rt"], cuts[j - 1], cuts[j]))
    use_sd <- ifelse((1:length(miu)) >= (j - 1), s1, s2)
    for (i in 1:length(miu)) fit[ind, i] <- dnorm(rt_profile[ind, "rt"], mean = miu[i], sd = use_sd[i]) * use_sd[i] * delta[i]
  }
  fit[is.na(fit)] <- 0
  return(fit)
}

#' Bi-gaussian mixture model estimation
#' @description
#' Estimates the optimal bi-gaussian parameters using the EM method. It accepts two internal computation of parameters for "moment"
#' and "EM" model input options.
#' @param rt_profile Dataframe that stores RTs and intensities of features.
#' @param moment_power The parameter for data transformation when fitting the bi-Gaussian or Gaussian mixture model in an EIC.
#' @param sigma_ratio_lim A vector of two. It enforces the belief of the range of the ratio between the left-standard deviation
#'  and the right-standard deviation of the bi-Gaussian function used to fit the data.
#' @param bw Bandwidth vector to use in the kernel smoother.
#' @param eliminate When a component accounts for a proportion of intensities less than this value, the component will be ignored.
#' @param max.iter Maximum number of iterations when executing the E step.
#' @param peak_estim_method The estimation method for the bi-Gaussian peak model. Two possible values: moment and EM.
#' @param BIC_factor The factor that is multiplied on the number of parameters to modify the BIC criterion. If larger than 1,
#'  models with more peaks are penalized more.
#' @importFrom dplyr filter arrange
#' @export
bigauss.mix <- function(rt_profile, moment_power = 1, do.plot = FALSE, sigma_ratio_lim = c(0.1, 10), bw = c(15, 30, 60), eliminate = .05, max.iter = 25, peak_estim_method, BIC_factor = 2) {
  results <- new("list")
  all.bw <- sort(bw)
  record.smoother <- setNames(vector("list", length(all.bw)), all.bw)  # record smoothed peaks and valleys
  record.bic <- all.bw  # record BIC for each bandwidth

  rt_profile_unfiltered <- rt_profile  # keep for kernel smoothing for tests but consider removing - gaussian estimation smooths on filtered dataset - not consistent
  rt_profile <- data.frame(rt_profile) |> dplyr::filter(intensity > 1e-5) |> dplyr::arrange(rt)  
  peaks_count <- Inf

  for (ind in length(all.bw):1)
  {
    # kernel smoothing, peak and valley detection
    bw <- all.bw[ind]
    this.smooth <- ksmooth(rt_profile_unfiltered[, "rt"], rt_profile_unfiltered[, "intensity"], kernel = "normal", bandwidth = bw)
    turns <- find.turn.point(this.smooth$y)
    pks <- this.smooth$x[turns$pks]
    vlys <- c(-Inf, this.smooth$x[turns$vlys], Inf)
    record.smoother[[as.character(bw)]] <- list(pks = pks, vlys = vlys)

    results[[ind]] <- NA
    record.bic[ind] <- Inf
    params <- matrix(numeric(0), nrow = 0, ncol = 4, dimnames=list(NULL, c("miu", "s1", "s2", "delta")))

    if (length(pks) != peaks_count) {
      peaks_count <- length(pks)
      dx <- compute_dx(rt_profile[, "rt"], apply_mask = FALSE)

      # initiation
      initiation_params <- compute_initiation_params(rt_profile, pks, vlys, dx)
      params <- rbind(params, cbind(pks, initiation_params[,'s1'], initiation_params[,'s2'], initiation_params[,'delta']))
      params[is.na(params)] <- 1e-10  # why?

      this.change <- Inf
      counter <- 0

      while (this.change > 0.1 && counter <= max.iter) {
        counter <- counter + 1
        old.m <- params[,'miu']

        # E step
        fit <- compute_e_step(rt_profile, params[,'miu'], params[,'s1'], params[,'s2'], params[,'delta'])
        
        # Elimination step
        fit <- fit / apply(fit, 1, sum)
        fit2 <- fit * rt_profile[, "intensity"]
        perc.explained <- apply(fit2, 2, sum) / sum(rt_profile[, "intensity"])
        max.erase <- max(1, round(length(perc.explained) / 5))
        to.erase <- which(perc.explained <= min(eliminate, perc.explained[order(perc.explained, na.last = FALSE)[max.erase]]))

        if (length(to.erase) > 0) {
          old.m <- old.m[-to.erase]
          params <- params[-to.erase, , drop = FALSE]
          fit <- fit[, -to.erase]
          if (is.null(ncol(fit))) fit <- matrix(fit, ncol = 1)
          fit <- fit / apply(fit, 1, sum)
        }

        # M step
        for (i in 1:length(params[,1])) {  
          this.y <- rt_profile[, "intensity"] * fit[, i]
          if (peak_estim_method == "moment") {
            this.fit <- bigauss.esti(rt_profile[, "rt"], this.y, moment_power = moment_power, do.plot = FALSE, sigma_ratio_lim = sigma_ratio_lim)
          } else {
            this.fit <- bigauss.esti.EM(rt_profile[, "rt"], this.y, do.plot = FALSE, sigma_ratio_lim = sigma_ratio_lim)
          }
          params[i, ] <- c(this.fit[1], this.fit[2], this.fit[3], this.fit[4])
        }

        params[is.na(params[, 'delta']), 4] <- 0  # why?
        this.change <- sum((old.m - params[,'miu'])^2)  # amount of change - sum of squared differences  
      }
      # E step again
      fit <- compute_e_step(rt_profile, params[,'miu'], params[,'s1'], params[,'s2'], params[,'delta'])

      if (do.plot) {
        par(mfrow = c(ceiling(length(all.bw) / 2), 2), mar = c(1, 1, 1, 1))
        plot_rt_profile(rt_profile, bw, fit, params[,1])
      }

      area <- params[,'delta'] * (params[,'s1'] + params[,'s2']) / 2
      rss <- sum((rt_profile[, "intensity"] - apply(fit, 1, sum))^2)
      l <- length(rt_profile[, "rt"]) 
      bic <- l * log(rss / l) + 4 * length(params[,'miu']) * log(l) * BIC_factor
      results[[ind]] <- cbind(params, area)
      record.bic[ind] <- bic
    }
  }

  sel <- order(record.bic, -all.bw)[1]
  record <- new("list")
  record$param <- results[[sel]]
  record$record.smoother <- record.smoother
  record$all.param <- results
  record$bic <- record.bic
  return(record)
}

#' Reevaluate parameters of chromatographic gaussian curves.
#' @description
#' Estimates Gaussian Mixture Model parameters with Expectation-Maximization (EM) algorithm.
#' Returns a matrix with columns: miu, sigma, scale, each row corresponds to a peak
#' @param that.curve Dataframe that stores RTs and intensities of features.
#' @param pks A vector of sorted RT-peak values at which the kernel estimate was computed.
#' @param vlys A vector of sorted RT-valley values at which the kernel estimate was computed.
#' @param ignore In fitting mixture of bi-Gaussian (or Gaussian) model of an EIC, when a component accounts for a
#' proportion of intensities less than this value, the component will be ignored.
#' @param max.iter Maximum number of iterations when reevaluating gaussian curves.
#' @param aver_diff Average retention time difference across RTs of all features.
#' @importFrom dplyr between
#' @return A list. The items are as follows going from first to last:
#' \itemize{
#'   \item miu - float - mean value of the gaussian curve
#'   \item sigma - float - standard deviation of the gaussian curve
#'   \item scale - float - estimated total signal strength (total area of the estimated normal curve)
#' }
#' @export
normix <- function(rt_profile, pks, vlys, ignore = 0.1, max.iter = 50, aver_diff) {
  x <- rt_profile[, 'rt']
  y <- rt_profile[, 'intensity']
  colnames <- c("miu", "sigma", "scale")
  params <- matrix(numeric(0), nrow=0, ncol=3, dimnames=list(NULL, colnames))

  if (length(pks) == 0) { # no peaks, does it happen?
    return(params)
  } 

  if (length(pks) == 1) {  # check for 1 peak case
    mu_sc_std <- compute_mu_sc_std(rt_profile, aver_diff)
    params <- rbind(params, cbind(mu_sc_std$miu, mu_sc_std$sigma, mu_sc_std$sc))
    return(params)
  }

  pks <- sort(pks)
  vlys <- sort(vlys)
  l <- length(pks)

  # predict initial parameters for each peak by rt and intensity values between neighbouring valleys
  for (m in 1:l)  {
    indices <- dplyr::between(x, max(vlys[vlys <= pks[m]]), min(vlys[vlys >= pks[m]]))

    if (length(x[indices]) == 0 | length(y[indices]) == 0) {  
      params <- rbind(params, cbind(NA, NA, 1))  # no data points in this region, why is scale=1?????
    } else {
      rt_profile_filt <- data.frame(rt = x[indices], intensity = y[indices])
      mu_sc_std <- compute_mu_sc_std(rt_profile_filt, aver_diff)
      params <- rbind(params, cbind(mu_sc_std$miu, mu_sc_std$sigma, mu_sc_std$sc))
    }
  }
  
  # erase invalid peaks, return record if no peaks left # maybe check to see if compute_mu_sc_std can return NA or 0
  to.erase <- which(is.na(params[,'miu']) | is.na(params[,'sigma']) | params[,'sigma'] == 0 | is.na(params[,'scale']))
  if (length(to.erase) > 0) {
    l <- l - length(to.erase)
    params <- params[-to.erase, , drop = FALSE]
    if (l == 0) return(params)
  }

  diff <- 1000
  counter <- 0
  # expectation-maximization loop
  while (diff > 0.05 & counter < max.iter) {
    counter <- counter + 1
    if (l == 1) {
      mu_sc_std <- compute_mu_sc_std(rt_profile, aver_diff)
      params <- matrix(cbind(mu_sc_std$miu, mu_sc_std$sigma, mu_sc_std$sc), nrow=1, ncol=3, dimnames=list(NULL, colnames))
      return(params)
    }
    
    miu.previous <- params[,'miu']    

    fit <- t(sapply(1:l, function(m) dnorm(x, mean = params[m,'miu'], sd = params[m,'sigma']) * params[m,'scale'])) # estimated Gaussian distributions (y values) for each peak (component) at each RT point
    total.fit <- y * 0
    for (m in 1:l) total.fit <- total.fit + fit[m, ]  # total estimated intensity (Gaussian curve) at each RT point      
    fit.weighted <- t(apply(fit, 1, function(x) x / total.fit)) # normalize by total estimated intensity to get weights for each component at each RT point
    
    if (any(is.na(fit.weighted))) break

    params <- matrix(numeric(0), nrow=0, ncol=3, dimnames=list(NULL, colnames))
    for (m in 1:l)  # estimate new parameters for each peak with weighted intensities
    {
      rt_profile_weighted <- data.frame(rt = x, intensity = y * fit.weighted[m, ])  # estimated Gaussian parameters for fitted weighted intensities at each RT point
      mu_sc_std <- compute_mu_sc_std(rt_profile_weighted, aver_diff)
      params <- rbind(params, cbind(mu_sc_std$miu, mu_sc_std$sigma, mu_sc_std$sc))  
      if (params[m,'sigma'] == 0) params[m,'scale'] <- NA  # why???
    }

    diff <- sum((miu.previous - params[,'miu'])^2)  # squared Euclidean distance between old and new mean values
    weights <- colSums(apply(fit.weighted, 1, function(row) row*y))  # total sum of weighted intensities at each RT point for each peak (component)
    weights[is.na(params[,'scale'])] <- 0  # filter if scale is NA
    weights <- weights / sum(weights)  # normalise by total sum of weighted intensities

    max.erase <- max(1, round(l / 5))
    to.erase <- which(weights <= min(ignore, weights[order(weights, na.last = FALSE)[max.erase]]))
    if (length(to.erase) > 0) {  # erase invalid peaks + their weights, return record if no peaks left
      l <- l - length(to.erase) 
      params <- params[-to.erase, , drop = FALSE]
      if (l == 0) {
        params <- matrix(numeric(0), nrow = 0, ncol = 3, dimnames=list(NULL, colnames))
        return(params)
      }
      diff <- 1000
    }
  }
  return(params)
}

#' Estimates parameters of a gaussian curve with normix and selects the best model using BIC criterion.
#' @description
#' Returns named list containing the parameters of the best model and records of all models.
#' @param rt_profile Dataframe that stores RTs and intensities of features.
#' @param x Vector of RTs that lay in the same RT cluster.
#' @param y Intensities that belong to x.
#' @param bw Bandwidth vector to use in the kernel smoother.
#' @param eliminate When a component accounts for a proportion of intensities less than this value, the component will be ignored.
#' @param max.iter Maximum number of iterations when executing the E step.
#' @param BIC_factor The factor that is multiplied on the number of parameters to modify the BIC criterion. If larger than 1,
#' @param aver_diff Average retention time difference across RTs of all features.
#' @export
normix.bic <- function(rt_profile, do.plot = FALSE, bw = c(15, 30, 60), eliminate = .05, max.iter = 50, BIC_factor = 2, aver_diff) {
  rt_profile <- rt_profile |> dplyr::filter(intensity > 1e-5) |> dplyr::arrange(rt)
  x <- rt_profile[, 'rt']
  y <- rt_profile[, 'intensity']

  results <- new("list")
  all.bw <- sort(bw) # order bandwidths from small to large
  record.bic <- all.bw  # record BIC for each bandwidth
  record.smoother <- setNames(vector("list", length(all.bw)), all.bw)  # record smoothed peaks and valleys
  
  peaks_count <- Inf

  for (ind in length(all.bw):1)
  { 
    # kernel smoothing, peak and valley detection
    bw <- all.bw[ind]
    this.smooth <- ksmooth(x, y, kernel = "normal", bandwidth = bw)
    turns <- find.turn.point(this.smooth$y)
    pks <- this.smooth$x[turns$pks]
    vlys <- c(-Inf, this.smooth$x[turns$vlys], Inf)
    record.smoother[[as.character(bw)]] <- list(pks = pks, vlys = vlys)

    results[[ind]] <- NA 
    record.bic[ind] <- Inf

    if (length(pks) != peaks_count) {  # if a different number of peaks is found after smoothing, reestimate gaussian mixture
      peaks_count <- length(pks)
      
      gaussian.est <- normix(rt_profile, pks = pks, vlys = vlys, ignore = eliminate, max.iter = max.iter, aver_diff = aver_diff)  # estimate gaussian mixture parameters
      if (nrow(gaussian.est) == 0) next  # on to next bandwidth, no valid model returned. maybe use all(is.na(gaussian.est)) to check for empty matrix
            
      total.fit <- x * 0  # compute total fitted values
      for (i in 1:nrow(gaussian.est)) total.fit <- total.fit + dnorm(x, mean = gaussian.est[i, 'miu'], sd = gaussian.est[i, 'sigma']) * gaussian.est[i, 'scale']

      # compute BIC
      rss <- sum((y - total.fit)^2)  # is BIC best criterion? 
      bic <- length(x) * log(rss / length(x)) + 3 * nrow(gaussian.est) * log(length(x)) * BIC_factor
      results[[ind]] <- gaussian.est
      record.bic[ind] <- bic

      if (do.plot) {
        par(mfrow = c(ceiling(length(all.bw) / 2), 2), mar = c(1, 1, 1, 1))
        plot_normix_bic(x, y, bw, gaussian.est)
      }
    } 
  }

  sel <- order(record.bic, -all.bw)[1]
  record <- new("list")
  record$record.smoother <- record.smoother
  record$all.param <- results
  record$bic <- record.bic
  record$param <- results[[sel]]  # is a list/matrix
  return(record)
}

#' Generate feature table from noise-removed LC/MS profile.
#' @description
#' Each LC/MS profile is first processed by the function remove_noise() to remove noise and reduce data size. A matrix containing m/z
#' value, retention time, intensity, and group number is output from remove_noise(). This matrix is then fed to the function
#' prof.to.features() to generate a feature list. Every detected feature is summarized into a single row in the output matrix from this function.
#' @param profile The matrix output from remove_noise(). It contains columns of m/z value, retention time, intensity and group number.
#' @param bandwidth A value between zero and one. Multiplying this value to the length of the signal along the time axis helps
#'  determine the bandwidth in the kernel smoother used for peak identification.
#' @param min_bandwidth The minimum bandwidth to use in the kernel smoother.
#' @param max_bandwidth The maximum bandwidth to use in the kernel smoother.
#' @param sd_cut A vector of two. Features with standard deviation outside the range defined by the two numbers are eliminated.
#' @param sigma_ratio_lim A vector of two. It enforces the belief of the range of the ratio between the left-standard deviation
#'  and the right-standard deviation of the bi-Gaussian function used to fit the data.
#' @param shape_model The mathematical model for the shape of a peak. There are two choices - "bi-Gaussian" and "Gaussian".
#'  When the peaks are asymmetric, the bi-Gaussian is better. The default is "bi-Gaussian".
#' @param peak_estim_method The estimation method for the bi-Gaussian peak model. Two possible values: moment and EM.
#' @param do.plot Whether to generate diagnostic plots.
#' @param moment_power The parameter for data transformation when fitting the bi-Gaussian or Gaussian mixture model in an EIC.
#' @param component_eliminate In fitting mixture of bi-Gaussian (or Gaussian) model of an EIC, when a component accounts for a
#'  proportion of intensities less than this value, the component will be ignored.
#' @param BIC_factor the factor that is multiplied on the number of parameters to modify the BIC criterion. If larger than 1,
#'  models with more peaks are penalized more.
#' @return A matrix is returned. The columns are: m/z value, retention time, spread (standard deviation of the estimated normal
#'  curve), and estimated total signal strength (total area of the estimated normal curve).
#' @export
prof.to.features <- function(profile,
                             bandwidth,
                             min_bandwidth,
                             max_bandwidth,
                             sd_cut,
                             sigma_ratio_lim,
                             shape_model,
                             peak_estim_method,
                             moment_power,
                             component_eliminate,
                             BIC_factor,
                             do.plot) {
  validate_model_method_input(shape_model, peak_estim_method)

  profile <- preprocess_profile(profile)  # returns dframe with columns: mz, rt, intensity, group_number

  bws <- preprocess_bandwidth(min_bandwidth, max_bandwidth, profile)
  min_bandwidth <- bws$min_bandwidth
  max_bandwidth <- bws$max_bandwidth

  base.curve <- data.frame('rt'=sort(unique(profile$rt))) # sorted unique rt values  
  all_diff_mean_rts <- compute_delta_rt(base.curve$rt) # computes diff of mean values from consecutive values 
  aver_diff <- mean(diff(base.curve$rt))  

  keys <- c("mz", "rt", "sd1", "sd2", "area")
  peak_parameters <- matrix(0, nrow = 0, ncol = length(keys), dimnames = list(NULL, keys))
  
  feature_groups <- split(profile, profile$group_number)

  # loop over each group
  for (i in seq_along(feature_groups))
  {
    # init variables
    feature_group <- feature_groups[[i]] |> dplyr::arrange_at("rt") 

    num_features <- nrow(feature_group)
    # The estimation procedure for a single peak
    # Defines the dataframe containing median_mz, median_rt, sd1, sd2, and area
    if (num_features < 2) {
      time_weights <- all_diff_mean_rts[which(base.curve$rt %in% feature_group$rt)]
      rt_peak_shape <- c(feature_group[,'mz'], feature_group[,'rt'], NA, NA, feature_group[,'intensity'] * time_weights)
      peak_parameters <- rbind(peak_parameters, rt_peak_shape)
    } else {
      # find bandwidth for these particular range
      rt_range <- range(feature_group$rt)
      bw <- min(max(bandwidth * (max(rt_range) - min(rt_range)), min_bandwidth), max_bandwidth)
      bw <- seq(bw, 2 * bw, length.out = 3)
      if (bw[1] > 1.5 * min_bandwidth) {
        bw <- c(max(min_bandwidth, bw[1] / 2), bw)
      }

      rt_profile <- compute_chromatographic_profile(feature_group, base.curve)  # returns df with columns: rt, intensity
      if (shape_model == "Gaussian") {
        rt_peak_shape <- compute_gaussian_peak_shape(rt_profile, bw, component_eliminate, BIC_factor, aver_diff)
        } else {
        rt_peak_shape <- bigauss.mix(rt_profile, sigma_ratio_lim = sigma_ratio_lim, bw = bw, moment_power = moment_power, peak_estim_method = peak_estim_method, eliminate = component_eliminate, BIC_factor = BIC_factor)$param[, c(1, 2, 3, 5)]
      }
      if (is.null(nrow(rt_peak_shape))) {  # only one peak is found
        peak_parameters <- rbind(peak_parameters, c(median(feature_group[,'mz']), rt_peak_shape))
      }
      
      else {     
        for (m in 1:nrow(rt_peak_shape))  # multiple peaks
        {
          rt_diff <- abs(feature_group[, "rt"] - rt_peak_shape[m, 'miu'])
          peak_parameters <- rbind(peak_parameters, c(mean(feature_group[which(rt_diff == min(rt_diff)), 'mz']), rt_peak_shape[m, ]))
        }
      }
    }
  }
  peak_parameters <- peak_parameters[order(peak_parameters[, "mz"], peak_parameters[, "rt"]), ]
  peak_parameters <- peak_parameters[which(apply(peak_parameters[, c("sd1", "sd2")], 1, min) > sd_cut[1] & apply(peak_parameters[, c("sd1", "sd2")], 1, max) < sd_cut[2]), ]
  rownames(peak_parameters) <- NULL

  if (do.plot) {
    plot_peak_summary(feature_groups, peak_parameters)
  }

  return(tibble::as_tibble(peak_parameters))
}
