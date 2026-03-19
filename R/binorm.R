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
preprocess_bandwidth.new <- function(min_bandwidth, max_bandwidth, profile) {
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

#' Reverse cumulative sum
#' @description
#' Computes vector of cumulative sums on reversed input. Returns cumulative sum vector going from the sum of all elements to one.
#' @param x float - vector of numerical values
#' @return Returns a vector
#' @export
rev_cum_sum.new <- function(x) {
  x <- rev(x)
  return(rev(cumsum(x)))
}

#' Computes initial bound of set of values.
#' @param x Cumulative intensity values.
#' @param left_sigma_ratio_lim Left-standard deviation of the bi-Gaussian function.
#' @return Returns end bound.
#' @export
compute_start_bound.new <- function(x, left_sigma_ratio_lim) {
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
compute_end_bound.new <- function(x, right_sigma_ratio_lim) {
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
compute_bounds.new <- function(x, sigma_ratio_lim) {
  start <- compute_start_bound.new(x, sigma_ratio_lim[1])
  end <- compute_end_bound.new(x, sigma_ratio_lim[2])
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
compute_dx.new <- function(x, apply_mask=TRUE) {
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
compute_chromatographic_profile.new <- function(profile, base.curve) {
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
compute_scale.new <- function(y, d) {
  dy_ratio <- d^2 * log(y / d)
  dy_ratio[is.na(dy_ratio)] <- 0
  dy_ratio[is.infinite(dy_ratio)] <- 0
  
  scale <- exp(sum(dy_ratio) / sum(d^2))
  scale <- sum(y * d) / sum(d^2)
  return(scale)
}

#' Estimate the parameters of Bi-Gaussian curve by Method of Moments
#' @param x Vector of RTs that lay in the same RT cluster.
#' @param y Fitted peak values scaled by observed intensities that correspond to the rt values in x.
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
bigauss.esti.new <- function(x, y, moment_power = 1, do.plot = FALSE, sigma_ratio_lim = c(0.3, 3)) {
  # even producing a dataframe with x and y as columns without actually using it causes the test to run forever
  sel <- which(y > 1e-10)
  if (length(sel) < 2) return (c(median(x), 1, 1, 0))
  
  x <- x[sel]
  y <- y[sel]
  
  y.0 <- y
  max.y.0 <- max(y.0, na.rm = TRUE)
  y <- (y / max.y.0)^moment_power
  
  dx <- compute_dx.new(x)
  # why ?
  y.cum <- cumsum(y * dx)
  x.y.cum <- cumsum(y * x * dx)
  xsqr.y.cum <- cumsum(y * x^2 * dx)
  
  y.cum.rev <- rev_cum_sum.new(y * dx)
  x.y.cum.rev <- rev_cum_sum.new(x * y * dx)
  xsqr.y.cum.rev <- rev_cum_sum.new(y * x^2 * dx)
  
  bounds <- compute_bounds.new(y.cum, sigma_ratio_lim)
  end <- bounds$end
  start <- bounds$start
  
  # miu estimation
  if (end <= start) {
    miu <- min(mean(x[start:end]), x[max(which(y.cum.rev > 0))])
  } else {
    m.candi <- x[start:end] + diff(x[start:(end + 1)]) / 2  # candidate miu values
    rec <- matrix(numeric(0), nrow = 0, ncol = 3)
    
    s1 <- sqrt((xsqr.y.cum[start:end] + m.candi^2 * y.cum[start:end] - 2 * m.candi * x.y.cum[start:end]) / y.cum[start:end])
    s2 <- sqrt((xsqr.y.cum.rev[start:end + 1] + m.candi^2 * y.cum.rev[start:end + 1] - 2 * m.candi * x.y.cum.rev[start:end + 1]) / y.cum.rev[start:end + 1])
    rec <- rbind(rec, cbind(s1, s2, y.cum[start:end] / y.cum.rev[start:end + 1]))
    # save estimated sigmas based on candidate miu
    
    d <- log(rec[,1] / rec[,2]) - log(rec[,3])  # compute log difference of variances - log(ratio)
    if (min(d, na.rm = TRUE) * max(d, na.rm = TRUE) < 0) {
      sel <- c(which(d == max(d[d < 0]))[1], which(d == min(d[d >= 0])))
      miu <- (sum(abs(d[sel]) * m.candi[sel])) / (sum(abs(d[sel])))
    } else {
      d <- abs(d)
      miu <- m.candi[which(d == min(d, na.rm = TRUE))[1]]
    }
  }
  
  sel1 <- which(x < miu)
  sel2 <- which(x >= miu)
  s1 <- sqrt(sum((x[sel1] - miu)^2 * y[sel1] * dx[sel1]) / sum(y[sel1] * dx[sel1]))  # as in paper
  s2 <- sqrt(sum((x[sel2] - miu)^2 * y[sel2] * dx[sel2]) / sum(y[sel2] * dx[sel2]))  # as in paper
  
  s1 <- s1 * sqrt(moment_power)
  s2 <- s2 * sqrt(moment_power)
  
  d <- dbinorm(x, A=1, miu = miu, sigma1 = s1, sigma2=s2)  # this density should integrate to 1
  y <- y.0
  delta <- compute_scale.new(y, d)
  
  to.return <- c(miu, s1, s2, delta)
  if (sum(is.na(to.return)) > 0) {  
    miu <- sum(x * y) / sum(y)
    s1 <- s2 <- sum(y * (x - miu)^2) / sum(y)  # what
    delta <- sum(y) / s1  # what
    to.return <- c(miu, s1, s2, delta)
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
compute_initiation_params_new <- function(rt_profile, pks, vlys, apex, dx) {
  miu <- s1 <- s2 <- delta <- pks
  for (i in 1:length(miu)) {
    ind.1 <- which(rt_profile[, "rt"] >= max(vlys[vlys < miu[i]]) & rt_profile[, "rt"] < miu[i])
    s1[i] <- sqrt(sum((rt_profile[ind.1, "rt"] - miu[i])^2 * rt_profile[ind.1, "intensity"] * dx[ind.1]) / sum(rt_profile[ind.1, "intensity"] * dx[ind.1]))
    
    ind.2 <- which(rt_profile[, "rt"] >= miu[i] & rt_profile[, "rt"] < min(vlys[vlys > miu[i]]))
    s2[i] <- sqrt(sum((rt_profile[ind.2, "rt"] - miu[i])^2 * rt_profile[ind.2, "intensity"] * dx[ind.2]) / sum(rt_profile[ind.2, "intensity"] * dx[ind.2]))
    
    # delta[i] <- (sum(rt_profile[ind.1, "intensity"] * dx[ind.1]) + sum(rt_profile[ind.2, "intensity"] * dx[ind.2])) / 
    # ((sum(dnorm(rt_profile[ind.1, "rt"], mean = miu[i], sd = s1[  i])) * s1[i] / 2) + (sum(dnorm(rt_profile[ind.2, "rt"], mean = miu[i], sd = s2[i])) * s2[i] / 2))
    delta[i] <- apex[i] * sqrt(pi / 2) * (s1[i] + s2[i])
    # why 
  }
  return (cbind(s1=s1, s2=s2, delta=delta))
}

#' Fits the rt values to binormal model, scaled for intensity
#' A should be the estimated delta, better investigate
#' miu, sigma1, sigma2 are all ok estimates
dbinorm <- function(rt, A, miu, sigma1, sigma2) {
  norm <- 2 / (sqrt(2 * pi) * (sigma1 + sigma2))  # normalisation constant for probability distribution
  sigma <- ifelse(rt < miu, sigma1, sigma2)
  # cat(sprintf("Height of curve = %g\n", H * norm))
  A * norm * exp(- (rt - miu)^2 / (2 * sigma^2))
}

#' Computes the expectation step of the EM method. no.
#' Computes the fitted values of the binormal distribution.
#' Each column should be the fitted distribution values to miu, sd and delta. 
#' @param miu A vector of sorted RT-peak values at which the kernel estimate was computed - the peak mode location.
#' @param rt_profile A matrix with two columns: "base.curve" (rt) and "intensity".
#' @param delta Parameter computed by the initiation step.
#' @param s1 Parameter computed by the initiation step.
#' @param s2 Parameter computed by the initiation step.
#' @return A matrix of the fitted values to the distribution. 
#' @export
compute_e_step_new <- function(rt_profile, miu, s1, s2, delta) {
  fit <- matrix(numeric(0), ncol = length(miu), nrow = length(rt_profile[, "rt"])) # this is the matrix of fitted values
  cuts <- c(-Inf, miu, Inf)
  for (j in 2:length(cuts)) {
    ind <- which(dplyr::between(rt_profile[, "rt"], cuts[j - 1], cuts[j]))
    for (i in 1:length(miu)) fit[ind, i] <- dbinorm(rt_profile[ind, "rt"], A=delta[i], miu = miu[i], sigma1 = s1[1], sigma2=s2[i])
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
bigauss.mix.new <- function(rt_profile, moment_power = 1, do.plot = FALSE, sigma_ratio_lim = c(0.1, 10), bw = c(15, 30, 60), eliminate = .05, max.iter = 25, peak_estim_method, BIC_factor = 2) {
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
    turns <- find.turn.point(this.smooth$y)  # OFF! for large bandwidths
    pks <- this.smooth$x[turns$pks]
    apex <- this.smooth$y[turns$pks]
    vlys <- c(-Inf, this.smooth$x[turns$vlys], Inf)
    record.smoother[[as.character(bw)]] <- list(pks = pks, vlys = vlys)
    
    results[[ind]] <- NA
    record.bic[ind] <- Inf
    params <- matrix(numeric(0), nrow = 0, ncol = 4, dimnames=list(NULL, c("miu", "s1", "s2", "delta")))
    
    if (length(pks) != peaks_count) {
      peaks_count <- length(pks)
      dx <- compute_dx.new(rt_profile[, "rt"], apply_mask = FALSE)
      
      # initiation
      initiation_params <- compute_initiation_params_new(rt_profile, pks, vlys, apex, dx)
      params <- rbind(params, cbind(pks, initiation_params[,'s1'], initiation_params[,'s2'], initiation_params[,'delta']))  # first estimate of miu is peak maximum loc
      params[is.na(params)] <- 1e-10  # why?
      
      this.change <- Inf
      counter <- 0
      
      while (this.change > 0.1 && counter <= max.iter) {
        counter <- counter + 1
        old.m <- params[,'miu']
        old.params <- params
        
        # E step - find fitted peak values scaled by intensity
        fit <- compute_e_step_new(rt_profile, params[,'miu'], params[,'s1'], params[,'s2'], params[,'delta'])  # matrix columns are the fitted values

        # Elimination step - why ?
        fit <- fit / apply(fit, 1, sum)  # Qj ??
        # normalise each values by sum of fitted values in every rt point, so that the sum of fitted values in every rt is 1
        fit2 <- fit * rt_profile[, "intensity"]  # multiply proportions by intensity
        perc.explained <- apply(fit2, 2, sum) / sum(rt_profile[, "intensity"])  # how much data is explained by the fitted scaled component compared to observed intensity - q_ij
        # sum column 'intensities multiplied by proportions' and divide by total intensity
        max.erase <- max(1, round(length(perc.explained) / 5))  # why this number. kinda arbitrary
        # eliminate components which don't explain data under certain percetnage threshold
        to.erase <- which(perc.explained <= min(eliminate, perc.explained[order(perc.explained, na.last = FALSE)[max.erase]]))  # indices of components to remove
                
        if (length(to.erase) > 0) {  # remove the components
          old.m <- old.m[-to.erase]
          params <- params[-to.erase, , drop = FALSE]
          fit <- fit[, -to.erase]
          if (is.null(ncol(fit))) fit <- matrix(fit, ncol = 1)
          fit <- fit / apply(fit, 1, sum)  # divide so that sum of fitted components at each rt is 1 
        }
        
        # M step: reestimate parameters of binorm
        for (i in 1:length(params[,1])) {  
          this.y <- rt_profile[, "intensity"] * fit[, i]  # probably wrong, need to be scaled by delta
          if (peak_estim_method == "moment") {
            this.params <- bigauss.esti.new(rt_profile[, "rt"], this.y, moment_power = moment_power, do.plot = FALSE, sigma_ratio_lim = sigma_ratio_lim)
          } else {
            this.params <- bigauss.esti.new.EM(rt_profile[, "rt"], this.y, do.plot = FALSE, sigma_ratio_lim = sigma_ratio_lim)
          }
          params[i, ] <- c(this.params[1], this.params[2], this.params[3], this.params[4])
        }

        params[is.na(params[, 'delta']), 4] <- 0  # why?
        this.change <- sum((old.m - params[,'miu'])^2)  # amount of change - sum of squared differences - higher penalisation for differences thanks to square
        # why not count change for all parameters?
      }
      # E step again
      fit <- compute_e_step_new(rt_profile, params[,'miu'], params[,'s1'], params[,'s2'], params[,'delta'])
      
      if (do.plot) {
        par(mfrow = c(ceiling(length(all.bw) / 2), 2), mar = c(1, 1, 1, 1))
        plot_rt_profile(rt_profile, bw, fit, params[,1])
      }
      
      area <- params[,'delta']
      rss <- sum((rt_profile[, "intensity"] - apply(fit, 1, sum))^2)
      l <- length(rt_profile[, "rt"]) 
      bic <- l * log(rss / l) + 4 * length(params[,'miu']) * log(l) * BIC_factor      
      # not sure, penalises mainly on length of params. why not just penalise based on rss ? 
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
  
  # cat(print(rt_profile), sep = "\n", file = "rt_profile.txt", append = TRUE)
  # print(record$param)
  # cat(record$param, sep = "\n", file = "param.txt", append = TRUE)
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
binorm <- function(profile,
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
  
  profile <- data.frame(profile)
  colnames(profile) <- c("mz", "rt", "intensity", "group_number")
  
  bws <- preprocess_bandwidth.new(min_bandwidth, max_bandwidth, profile)
  min_bandwidth <- bws$min_bandwidth
  max_bandwidth <- bws$max_bandwidth
  
  base.curve <- data.frame('rt'=sort(unique(profile$rt))) # sorted unique rt values  
  all_diff_mean_rts <- compute_delta_rt(base.curve$rt) # computes diff of mean values from consecutive values 
  aver_diff <- mean(diff(base.curve$rt))  
  
  keys <- c("mz", "rt", "sd1", "sd2", "area")
  peak_parameters <- matrix(0, nrow = 0, ncol = length(keys), dimnames = list(NULL, keys))
  
  feature_groups <- split(profile, profile$group_number)
  
  # loop over each group
  results <- vector("list", length(feature_groups))
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
      
      rt_profile <- compute_chromatographic_profile.new(feature_group, base.curve)  # returns df with columns: rt, intensity
      if (shape_model == "Gaussian") {
        rt_peak_shape <- compute_gaussian_peak_shape(rt_profile, bw, component_eliminate, BIC_factor, aver_diff)
      } else {
        rt_peak_shape <- bigauss.mix.new(rt_profile, sigma_ratio_lim = sigma_ratio_lim, bw = bw, moment_power = moment_power, peak_estim_method = peak_estim_method, eliminate = component_eliminate, BIC_factor = BIC_factor)$param  #[, c(1, 2, 3, 5)]
        results[[i]] <- list(
          run_id = i,
          timestamp = Sys.time(),
          profile = rt_profile,
          params = as_tibble(rt_peak_shape)
        )
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
  library(purrr)
  
  results_tbl <- tibble(
    run_id   = map_int(results, "run_id"),
    timestamp= map_chr(results, ~ as.character(.x$timestamp)),
    profile  = map(results, "profile"),
    params   = map(results, "params")
  )
  browser()
  out_file <- "r_output/all_results.feather"
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  write_feather(results_tbl, out_file, compression = "zstd")
  
  return(tibble::as_tibble(peak_parameters))
}
