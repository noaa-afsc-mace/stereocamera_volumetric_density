#' Dummy function
#' @param x Dummy
#' @export
dummy <- function(x) x

#' @param x input for range from camera in m
#' @param f_vol function that defines change in volume with range from camera
#' @param f_detect function that defines probability of detection with range from camera
#' @export
# define effective volume function
eff_vol_func<- function(x, f_vol, f_detect){
  f_vol(x)*f_detect(x)
}

#' @param x input for range from camera.
#' @param cov_factor covariance factor. Can be a single value or vector same length as x
#' @param f_vol function that defines change in volume with range from camera
#' @param f_detect_cov function that defines probability of detection with range from camera
#' @export
# define effective volume function
eff_vol_func_covariate<- function(x, cov_factor, f_vol, f_detect_cov){
  f_vol(x)*f_detect_cov(x, cov_factor)
}


#' Bin the data for fitting the model.  This involves
#' estimation the volume for each range bin as well as how many
#' animals are in the bin. The resulting density is used to estimate steady state local density, which is then used to compute the number
#' of animals that are expected in each range interval. The expected / observed is the basis for determining probability of detection.
#' @param targets vector of ranges for a given target
#' @param vol_func Volume function as returned by get_vol_func
#' @param nbins number of bins to use for density histogram
#' @param method method to determine the local density based on change in density by range.  Valid choices are 'median' and 'mean'
#' @param nvals an integer specifying how many of the highest density values to average for local density.
#' @param loc_dens a value that can be entered to provide a separate externally estimated local density. Default is NULL.
#' @param model_method flag specifying which model to use. Valid choices include 'logistic glm', ADD
#' @param formula optional formula to pass to the model, otherwise defaults are used
#' @param boot_n number of bootstrap iterations
#' @param percentiles two element vector between 0 and 1 indicating the lower and upper confidence intervals
#' @param plotting boolean flag to indicate whether to make plots
#' @export
#' @return A list, including vector of effective volumes, the model object, and the confidence intervals
bootstrap_effective_volume <- function(target_ranges, vol_func, nbins=25, method='median', nvals=5, loc_dens=NULL, model_method='logistic glm', formula=NULL, boot_n=500, percentiles=c(0.025,0.975), plotting=FALSE){


  # now for the bootstrap

  eff_vol=vector()
  loc_dens_vec=vector()
  resample_bin=matrix(nrow=boot_n,ncol=nbins)
  for (i in 1:boot_n){# bootstrap loop
    # do the resampling
    target_ranges_res=target_ranges[sample(length(target_ranges), size = length(target_ranges), replace = TRUE)]

    # prepare data based on resample
    prep_out=prep_detection_data(target_ranges=target_ranges_res,
                                       vol_func=vol_func, nbins=nbins, method=method, nvals=nvals, loc_dens=loc_dens, plotting=FALSE)
    # fit the model
    detection_out <- fit_density_function(prep_out$data,method=model_method,formula=formula, plotting=FALSE)

    est=detection_out$detect.function(prep_out$data$range)
    resample_bin[i,]=est
    # for sanity check, compute the local density
    loc_dens_vec[i]=prep_out$loc_dens

    # integrate
    eff_vol[i]=integrate(eff_vol_func,lower =0, upper =max(target_ranges_res),vol_func, detection_out$detect.function)$value

  }
  # get basic data for plotting
  prep_out=prep_detection_data(target_ranges=target_ranges,
                               vol_func=vol_func, nbins=nbins, method=method, nvals=nvals, loc_dens=NULL, plotting=FALSE)
  detection_out <- fit_density_function(prep_out$data,method=model_method,formula=formula, plotting=plotting)
  #plot(detection_data$range,detection_data$obs_count/detection_data$exp_count)
  upper=vector()
  lower=vector()
  for (i in 1:length(prep_out$data$range)){
    v=quantile(resample_bin[,i], percentiles, na.rm=TRUE)
    lower[i]=v[1]
    upper[i]=v[2]
  }
  if (plotting){
    lines(prep_out$data$range,lower,col="red",lwd=1,lty=2)
    lines(prep_out$data$range,upper,col="red",lwd=1,lty=2)
  }
  fit=y=detection_out$detect.function(prep_out$data$range)
  return(list(eff_vol=eff_vol, lower_CI=lower, upper_CI=upper, fit=fit))
}
