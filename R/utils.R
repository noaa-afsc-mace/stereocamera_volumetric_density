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
#' @param plotting boolean flag to indicate whether to make plots
#' @export
#' @return A data frame containing range, number observed and number expected. The number observed is constrained to be <=number expected.
bootstrap_effective_volume <- function(target_ranges, vol_func, nbins=25, method='median', nvals=5, loc_dens=NULL, model_method='logistic glm', formula=NULL, boot_n=500, plotting=FALSE){


  # now for the bootstrap

  eff_vol=vector()
  loc_dens_vec=vector()
  resample_bin=matrix(nrow=boot_n,ncol=nbins)
  for (i in 1:boot_n){# bootstrap loop
    # do the resampling
    target_ranges_res=target_ranges[sample(length(target_ranges), size = length(target_ranges), replace = TRUE)]

    # prepare data based on resample
    detection_data=prep_detection_data(target_ranges=target_ranges_res,
                                       vol_func=vol$vol_func, nbins=nbins, method=method, nvals=nvals, loc_dens=loc_dens, plotting=FALSE)
    # fit the model
    out <- fit_density_function(detection_data,method=model_method,formula=formula, dostepAIC=FALSE, plotting=FALSE)

    est=out$detect.function(detection_data$range)
    resample_bin[i,]=est
    # for sanity check, compute the local density
    sort_dens=sort(detection_data$dens,decreasing=TRUE)
    if (method=='mean'){loc_dens_res=mean(sort_dens[1:nvals])}
    else if (method=='median'){loc_dens_res=median(sort_dens[1:nvals])}
    else loc_dens_res=max(dens)
    loc_dens_vec[i]=loc_dens_res

    # estimate effective volume
    eff_vol_func<- function(x){
      vol$vol_func(x)*out$detect.function(x)
    }

    # save vectors into matrix for CI plotting

    # integrate
    eff_vol[i]=integrate(eff_vol_func,lower =0, upper =max(target_ranges_res))$value

  }
  if (plotting){
    # get basic data for plotting
    detection_data=prep_detection_data(target_ranges=target_ranges,
                                       vol_func=vol$vol_func, nbins=nbins, method=method, nvals=nvals, loc_dens=NULL, plotting=FALSE)

    out <- fit_density_function(detection_data,method=model_method,formula=formula, dostepAIC=FALSE, plotting=TRUE)
    #plot(detection_data$range,detection_data$obs_count/detection_data$exp_count)
    upper=vector()
    lower=vector()
    for (i in 1:length(detection_data$range)){
      v=quantile(resample_bin[,i], c(0.025,0.975),na.rm=TRUE)
      lower[i]=v[1]
      upper[i]=v[2]
    }
    lines(detection_data$range,lower,col="red",lwd=1,lty=2)
    lines(detection_data$range,upper,col="red",lwd=1,lty=2)
  }
  return(eff_vol)
}
