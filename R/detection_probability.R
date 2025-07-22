


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
#' @param plotting boolean flag to indicate whether to make plots
#' @export
#' @return A data frame containing range, number observed and number expected. The number observed is constrained to be <=number expected.
prep_detection_data <- function(target_ranges, vol_func, nbins=25, method='median', nvals=5, loc_dens=NULL, plotting=FALSE){
  bin_edges=seq(min(target_ranges),max(target_ranges),length.out=nbins+1)
  bin_mids=bin_edges[1:nbins]+(bin_edges[2]-bin_edges[1])/2

  vol=numeric(length=length(bin_mids))
  bin_counts=numeric(length=length(bin_mids))
  for (i in 1:length(bin_mids)){
    # this is integrating how much volume is in each bin to get to density
    v=integrate(vol_func,lower = bin_edges[i], upper = bin_edges[i+1])
    vol[i]=v$value
    ind=which(target_ranges>bin_edges[i] & target_ranges<=bin_edges[i+1])
    bin_counts[i]=length(ind)
  }
  dens=bin_counts/vol
  # determine the local_density
  if (is.null(loc_dens)){# user didn't specify
    sort_dens=sort(dens,decreasing=TRUE)
    if (method=='mean'){loc_dens=mean(sort_dens[1:nvals])}
    else if (method=='median'){loc_dens=median(sort_dens[1:nvals])}
    else loc_dens=max(dens)
  }
  # plot density values
  if(plotting){
    plot(bin_mids,dens,xlab='range from camera (m)',ylab='density')
    abline(h=loc_dens,col="red",lty=2)
    legend(max(bin_mids)*0.7, 0.95*max(dens), legend=c("estimated local density"),
           col="red", lty=2, cex=0.8)
  }
  # generate expected count per range interval
  exp_count=loc_dens*vol
  # create data frame
  density_data=data.frame(range=bin_mids,dens,exp_count,obs_count=bin_counts)
  # implement constraint for max observed value
  ind=which(density_data$obs_count>density_data$exp_count)
  density_data$obs_count[ind]=density_data$exp_count[ind]

  return(density_data)
}

#' density function fitting. Models specified with a method flag.
#' @param density_data vector of ranges for a given target
#' @param method flag specifying which model to use. Valid choices include 'logistic glm', ADD
#' @param formula optional formula to pass to the model, otherwise defaults are used
#' @param stepAIC boolean flag to do backwards step AIC selection or not
#' @param plotting boolean flag to indicate whether to make plots
#' @export
#' @return A R model output structure

fit_density_function = function(density_data,
                                method=c('logistic glm','logistic gam'),
                                formula=NULL,
                                dostepAIC=TRUE,
                                plotting=FALSE){
  method <- match.arg(method)
  # do the binning bit
  if (method=='logistic glm'){
    if(is.null(formula))
      formula <- cbind(obs_count, exp_count - obs_count) ~ range + I(range^2) + I(range^3) +I(range^4)
    out <- glm(formula = formula, data = density_data,
               family = binomial(link="logit"))
  } else if(method=='logistic gam'){
    formula <- cbind(obs_count, exp_count - obs_count) ~ s(range)
    out <- gam::gam(formula = formula, data = density_data,
                  family = binomial(link="logit"))
  }
  if(dostepAIC) out <- step(out)
  detect.function <- function(range) {
    predict(out, newdata=data.frame(range=range), type='response')
  }

  if (plotting){
    plot(density_data$range,density_data$obs_count/density_data$exp_count,xlab='range from camera (m)',ylab='probability of detection')
    x=seq(0, max(density_data$range),length.out=100)
    y=detect.function(x)
    lines(x,y,col="red")
    legend(max(x)*0.7, 0.95*max(y), legend=c("probability function"),
           col="red", lty=1, cex=0.8)
  }
  return(list(model=out, detect.function=detect.function))
}

