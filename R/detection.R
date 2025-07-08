

#' Bin the data for fitting the model.  This involves
#' estimation the volume for each range bin as well as how many
#' animals are in the bin.
#' @param targets testing
#' @param vol_func Volume function as returned by
#' @param nbins to do
#' @export
#' @return A vector of densities for each bin
get_dens <- function(targets, vol_func, nbins=25){
  bin_edges=seq(min(targets$RANGE),max(targets$RANGE),length.out=nbins+1)
  bin_mids=bin_edges[1:nbins]+(bin_edges[2]-bin_edges[1])/2

  vol=numeric(length=length(bin_mids))
  bin_counts=numeric(length=length(bin_mids))
  for (i in 1:length(bin_mids)){
    # this is integrating how much volume is in each bin to get to density
    v=integrate(vol_func,lower = bin_edges[i], upper = bin_edges[i+1])
    vol[i]=v$value
    ind=which(targets$RANGE>bin_edges[i] & targets$RANGE<=bin_edges[i+1])
    bin_counts[i]=length(ind)
  }
  dens=bin_counts/vol
  return(list(bin_counts=bin_counts, bin_mids=bin_mids, vol=vol, dens=dens))
}

#' single logistic function
#' @param parm Vector of parameters for single logistic (scale,r50,slope)
#' @param x vector of range values
#' @export
func_single_logistic = function(parm,x){parm[1]/(1+9^((parm[2]-x)/parm[3]))}

#' single logistic evaluation
#' @param parm Vector of parameters for single logistic (scale,r50,slope)
#' @param vol vector of volumes of each range interval
#' @param rng vector of midpoints of each range interval
#' @param fishcount vector of numbers of fish in each range interval
#' @export
eval_single_logistic = function(parm,vol,rng,fishcount){
  like=numeric(length=length(rng))
  for (i in 1:length(rng)){
    # expected number of fish per m3
    d_hat=func_single_logistic(parm,rng[i])
    # poisson likelihood?
    # ll=likpois(round(fishcount[i]/vol[i]), round(d_hat), log = TRUE)*fishcount[i]*-1
    # normal
    ll=dnorm(fishcount[i]/vol[i], d_hat, log = TRUE)*fishcount[i]*-1
    # key part here is to multiply the likelihood by the number of observations in the bin.
    # this way bins with more fish in them drive the fit
    like[i]=ll
  }
  return(sum(like))
}

#' single logistic fitting - convenience function that rolls all three steps of fitting a function into one
#' @param target_ranges vector of ranges for a given target
#' @param vol_func volume function as derived from the get_volume method
#' @param nbins number of bins to use for density histogram
#' @param plotting boolean flag to indicate whether to make plots
#' @param start_params Vector of parameters for single logistic (scale,r50,slope)
#' @export
fit_single_logistic = function(target_ranges,vol_func,nbins,plotting=FALSE,start_params=NULL){
  # do the binning bit
  range_bin_edges=seq(min(target_ranges),max(target_ranges),length.out=nbins+1)
  range_bin_mids=range_bin_edges[1:nbins]+(range_bin_edges[2]-range_bin_edges[1])/2

  vol=numeric(length=length(range_bin_mids))
  bin_counts=numeric(length=length(range_bin_mids))
  for (i in 1:length(range_bin_mids)){
    # this is integrating how much volume is in each bin to get to density
    v=integrate(vol_func,lower = range_bin_edges[i], upper = range_bin_edges[i+1])
    vol[i]=v$value
    ind=which(target_ranges>range_bin_edges[i] & target_ranges<=range_bin_edges[i+1])
    bin_counts[i]=length(ind)
  }
  dens=bin_counts/vol
  if (plotting){
    plot(range_bin_mids,dens)
  }

  #run minimization
  # good starting params for single logistic
  if (is.null(start_params)){
    dens_sort=sort(dens,decreasing=TRUE)
    scale=mean(dens_sort[1:5])
    r50 = max(target_ranges)/2
    slope=-1
    start_params=c(scale,r50,slope)}


  out <- optim(par=start_params, fn = eval_single_logistic, vol=vol,
               rng=range_bin_mids,fishcount=bin_counts)
  if (plotting){
    x=seq(0, max(range_bin_edges),length.out=100)
    y=func_single_logistic(out$par,x)
    lines(x,y,col="red")
  }
  return(out)
}

#' double logistic function
#' @param parm Vector of parameters for single logistic (scale,r50,slope)
#' @param x vector of range values
#' @export
func_double_logistic = function(parm,x){(parm[1]*(1/(1+9^((parm[2]-x)/parm[3])))*(1/(1+9^((parm[4]-x)/parm[5]))))}


#' double logistic evaluation
#' @param parm Vector of parameters for double logistic scale,r50_asc,slope_asc,r50_desc,slope_desc)
#' @param vol vector of volumes of each range interval
#' @param rng vector of midpoints of each range interval
#' @param fishcount vector of numbers of fish in each range interval
#' @export
eval_double_logistic = function(parm,vol,rng,fishcount){
  like=numeric(length=length(rng))
  for (i in 1:length(rng)){
    # expected number of fish per m3
    # this functin has a scale (param 1) and an rising and descending logistic component
    d_hat=func_double_logistic(parm,rng[i])
    # poisson likelihood?
    # ll=likpois(round(fishcount[i]/vol[i]), round(d_hat), log = TRUE)*fishcount[i]*-1
    # normal
    ll=dnorm(fishcount[i]/vol[i], d_hat, log = TRUE)*fishcount[i]*-1
    # key part here is to multiply the likelihood by the number of observations in the bin.
    # this way bins with more fish in them drive the fit
    like[i]=ll
  }
  return(sum(like))
}

#' double logistic fitting - convenience function that rolls all three steps of fitting a function into one
#' @param target_ranges vector of ranges for a given target
#' @param vol_func volume function as derived from the get_volume method
#' @param nbins number of bins to use for density histogram
#' @param plotting boolean flag to indicate whether to make plots
#' @param start_params Vector of parameters for single logistic (scale,r50,slope)
#' @export
fit_double_logistic = function(target_ranges,vol_func,nbins,plotting=FALSE,start_params=NULL){
  # do the binning bit
  range_bin_edges=seq(min(target_ranges),max(target_ranges),length.out=nbins+1)
  range_bin_mids=range_bin_edges[1:nbins]+(range_bin_edges[2]-range_bin_edges[1])/2

  vol=numeric(length=length(range_bin_mids))
  bin_counts=numeric(length=length(range_bin_mids))
  for (i in 1:length(range_bin_mids)){
    # this is integrating how much volume is in each bin to get to density
    v=integrate(vol_func,lower = range_bin_edges[i], upper = range_bin_edges[i+1])
    vol[i]=v$value
    ind=which(target_ranges>range_bin_edges[i] & target_ranges<=range_bin_edges[i+1])
    bin_counts[i]=length(ind)
  }
  dens=bin_counts/vol
  if (plotting){
    plot(range_bin_mids,dens)
  }

  #run minimization
  # good starting params for single logistic
  if (is.null(start_params)){
    dens_sort=sort(dens,decreasing=TRUE)
    scale=mean(dens_sort[1:5])
    r501 = max(target_ranges)/3
    slope1=1
    r502 = max(target_ranges)/2
    slope2=-1
    start_params=c(scale,r501,slope1,r502,slope2)}


  out <- optim(par=start_params, fn = eval_double_logistic, vol=vol,
               rng=range_bin_mids,fishcount=bin_counts)
  if (plotting){
    x=seq(0, max(range_bin_edges),length.out=100)
    y=func_double_logistic(out$par,x)
    lines(x,y,col="red")
  }
  return(out)
}

#' normal function
#' @param parm Vector of parameters for single logistic (scale,r50,slope)
#' @param x vector of range values
#' @export
func_normal = function(parm,x){parm[1]*exp(-((x-parm[2])^2/(2*parm[3]^2)))}


#' normal evaluation
#' @param parm vector of 3 parameters: scale, mean, and standard deviation of a normal curve
#' @param vol vector of volumes of each range interval
#' @param rng vector of midpoints of each range interval
#' @param fishcount vector of numbers of fish in each range interval
#' @export
eval_normal = function(parm,vol,rng,fishcount){
  like=numeric(length=length(rng))
  for (i in 1:length(rng)){
    # expected number of fish per m3
    # this function follows a normal distribution but is scaled
    d_hat=func_normal(parm,rng[i])
    # poisson likelihood?
    # ll=likpois(round(fishcount[i]/vol[i]), round(d_hat), log = TRUE)*fishcount[i]*-1
    # normal
    ll=dnorm(fishcount[i]/vol[i], d_hat, log = TRUE)*fishcount[i]*-1
    # key part here is to multiply the likelihood by the number of observations in the bin.
    # this way bins with more fish in them drive the fit
    like[i]=ll
  }
  return(sum(like))
}
