

#' Bin the data for fitting the model.  This involves
#' estimation the volume for each range bin as well as how many
#' animals are in the bin.
#' @param targets to do
#' @param vol_func Volume function as returned by \link{\code{get_vol_func}}
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
  return(list(bin_counts=bin_counts, bin_mids=bin_mids, vol=vol))
}




#' logistic
#' @param parm Vector of parameters for logistic
#' @param vol to do
#' @param rng to do
#' @param fishcount to do
#' @export
detect_single_logistic = function(parm,vol,rng,fishcount){
  like=numeric(length=length(rng))
  for (i in 1:length(rng)){
    # expected number of fish per m3
    d_hat=parm[1]/(1+9^((parm[2]-rng[i])/parm[3]))
    # poisson likelihood?
    # ll=likpois(round(fishcount[i]/vol[i]), round(d_hat), log = TRUE)*fishcount[i]*-1
    # normal?
    ll=liknorm(fishcount[i]/vol[i], d_hat, log = TRUE)*fishcount[i]*-1
    # key part here is to multiply the likelihood by the number of observations in the bin.
    # this way bins with more fish in them drive the fit
    like[i]=ll
  }
  return(sum(like))
}
