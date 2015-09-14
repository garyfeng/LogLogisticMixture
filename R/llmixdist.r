# Gary Feng, 
# Sept 2010, Potsdam, Germany
# Maxmum Likelihood Estimation of the Log-logistic Mixture model 
# 
# For an application in eye movement research, see Feng (2012)
# This is a mixture of 2 log-logistic distributions, with a shift parameter and a mixing prob. parameter
#	Note: The MLE seems to be pretty sensitive to the initial values
#		Testing with 50% of variation in the initial parameters and different sameple sizes.
#	Also, the LL for the "true" model may not be the maximum compared to other possible models,
#		due to random errors. I often find the estimated para are actually better than the "true" para. 
#		An issue of sample size?
# to-do:
#	---1). Tieing the location parameters
#	---2). adding the express saccade prob mass estimation
#	3). use LeastSquare or Chi-Square methods to estimate the init parameters 
#	---4). detect fit failure and automatically re-estimate using different init-para
#	---5). some ways to optimize globally, using other optim algorithms? Grid search?
#	-X-6). bayesian estimation?? <-- way to slow and not practical for large datasets.
#	---7). Add AIC/BIC

###############################
# get pdf & haz
###############################
require(muhaz)

#' Calculate the pdf, haz, and other info for a given rt data
#' 
#' @param fixdur A numeric vector of RT, in milliseconds; hopefully no NAs. 
#' @param binwidth The bin width for calculating the pdf and hazard rates; default to 25 msec
#' @return a data frame containing the bins, pdf, haz, counts, etc. 
#' 
#' @export
getDistr <-function(fixdur, binwidth=2) {
  fixdur<-round(fixdur/binwidth)*binwidth
  status<-rep(1,length(fixdur))
  out1<-kphaz.fit(fixdur,status)
  
  # get the density fun
  fixhist <-hist(fixdur, c(0, out1$time, max(fixdur)), plot=F)
  
  # the lengths don't match.
  n<-length(fixhist$density);
  #out1$time2 <-fixhist$breaks[2:n]
  out1$pdf <-fixhist$density[2:n]
  out1$counts <-fixhist$counts[2:n]
  
  time2<-seq(0, max(fixdur), binwidth);
  fixhist2 <-hist(fixdur, breaks=time2, plot=F)
  out1$time2<-fixhist2$breaks[2:length(fixhist2$breaks)-1];
  out1$pdf2<-fixhist2$density;
  out1$counts2 <-fixhist2$counts

  return (out1)
}

###############################
# global vars
###############################
lx.default<-80; sx.default<-4; # for the Px component

###############################
# Log-logistic functions
###############################
# functions for Log-logistic density, cumulative prob., and random number generation

#' A random log-logistic variable generator, by reversing the CDF
#' @param n The size of the return vector
#' @param location The location parameter, default = 0, which always returns 0
#' @param scale The scale parameter, defaulting to 1
#' @return A vector of random values following the specified log-logistic distribution
#' @export
rllogis <-function (n, location = 0, scale = 1) 
{
  u <-runif(n)
	location * (u/(1-u))^(1/scale)
}

#' The log-logistic (cumulative) probility function
#' @param q A vector of the x values to evaluate the probability function
#' @param location The location parameter of the log-logistic distribution
#' @param scale The scale parameter of the log-logistic distribution
#' @param d The translational or shift parameter, i.e., (q-d) follows the standard log-logistic distribution
#' @param lower.tail Whether to output the probability or survival function
#' @param log.p Whether to output the log of the probability (or survival) function
#' @export
pllogis <- function (q, location, scale, d, lower.tail = TRUE, log.p = FALSE) 
{
  #if (q<0) return (0)
	#if (q<d) return (0)
	p <- 1/(1+((q-d)/location)^(-scale))
	# make sure it's valid prob.
	p[q<0] <-0;
	p[q<d] <-0;
	if (!lower.tail) p <- 1-p;
	if (log.p) p <-log(p);
	return (p)
}

#' The probability density function of the log-logistic distribution
#' @param x A vector of the x values to evaluate the probability function
#' @param location The location parameter of the log-logistic distribution
#' @param scale The scale parameter of the log-logistic distribution
#' @param d The translational or shift parameter, i.e., (x-d) follows the standard log-logistic distribution
#' @export
dllogis <-function (x, location, scale, d) {
    if (any(x <= 0)) 
        stop("x elements must be strictly positive")
    x[x<=d] <- d+.Machine$double.xmin
		# .Machine.double.xmin
	a=scale;
	b=1/location;
	p<-a * b *(b * (x-d))^(a-1) * (1+(b * (x-d))^a)^-2;
	return (p);
}

#' A function that returns the probability density of a 2-log-logistic-mixture model
#' 
#' @param x A vector of the x values to evaluate the probability function
#' @param location1 The location parameter of the first log-logistic distribution
#' @param scale1 The scale parameter of the first log-logistic distribution
#' @param location2 The location parameter of the second log-logistic distribution
#' @param scale2 The scale parameter of the second log-logistic distribution
#' @param mixp The mixture probability of component 1
#' @param d The translational or shift parameter of the second log-logistic distribution
#' @return A vector of the probability densities for the 2ll mixture
#' @export
# had to make sure it doesn't return 0
dll2 <- function (x, location1, scale1, location2, scale2, mixp, d)
{
  if (any(x <= 0)) 
        stop("x elements must be strictly positive")
	p<- mixp * dllogis(x, location1, scale1, 0) + (1-mixp) * dllogis(x, location2, scale2, d)
	#mixp * dlogis(log(x), location1, scale1)/x 
	#+ (1-mixp) * dlogis(log(x), location2, scale2)/(x)
	
	# .Machine$double.xmin
	p[p<=0]<-.Machine$double.xmin
	p[p>1] <-.Machine$double.xmin
	return(p)
}

#' A function that returns the density for 2-log-logistic-mixture + eXpress component model
#' 
#' The "eXpress" component often appears in the data may be a result of noise or corrective 
#'   responses that have little to do with the present task.
#' Note the mixture probabilities are set up in the following way:
#'   px + (1-px)(mixp + 1-mixp) === 100%
#'   
#' @param x A vector of the x values to evaluate the probability function
#' @param location1 The location parameter of the first log-logistic distribution
#' @param scale1 The scale parameter of the first log-logistic distribution
#' @param location2 The location parameter of the second log-logistic distribution
#' @param scale2 The scale parameter of the second log-logistic distribution
#' @param mixp The mixture probability of component 1; (1-mixp) is the proportion of component 2
#' @param d The translational or shift parameter of the second log-logistic distribution
#' @param px The mixing probability of the eXpress component; (1-px) are the proportion of RTs for the 2ll mixture
#' @param lx The location parameter for the eXpress component, which can be set by the global lx.default
#' @param sx The scale parameter for the eXpress component, which can be set by the global sx.default
#' @return A vector of the probability densities for the ll2x mixture
#' @export
# Px is free, location=lx=a=80, shape=sx=b=6 by default; SD(80/6)=26.84; SD(80/8)=19.22;
# px_mean =a * pi/b/sin(pi/b)
# px_sd = sqrt(a^2*pi/b*(2/sin(2*pi/b) - pi/b*(sin(pi/b))^-2))
dll2x <- function (x, location1, scale1, location2, scale2, mixp, d, px, lx=lx.default, sx=sx.default)
{
    if (any(x <= 0)) 
        stop("x elements must be strictly positive")
	# p<- mixp * dllogis(x, location1, scale1, 0) + (1-mixp) * dllogis(x, location2, scale2, d)
	p<- px * dllogis(x, lx, sx, 0) + (1-px)*(mixp * dllogis(x, location1, scale1, 0) + (1-mixp) * dllogis(x, location2, scale2, d));
	#mixp * dlogis(log(x), location1, scale1)/x 
	#+ (1-mixp) * dlogis(log(x), location2, scale2)/(x)
	
	# .Machine$double.xmin
	p[p<=0]<-.Machine$double.xmin
	p[p>1] <-.Machine$double.xmin
	return(p)
	
}

#' A function that returns the (cumulative) probability for 2-log-logistic-mixture
#' 
#' @param x A vector of the x values to evaluate the probability function
#' @param location1 The location parameter of the first log-logistic distribution
#' @param scale1 The scale parameter of the first log-logistic distribution
#' @param location2 The location parameter of the second log-logistic distribution
#' @param scale2 The scale parameter of the second log-logistic distribution
#' @param mixp The mixture probability of component 1
#' @param d The translational or shift parameter of the second log-logistic distribution
#' @return A vector of the probability densities for the 2ll mixture
#' @export

pll2 <- function (x, location1, scale1, location2, scale2, mixp, d)
{
    if (any(x < 0)) 
        stop("x elements must be strictly positive")
	p<- mixp * pllogis(x, location1, scale1, 0) + (1-mixp) * pllogis(x, location2, scale2, d)
	return(p)
	
}

#' A function that returns the (cumulative) probability for 2-log-logistic-mixture + eXpress component model
#' 
#' The "eXpress" component often appears in the data may be a result of noise or corrective 
#'   responses that have little to do with the present task.
#' Note the mixture probabilities are set up in the following way:
#'   px + (1-px)(mixp + 1-mixp) === 100%
#'   
#' @param x A vector of the x values to evaluate the probability function
#' @param location1 The location parameter of the first log-logistic distribution
#' @param scale1 The scale parameter of the first log-logistic distribution
#' @param location2 The location parameter of the second log-logistic distribution
#' @param scale2 The scale parameter of the second log-logistic distribution
#' @param mixp The mixture probability of component 1; (1-mixp) is the proportion of component 2
#' @param d The translational or shift parameter of the second log-logistic distribution
#' @param px The mixing probability of the eXpress component; (1-px) are the proportion of RTs for the 2ll mixture
#' @param lx The location parameter for the eXpress component, which can be set by the global lx.default
#' @param sx The scale parameter for the eXpress component, which can be set by the global sx.default
#' @return A vector of the probability densities for the ll2x mixture
#' @export

# probability for 2-log-logistic-mixture + eXpress saccade
# Px is free, location=lx=a=80, shape=sx=b=6 by default; SD(80/6)=26.84; SD(80/8)=19.22;
# px_mean =a * pi/b/sin(pi/b)
# px_sd = sqrt(a^2*pi/b*(2/sin(2*pi/b) - pi/b*(sin(pi/b))^-2))
pll2x <- function (x, location1, scale1, location2, scale2, mixp, d, px, lx=lx.default, sx=sx.default)
{
    if (any(x <= 0)) 
        stop("x elements must be strictly positive")
	# p<- mixp * dllogis(x, location1, scale1, 0) + (1-mixp) * dllogis(x, location2, scale2, d)
	p<- px * pllogis(x, lx, sx, 0) + (1-px)*(mixp * pllogis(x, location1, scale1, 0) + (1-mixp) * pllogis(x, location2, scale2, d));
	#mixp * dlogis(log(x), location1, scale1)/x 
	#+ (1-mixp) * dlogis(log(x), location2, scale2)/(x)
	
	# .Machine$double.xmin
	p[p<=0]<-.Machine$double.xmin
	p[p>1] <-.Machine$double.xmin
	return(p)
	
}


#' A function that returns the survival function for 2-log-logistic-mixture
#' 
#' @param x A vector of the x values to evaluate the probability function
#' @param location1 The location parameter of the first log-logistic distribution
#' @param scale1 The scale parameter of the first log-logistic distribution
#' @param location2 The location parameter of the second log-logistic distribution
#' @param scale2 The scale parameter of the second log-logistic distribution
#' @param mixp The mixture probability of component 1
#' @param d The translational or shift parameter of the second log-logistic distribution
#' @return A vector of the probability densities for the 2ll mixture
#' @export
sll2 <- function (x, location1, scale1, location2, scale2, mixp, d)
{
	return(1-pll2(x, location1, scale1, location2, scale2, mixp, d))
}

#' A function that returns the survival function for 2-log-logistic-mixture with eXpress saccade component
#' 
#' The "eXpress" component often appears in the data may be a result of noise or corrective 
#'   responses that have little to do with the present task.
#' Note the mixture probabilities are set up in the following way:
#'   px + (1-px)(mixp + 1-mixp) === 100%
#'   
#' @param x A vector of the x values to evaluate the probability function
#' @param location1 The location parameter of the first log-logistic distribution
#' @param scale1 The scale parameter of the first log-logistic distribution
#' @param location2 The location parameter of the second log-logistic distribution
#' @param scale2 The scale parameter of the second log-logistic distribution
#' @param mixp The mixture probability of component 1; (1-mixp) is the proportion of component 2
#' @param d The translational or shift parameter of the second log-logistic distribution
#' @param px The mixing probability of the eXpress component; (1-px) are the proportion of RTs for the 2ll mixture
#' @param lx The location parameter for the eXpress component, which can be set by the global lx.default
#' @param sx The scale parameter for the eXpress component, which can be set by the global sx.default
#' @return A vector of the probability densities for the ll2x mixture
#' @export
# use default lx, sx parameters
sll2x <- function (x, location1, scale1, location2, scale2, mixp, d, px)
{
	return(1-pll2x(x, location1, scale1, location2, scale2, mixp, d, px))
}


#' A function that returns the hazard rate for 2-log-logistic-mixture
#' 
#' @param x A vector of the x values to evaluate the probability function
#' @param location1 The location parameter of the first log-logistic distribution
#' @param scale1 The scale parameter of the first log-logistic distribution
#' @param location2 The location parameter of the second log-logistic distribution
#' @param scale2 The scale parameter of the second log-logistic distribution
#' @param mixp The mixture probability of component 1
#' @param d The translational or shift parameter of the second log-logistic distribution
#' @return A vector of the probability densities for the 2ll mixture
#' @export
hll2 <- function (x, location1, scale1, location2, scale2, mixp, d)
{
	h<-dll2(x, location1, scale1, location2, scale2, mixp, d)/sll2(x, location1, scale1, location2, scale2, mixp, d)
	return(h)
}

#' A function that returns the hazard rate for 2-log-logistic-mixture with eXpress saccade component
#' 
#' The "eXpress" component often appears in the data may be a result of noise or corrective 
#'   responses that have little to do with the present task.
#' Note the mixture probabilities are set up in the following way:
#'   px + (1-px)(mixp + 1-mixp) === 100%
#'   
#' @param x A vector of the x values to evaluate the probability function
#' @param location1 The location parameter of the first log-logistic distribution
#' @param scale1 The scale parameter of the first log-logistic distribution
#' @param location2 The location parameter of the second log-logistic distribution
#' @param scale2 The scale parameter of the second log-logistic distribution
#' @param mixp The mixture probability of component 1; (1-mixp) is the proportion of component 2
#' @param d The translational or shift parameter of the second log-logistic distribution
#' @param px The mixing probability of the eXpress component; (1-px) are the proportion of RTs for the 2ll mixture
#' @param lx The location parameter for the eXpress component, which can be set by the global lx.default
#' @param sx The scale parameter for the eXpress component, which can be set by the global sx.default
#' @return A vector of the probability densities for the ll2x mixture
#' @export
hll2x <- function (x, location1, scale1, location2, scale2, mixp, d, px)
{
	h<-dll2x(x, location1, scale1, location2, scale2, mixp, d, px)/sll2x(x, location1, scale1, location2, scale2, mixp, d, px)
	return(h)
}
