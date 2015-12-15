
##############################
# MLE of log-logistic mixture from data
# Gary Feng, Copyleft
# Sept 2010
##############################
# source("ll2sim_3.r")
#source("R/llmixdist.r")	# functions for generating/evaluating loglogistic and their mixtures
#source("R/parsepara.r")	# parameter parser

###############################
# fit 2.5 component loglogistic mixtures
# 	that is, a semi-fixed eXpress saccade component, controlled by fixed parameters (global var) lx=80 and sx=3 
#	P component typically with location=180, shape =8 for reading
#	C component typically with location=180 (tied to that the of the P component), shape=4, and a delay of d=100
#	The full model estimates 7 parameters: l1/l2=locations, s1/s2=shapes, 
#		pp=probs of P =(1-prob. of C), px=prob of eXpress sacc, and d=delay of C component. 
#		Note the mixture model is specified as:
#			f(t) = Px fx(t) + (1-Px) (Pc fc(t) + Pp fp(t)), where fx(t) is the pdf of the eXpress saccades, which are fixed
#	  The order of parameters are: {l1, s1, l2, s2, pp, d, px}, 
#		e.g., the following specifies a model with fixed d=100, tied l1/l2, and the rest are free parameters
#		"?150/180/130, ?8/12/6, =1, ?6/7/2, ?0.7|0.01|0.98, 100, ?0.05|0.15|0.01"
# Input:
#	data, is a list of vectors; i.e., you put fixdur/srt of each condition as a vector (in msec), and
#		put all conditions into a list, and name the vectors in the list meaningfully. e.g. you can do:
#		data<-list(preview100=data_for_preview_100, control=data_for_control);
# 	strpara: optional parameter, a string that specifies free/fix/tied parameters with upper/lower bounds
#		for the format of strpara, see the head of "parsepara.r" and example code; if missing, then no MLE
#	And the rest are control optional paras that should be self-explanatory. 
# Returns:
#	A list with MLE fits for each condition in the data list, along with AIC, BIC, etc.
#################################

#' MLE of the log-logistic mixture model
#' 
#' @param data A a vector or a list of vectors of response time, in milliseconds
#' @param strpara A string of parameters, e.g.
#'   "?150/180/130, ?8/12/6, =1, ?6/7/2, ?0.7|0.01|0.98, ?100|150|50, ?0.05|0.15|0.01"
#'   For details, see ____
#' @param maxtime Optional parameter for the maximal RT for the plot (default = 1000ms). 
#'   For plotting only. No data is ever censored. 
#' @param ylim Optional parameter for the max density (default = 0.010) when using milliseconds
#' @param maintitle Optional title for the plot
#' @param plot_pdf Optional; if TRUE then plot the probability density function
#' @param plot_haz Optional; if TRUE the plot the hazard rate, which is on the same scale as the pdf
#' @param do-mle: Optional, TRUE to perform MLE, FALSE to return only empirical pdf and hazard rate, etc.
#' @return list(fit=fit, stat=stat)
#' @export

fitll2 <- function(data, 
                   strpara='skip_mle', maxtime=1000, ylim=0.010, 
                   maintitle="", plot_pdf=TRUE, plot_haz=TRUE, 
                   do_mle=TRUE) {
	if (!do_mle) strpara<-'skip_mle';
	if (missing(data)) {
		cat("Data should be list of one or more vectors of SRT/fixation duration\n");
		return()
	} else if (class(data)!=  "numeric" & class(data)!=  "list") {
		cat("Data should be a vector or a list of vectors of SRT/fixation duration\n");
		return()
	} else if (missing(strpara)) {
		cat("Please supply the correct parameters; see 'parsepara.r' for instructions\n");
		return()
	} else if (strpara=='skip_mle') {
		cat("\nNo MLE parameter specified. Will only show empirical data. \nTo perform MLE, please supply the correct parameters; see 'parsepara.r' for instructions\n");
		do_mle=FALSE;
	}
	# handle the case when data is a vector (single condition)
	if (class(data) == "list") {datalist<-data;} 
	else {datalist<-list(data=data);}
	
	# set up things
	fit <-list(); stat <-list();colors <-c('red', 'blue','black', 'green', 'brown', 'grey');
	if(do_mle){
		# now parse the LL parameters
		para<-parsepara(strpara); #parse the parameter, return a list of named lists, see parsepara.r
		if (is.null(para)) {
			cat("Please supply the correct parameters; see 'parsepara.r' for instructions\n");
			return()
		}
		strparalist<-paste(para$paralist, collapse=", "); cat("para=",strparalist,"\n");
		#llfun<-eval(parse(text = paste("function(par) {-sum(log(dll2x(data, ", strparalist, ")));}")));
		initpara<-para$initpara; upperpara<-para$upper; lowerpara<-para$lower;
		# now do the estimates
		# TO-DO: can use lapply() instead of the for loop
		counter<-1;
		for (i in names(datalist)) {
			data<-datalist[[i]]; data<-data[data>0]; data<-data[is.finite(data)];
			# optimize based on the (negative) log likelihood. 
			mleFit<- optim(fn = eval(parse(text = 
			                                 paste("function(par) {-sum(log(dll2x(data, ", strparalist, ")));}"))), 
			               par = initpara, method = "L-BFGS-B", upper=upperpara, lower=lowerpara, hessian = T)
			mleFit$initpara<-initpara;
			mleFit$npara <-length(initpara); # only non identical upper/lower limites
			mleFit$n <- length(data);
			# calculate the AIC and BIC
			mleFit$BIC <- log1p(mleFit$n)*mleFit$npara + 2*mleFit$value; #using log1p to prevent n=0; shouldn't be an issue
			mleFit$AIC <- 2*mleFit$npara + 2*mleFit$value;
			fit[[i]]  <-mleFit;
			# calc the pdf, hazard rate, etc.
			stat[[i]] <-getDistr(data, binwidth=20);
			counter<-counter+1;
		}
	} else {	# will only do empirical pdf/haz, fit will be an empty list
		cat('\nwill only do empirical pdf/haz\n')
		for (i in names(datalist)) {
			data<-datalist[[i]]; data<-data[data>0]; data<-data[is.finite(data)];
			stat[[i]] <-getDistr(data, binwidth=20);
			cat('getDistr ',i,'\n')
		}
	}
	
	# TO-DO: move the plotting function out
	if (plot_pdf) {
		# set up the pdf
		plot(0,0, type='n', xlim=c(0,min(maxtime, 1000)), ylim=c(0,ylim),frame.plot=FALSE, 
		     xlab="Fixation Duration", ylab="Probability Density");
		title(main = maintitle)
		c<-1; time2<-seq(10, min(maxtime, 1000), by=1);
		counter<-1
		for (i in names(stat)) {
			lines(stat[[i]]$time, stat[[i]]$pdf, col=colors[c], lty='dashed'); 
			if (do_mle) {
				#lines(time2-10, dll2x(time2, a, b1, a, b2, mixp, d, px), col=colors[c]); # shared L parameter
				lines(time2-10, eval(parse(text = paste("dll2x(time2, ", 
				                                        gsub("par", "fit[[i]]$par", strparalist), ")"))), 
				      col=colors[c]); # paralist replaced with fit[[i]]$par
				text(min(maxtime, 1000)*0.3, (ylim-c*0.0015), 
				     paste(c("parameters=", unlist(para$paralist)), collapse=", "), col=colors[c], pos=4, cex=0.8);
				text(min(maxtime, 1000)*0.3, (ylim-c*0.0015-0.0005), 
				     paste(c(i, unlist(round(fit[[i]]$par, 2))), collapse=", "), col=colors[c], pos=4, cex=0.8);
				text(min(maxtime, 1000)*0.3, (ylim-c*0.0015-0.0010), 
				     paste(c('n=', 'LL=', 'AIC=', 'BIC='), 
				           c(unlist(round(c(fit[[i]]$n, fit[[i]]$value, fit[[i]]$AIC, fit[[i]]$BIC), 2))), collapse=", "), 
				     col=colors[c], pos=4, cex=0.8);
			} else {
				text(min(maxtime, 1000)*0.3, (ylim-c*0.0015), i, col=colors[c], pos=4, cex=0.8);
			}
			c=c+1;counter<-counter+1;
		}
	}
	if (plot_haz) {
		c<-1;
		#set up the haz; missing the eXpress saccade component
		plot(0,0, type='n', xlim=c(0,min(maxtime, 1000)), ylim=c(0,ylim*2.5), frame.plot=FALSE, 
		     xlab="Fixation Duration", ylab="Hazard Rate");
		title(main = maintitle)
		counter<-1
		for (i in names(stat)) {
			lines(stat[[i]]$time, stat[[i]]$haz, col=colors[c], lty='dashed'); 
			if (do_mle) {
				lines(time2-10, eval(parse(text = paste("hll2x(time2, ", 
				                                        gsub("par", "fit[[i]]$par", strparalist), ")"))), 
				      col=colors[c]); # paralist replaced with fit[[i]]$par
			} else {
				text(min(maxtime, 1000)*0.3, (ylim*2.5-c*0.0015), i, col=colors[c], pos=4, cex=0.8);
			}
			c=c+1;counter<-counter+1;
		}
	
	}
	
	return (list(fit=fit, stat=stat))
}
