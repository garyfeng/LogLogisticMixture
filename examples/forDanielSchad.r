
source("ll2sim_3.r")	# MLE estimate and a whole load of junk
source("parsepara.r")	# look here for how to specify fixed, free, and tied parameters

# load data
ds<-load("gary.rda")
ds<-list(mindless=i1st$dur[i1st$condition=='mindless reading'], mindfull=i1st$dur[i1st$condition=='mindful reading'],control=i1st$dur[i1st$condition=='normal reading'])
# set the parameters for the express saccade, which will be fixed parameters
lx.default<-100; sx.default<-4; # for the Px component

# func to do MLE
fitll2_2 <- function(data, strpara, maxtime=1000, ylim=0.010) {
	if (missing(data)) {
		cat("Data should be a vector of SRT/fixation duration\n");
		return()
	} else if (class(data)!=  "numeric" & class(data)!=  "list") {
		cat("Data should be a vector or a list of vectors of SRT/fixation duration\n");
		return()
	} else if (missing(strpara)) {
		cat("Please supply the correct parameters; see 'parsepara.r' for instructions\n");
		return()
	}
	if (class(data) == "list") {datalist<-data;} 
	else {datalist<-list(data=data);}
	# now parse the LL parameters
	para<-parsepara(strpara); #parse the parameter, return a list of named lists, see parsepara.r
	if (is.null(para)) {
		cat("Please supply the correct parameters; see 'parsepara.r' for instructions\n");
		return()
	}
	strparalist<-paste(para$paralist, collapse=", "); cat("para=",strparalist,"\n");
	llfun<-eval(parse(text = paste("function(par) {-sum(log(dll2x(data, ", strparalist, ")));}")));
	initpara<-para$initpara; upperpara<-para$upper; lowerpara<-para$lower;
	# set up things
	fit <-list(); stats <-list();colors <-c('red', 'blue','black', 'green', 'brown', 'grey');
	# now do the estimates
	for (i in names(datalist)) {
		data<-datalist[[i]]; data<-data[data>0]; data<-data[is.finite(data)];
		mleFit<- optim(fn = eval(parse(text = paste("function(par) {-sum(log(dll2x(data, ", strparalist, ")));}"))), par = initpara, method = "L-BFGS-B", upper=upperpara, lower=lowerpara, hessian = T)
		mleFit$initpara<-initpara;
		mleFit$npara <-length(initpara); # only non identical upper/lower limites
		mleFit$n <- length(data);
		mleFit$BIC <- log1p(mleFit$n)*mleFit$npara + 2*mleFit$value; #using log1p to prevent n=0; shouldn't be an issue
		mleFit$AIC <- 2*mleFit$npara + 2*mleFit$value;
		fit[[i]]  <-mleFit;
		stat[[i]] <-getDistr(data, binwidth=20);
	}

	if (T) {
		# set up the pdf
		plot(0,0, type='n', xlim=c(0,min(maxtime, 1000)), ylim=c(0,ylim),frame.plot=FALSE, xlab="Fixation Duration", ylab="Probability Density");
		c<-1; time2<-seq(10, min(maxtime, 1000), by=1);
		for (i in names(fit)) {
			lines(stat[[i]]$time, stat[[i]]$pdf, col=colors[c], lty='dashed'); 
			#lines(time2-10, dll2x(time2, a, b1, a, b2, mixp, d, px), col=colors[c]); # shared L parameter
			lines(time2-10, eval(parse(text = paste("dll2x(time2, ", gsub("par", "fit[[i]]$par", strparalist), ")"))), col=colors[c]); # paralist replaced with fit[[i]]$par
			text(min(maxtime, 1000)*0.3, (ylim-c*0.0015), paste(c("parameters=", unlist(para$paralist)), collapse=", "), col=colors[c], pos=4, cex=0.8);
			text(min(maxtime, 1000)*0.3, (ylim-c*0.0015-0.0005), paste(c(i, unlist(round(fit[[i]]$par, 2))), collapse=", "), col=colors[c], pos=4, cex=0.8);
			text(min(maxtime, 1000)*0.3, (ylim-c*0.0015-0.0010), paste(c('n=', 'LL=', 'AIC=', 'BIC='), c(unlist(round(c(fit[[i]]$n, fit[[i]]$value, fit[[i]]$AIC, fit[[i]]$BIC), 2))), collapse=", "), col=colors[c], pos=4, cex=0.8);
			c=c+1;
		}
		c<-1;
		#set up the haz; missing the eXpress saccade component
		plot(0,0, type='n', xlim=c(0,min(maxtime, 1000)), ylim=c(0,ylim*2), frame.plot=FALSE, xlab="Fixation Duration", ylab="Hazard Rate");
		for (i in names(fit)) {
			lines(stat[[i]]$time, stat[[i]]$haz, col=colors[c], lty='dashed'); 
			a<-fit[[i]]$par[1]; b1<-fit[[i]]$par[2]; b2<-fit[[i]]$par[3]; 
			mixp<-fit[[i]]$par[4]; d<-fit[[i]]$par[5]; px<-fit[[i]]$par[6]; 
			#lines(time2-10, hll2x(time2, a, b1, a, b2, mixp, d, px), col=colors[c]); # No Px parameter
			lines(time2-10, eval(parse(text = paste("hll2x(time2, ", gsub("par", "fit[[i]]$par", strparalist), ")"))), col=colors[c]); # paralist replaced with fit[[i]]$par
			c=c+1
		}
	
	}
	
	return (mleFit)
}

# ready to play
dev.new(); par(mfrow=c(2,1));
full		<-fitll2_2(ds,  "?190, ?9, ?190, ?4/7/2, ?0.7|0.1|0.9, ?80, ?0.05|0.15|0.01", ylim=0.008)

dev.new(); par(mfrow=c(2,1));
restricted	<-fitll2_2(ds,  "190, 10.2, 190, ?4/7/2, 0.3, 80, 0.03", ylim=0.008)