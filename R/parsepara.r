##############
# Parse parameters for the MLE estimation routine
# Gary Feng, copyleft, Sept 2010, Potsdam, Germany
##############
# Format of the input parameter specification:
#    e.g., str<-"2, ?0.05,3, ?200//400||100, =1"
#  Any numeric values --> fixed parameters, will not generate parameter list
#  ?## --> free parameter, ## is the initial value, upper and lower limites set to Inf, -Inf
#  ?##1/##2/##3 --> bounded parameter, ##1 is the init value, ##2 and ##3 will be the upper/lower bounds, auto re-ordered
#  =# --> tied parameter, pointing to parameter number #. 
#  @counter*2+180 --> a formula: CAUTION, make sure you know what you are doing.
#		the string after @ will be treated as a formula and will be evaluated at the run time, 
#		make sure it is a valid expression. So far the only "variable" it can use is "counter", which is an internal
#		variable in fitll2() that keeps track of which subset (list) of the data is currently being fitted (starting 1)
##############
# Returns
#  a list of:
#   -- input: original string
#   -- paralist: list of parameters, with free/tied parameters replaced with "para[#]"; 
#		Note the # of items is the same as in the input; i.e., fixed parameters are in here, in the same order. 
#   -- initpara/upperpara/lowerpara: bounds and init paras, with # equal to the # of free parameters. 
##############

parsepara <-function (str) {
	res<-sapply(strsplit(str, "[ ]*,[ ]*")[[1]], as.numeric); #names(res); res
	initpara<-c(); upperpara<-c(); lowerpara<-c(); paralist<-c(); c<-1;
	for (n in names(res)) {
		#cat(n, "\n");
		if(!is.na(res[n])) {
			paralist<-c(paralist, res[n]);#cat("paralist= ",paralist, "\n");
		} else if(substring(n, 1,1)=='?') {
			dump<- substring(n, 2); # get the rest of the strings
			#cat("dump= ",dump, "\n");
			args<- sapply(strsplit(dump, "[/|\\]+")[[1]], as.numeric)
			#cat("args= ",args, "\n");
			if(any(is.na(args))) cat("Error! Some parameters are not numeric: ", args)
			if(length(args)==3) {
				initpara<-c(initpara, args[1]); #cat("initpara= ",initpara, "\n");
				upperpara<-c(upperpara,max(args[2],args[3])); #cat("upperpara= ",upperpara, "\n");;
				lowerpara<-c(lowerpara,min(args[2],args[3])); #cat("lowerpara= ",lowerpara, "\n");;
			} else if(length(args)==1) {
				initpara<-c(initpara, args[1]); #cat("initpara= ",initpara, "\n");
				upperpara<-c(upperpara,Inf); #cat("upperpara= ",upperpara, "\n");;
				lowerpara<-c(lowerpara,-Inf); #cat("lowerpara= ",lowerpara, "\n");;
			} else {
				cat("Error! Format= '?initpara/upper/lower': ", args);
				return()
			}
			paralist<-c(paralist, paste("par[",c,"]", sep='')); #cat("paralist= ",paralist, "\n");
			c<-c+1;
		} else if(substring(n, 1,1)=='=') {
			dump<- substring(n, 2); # get the rest of the strings
			args<- as.numeric(dump);
			if (is.na(args))  {cat("Error! format= 'tieto3'. ", args); return();}
			#if (!is.na(args)) {cat("parameter tied to ", args, "\n");}
			paralist<-c(paralist, paste("par[",args,"]", sep=''));#cat("paralist= ",paralist, "\n");
		} else if(substring(n, 1,1)=='@') {	# formula
			dump<- substring(n, 2); # get the rest of the strings
			# if it's a formula, we do nothing now, just stick it into theparalist
			paralist<-c(paralist, dump);#cat("paralist= ",paralist, "\n");
		} else {
			# the above should handle everything; we got here by error
			cat("Error! Unrecognized input\n");
			return();
		}
	}
	out<-list(input=str, paralist=paralist, initpara=initpara, upper=upperpara, lower=lowerpara);
	return(out)
}
# test
# str<-"2, ?0.05,3, ?200//400||100, =1"
# parsepara(str)