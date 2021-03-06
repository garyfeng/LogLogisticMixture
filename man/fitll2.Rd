% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/fitll2.r
\name{fitll2}
\alias{fitll2}
\title{MLE of the log-logistic mixture model}
\usage{
fitll2(data, strpara = "skip_mle", maxtime = 1000, ylim = 0.01,
  maintitle = "", plot_pdf = TRUE, plot_haz = TRUE, do_mle = TRUE)
}
\arguments{
\item{data}{A a vector or a list of vectors of response time, in milliseconds}

\item{strpara}{A string of parameters, e.g.
"?150/180/130, ?8/12/6, =1, ?6/7/2, ?0.7|0.01|0.98, ?100|150|50, ?0.05|0.15|0.01"
For details, see ____}

\item{maxtime}{Optional parameter for the maximal RT for the plot (default = 1000ms).
For plotting only. No data is ever censored.}

\item{ylim}{Optional parameter for the max density (default = 0.010) when using milliseconds}

\item{maintitle}{Optional title for the plot}

\item{plot_pdf}{Optional; if TRUE then plot the probability density function}

\item{plot_haz}{Optional; if TRUE the plot the hazard rate, which is on the same scale as the pdf}

\item{do-mle:}{Optional, TRUE to perform MLE, FALSE to return only empirical pdf and hazard rate, etc.}
}
\value{
list(fit=fit, stat=stat)
}
\description{
MLE of the log-logistic mixture model
}

