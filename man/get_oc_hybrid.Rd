% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_oc_hybrid.R
\name{get_oc_hybrid}
\alias{get_oc_hybrid}
\title{Generate operating characteristics for single-agent dose-finding studies using the Hybrid design}
\usage{
get_oc_hybrid(trueMTD, trueDLTvec, dose, target, ncohort, cohortsize,
                     eps1=0.05, eps2=0.05,a=1, b=1, cutoff.eli=0.95,
                     tox.control=TRUE, cut.tox=0.8, esc.control=FALSE, cut.esc=0.5,
                     ntrial, seednum=10000)
}
\arguments{
\item{trueMTD}{the dosage of true maximum tolerated dose (MTD)}

\item{trueDLTvec}{a vector of true dose-limiting toxicity (DLT) rates at each dose level}

\item{dose}{a vector containing the numerical dosage of each dose level}

\item{target}{target toxicity rate}

\item{ncohort}{the total number of cohorts}

\item{cohortsize}{the cohort size}

\item{eps1}{mTPI design parameter epsilon1. Default: 0.05}

\item{eps2}{mTPI design parameter epsilon2. Default: 0.05}

\item{a}{Beta prior shape parameter 1. Default: 1}

\item{b}{Beta prior shape parameter 2. Default: 1}

\item{cutoff.eli}{Posterior probability cutoff of eliminating dose due to unacceptable toxicity. Default: 0.95}

\item{tox.control}{indicator of whether to perform toxicity control. If \code{TRUE}, change "stay"
to "deescalation" if the posterior probability of DLT rate greater than
\code{target+eps2} is greater than the toxicity control cutoff \code{cut.tox}}

\item{cut.tox}{toxicity control cutoff. Default: 0.8}

\item{esc.control}{indicator of whether to perform escalation control. If \code{TRUE}, change decision
of "escalation" to "stay" if the posterior probability of DLT rate less than
\code{target-eps1} is greater than the escalation control cutoff \code{cut.esc}}

\item{cut.esc}{escalation control cutoff. Default: 0.5}

\item{ntrial}{the total number of trials to be simulated}

\item{seednum}{the random seed for simulation}
}
\value{
This function returns the operating characteristics of the Hybrid design as a list,
       including:
       (1) Percentage of correct selection of the true MTD in all simulated trials,
       (2) Percentage of selecting a dose above MTD in all simulated trials,
       (3) Percentage of selecting a dose below MTD in all simulated trials,
       (4) Average number of patients treated at MTD in all simulated trials.
}
\description{
Obtain the operating characteristics of the Hybrid design for single-agent dose-finding studies by simulation.
}
\examples{
\donttest{
get_oc_hybrid(trueMTD=12, trueDLTvec=c(0.15,0.20,0.25,0.30,0.35), dose=c(3, 6, 12, 18, 24),
              target=0.25, ncohort=10, cohortsize=3, eps1=0.05, eps2=0.05, a=1, b=1,
              cutoff.eli=0.95, tox.control=TRUE, cut.tox=0.8, esc.control=FALSE, cut.esc=0.5,
              ntrial=100, seednum=10000)
}

}
