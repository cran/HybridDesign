#' Generate modified mTPI Design Decision Boundary
#'
#' Generate dose escalation and deescalation boundaries of modified Toxicity Probability Interval (mTPI) design
#' with overdose control.
#'
#' @param target target toxicity rate
#' @param ncohort the total number of cohorts
#' @param cohortsize the cohort size
#' @param eps1 modified Toxicity Probability Interval (mTPI) design parameter epsilon1. Default: 0.05
#' @param eps2 modified Toxicity Probability Interval (mTPI) design parameter epsilon2. Default: 0.05
#' @param a Beta prior shape parameter 1. Default: 1
#' @param b Beta prior shape parameter 2. Default: 1
#' @param cutoff.eli Posterior probability cutoff of eliminating dose due to unacceptable toxicity. Default: 0.95
#' @param tox.control indicator of whether to perform toxicity control. If \code{TRUE}, change "stay"
#'                    to "deescalation" if the posterior probability of DLT rate greater than
#'                    \code{target+eps2} is greater than the toxicity control cutoff \code{cut.tox}
#' @param cut.tox toxicity control cutoff. Default: 0.8
#' @param esc.control indicator of whether to perform escalation control. If \code{TRUE}, change decision
#'                    of "escalation" to "stay" if the posterior probability of DLT rate less than
#'                    \code{target-eps1} is greater than the escalation control cutoff \code{cut.esc}
#' @param cut.esc escalation control cutoff. Default: 0.5
#'
#' @return This function returns the table of escalation and deescalation boundaries.

#' @examples
#' get_boundary_mtpi(target=0.30, ncohort=10, cohortsize=3)
#'
#' @export

get_boundary_mtpi <- function(target, ncohort, cohortsize, eps1=0.05, eps2=0.05, a=1, b=1, cutoff.eli=0.95,
                              tox.control=FALSE, cut.tox=0.8, esc.control=FALSE, cut.esc=0.5) {

  ### simple error checking
  if(target<0.05) {stop("Error: the target is too low! \n");}
  if(target>0.6)  {stop("Error: the target is too high! \n");}

  sampsize = ncohort * cohortsize
  decision_table = matrix(NA,nrow=sampsize+1,ncol=sampsize)
  for (ntr in 1:sampsize)
  {
    elim=0
    for (ntox in 0:ntr)
    {
      #if (1-stats::pbeta(pT, ntox+a, ntr-ntox+b)>xi) {elim=1;break;}
      q1 <- (1-stats::pbeta(eps2+target, ntox+a, ntr-ntox+b))/(1-eps2-target)
      q2 <- (stats::pbeta(eps2+target, ntox+a, ntr-ntox+b) - stats::pbeta(target-eps1, ntox+a, ntr-ntox+b))/(eps2+eps1)
      q3 <- (stats::pbeta(target-eps1, ntox+a, ntr-ntox+b)/(target-eps1))
      if (q3==max(q1,q2,q3)) decision_table[ntox+1,ntr] <- ifelse(tox.control,ifelse(1-stats::pbeta(eps2+target, ntox+a, ntr-ntox+b)>cut.tox,"S","E"),"E")
      if (q2==max(q1,q2,q3)) {
        temp1 <- ifelse(tox.control,ifelse(1-stats::pbeta(eps2+target, ntox+a, ntr-ntox+b)>cut.tox,"D","S"),"S")
        temp2 <- ifelse(esc.control,ifelse(stats::pbeta(target-eps1, ntox+a, ntr-ntox+b)>cut.esc,"E","S"),"S")
        if(temp1=="D") decision_table[ntox+1,ntr] <- "D"
        else if(temp2=="E") decision_table[ntox+1,ntr] <- "E"
        else decision_table[ntox+1,ntr] <- "S"
      }
      if (q1==max(q1,q2,q3)) decision_table[ntox+1,ntr] <- "D"
      if (1-stats::pbeta(target, ntox+a, ntr-ntox+b)>cutoff.eli) {elim=1;break;}
    }
    if (elim==1) decision_table[(ntox+1):(ntr+1),ntr] <- rep("DU",ntr-ntox+1)

  }
  colnames(decision_table) <- 1:sampsize
  rownames(decision_table) <- 0:sampsize

  boundary = matrix(NA, nrow=4, ncol=sampsize)
  boundary[1,] = 1:sampsize
  for(i in 1:sampsize)
  {
    boundary[2,i] = max(which(decision_table[,i]=="E"))-1
    #boundary[3,i] = min(which(decision_table[,i]=="D"))-1
    if (length(which(decision_table[,i]=="D"))) {boundary[3,i] = min(which(decision_table[,i]=="D"))-1}
    else if (length(which(decision_table[,i]=="DU"))) {boundary[3,i] = min(which(decision_table[,i]=="DU"))-1}
    if (length(which(decision_table[,i]=="DU"))) {boundary[4,i] = min(which(decision_table[,i]=="DU"))-1}
  }
  colnames(boundary) <- c(rep("", sampsize))
  rownames(boundary) <- c("Number of patients treated", "Escalate if # of DLT <=", "Deescalate if # of DLT >=", "Eliminate if # of DLT >=")


  out=list()
  if(cohortsize>1){
    out=list(
      boundary_tab=boundary[, (1:ncohort)*cohortsize],
      full_boundary_tab=boundary,decision_table=decision_table)
  }else
    out=list(full_boundary_tab=boundary[, (1:ncohort)*cohortsize],decision_table=decision_table)

  class(out) <- "mtpi"

  return(out);
}
