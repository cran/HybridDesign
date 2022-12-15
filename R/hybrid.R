#' @import testit
#' @import ResourceSelection
NULL

#' Implement Hybrid design with real data
#'
#' Obtain decision for the next dose level to be tested given current trial data.
#'
#' @usage hybrid(dose, nDLT, npts, currdose, nextdose=0, target, ncohort, cohortsize,
#'               eps1=0.05, eps2=0.05, a=1, b=1, cutoff.eli=0.95, tox.control=TRUE,
#'               cut.tox=0.8, esc.control=FALSE, cut.esc=0.5, regrule)
#' @param dose a vector containing the numerical dosage of each dose level
#' @param nDLT a vector containing the number of patients who experienced dose-limiting
#'              toxicity at each dose level
#' @param npts a vector containing the number of patients at each dose level
#' @param currdose the dosage at the current dose level
#' @param nextdose the dosage of next higher dose level; could be an intermediate dose
#' @param target the target toxicity rate
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
#' @param regrule indicator of whether to apply additional overdose control rule

#' @return This function returns the decision of implementing the Hybrid design with real trial data
#'         as a list, including:
#'         (1) dose transition boundaries of modified mTPI design,
#'         (2) decision table of modified mTPI design,
#'         (3) the decision given current data,
#'         (4) the summary table of tested dose levels

#' @examples
#' hybrid(dose=c(2,4,8,16,22,28,40), nDLT=c(0,0,0,0,1,0,2), npts=c(3,3,4,6,9,5,16), currdose=40,
#'        nextdose=54, target=0.3, ncohort=10, cohortsize=3, eps1=0.05, eps2=0.05, a=1, b=1,
#'        cutoff.eli=0.95, tox.control=TRUE, cut.tox=0.8, regrule=TRUE)
#'
#' @export

hybrid <- function(dose, nDLT, npts, currdose, nextdose=0, target, ncohort, cohortsize, eps1=0.05, eps2=0.05, a=1, b=1, cutoff.eli=0.95,
                   tox.control=TRUE, cut.tox=0.8, esc.control=FALSE, cut.esc=0.5, regrule=TRUE){

  mTPIout <- get_boundary_mtpi(target, ncohort, cohortsize, eps1, eps2, a, b, cutoff.eli, tox.control, cut.tox, esc.control, cut.esc)

  mTPIbound <- mTPIout$full_boundary_tab
  mTPIbound[4,-which(mTPIbound[4,]<10000)] <- 9999
  decision_table <-  mTPIout$decision_table

  currDLT =nDLT[which(dose==currdose)]

  if (currDLT >= mTPIbound[4,npts[which(dose==currdose)]]) {decision1 <- "Eliminate"}
  else if (currDLT >= mTPIbound[3,npts[which(dose==currdose)]] & currDLT < mTPIbound[4,npts[which(dose==currdose)]]) {decision1 <- "Deescalate"}
  else if (currDLT > mTPIbound[2,npts[which(dose==currdose)]] & currDLT < mTPIbound[3,npts[which(dose==currdose)]]){decision1 <- "Stay"}
  else if (currDLT <= mTPIbound[2,npts[which(dose==currdose)]]){decision1 <- "Escalate"}

  mTPIbound[4, -which(mTPIbound[4,]<9999)] <- NA


  ## Logstic regression model

  nlength <- length(dose)
  ptox <- NULL
  dosevec <- NULL

  for (i in 1:nlength){
    ptox <- append(ptox, c(rep(0, (npts[i]-nDLT[i])), rep(1, nDLT[i])))
    dosevec <- append(dosevec, rep(dose[i], npts[i]))
  }

  warnmsg <- 0
  BLRMrule <- 1
  pval <- 0

  model <- stats::glm(ptox ~ dosevec, family = "binomial")
  warnmsg <- has_warning(model <- stats::glm(ptox ~ dosevec, family = "binomial"))

  if (is.na(model$coeff[2])==FALSE &  is.na(model$coeff[1])==FALSE){
    if (abs(model$coeff[2])> 0.0001 & abs(model$coeff[1])> 0.0001){
      pval <- hoslem.test(ptox, stats::fitted(model))$p.value
    }
  }

  # Use mTPI when the logistic model does not fit well

  if (is.na(model$coeff[2])==TRUE | warnmsg==1 | is.na(pval)==TRUE){BLRMrule <- 0}
  else if (is.na(model$coeff[2])==TRUE | warnmsg==1 | pval < 0.05){BLRMrule <- 0}
  else if (is.na(model$coeff[2])==FALSE & model$coeff[2]<0){BLRMrule <- 0}

  if (!regrule){BLRMrule <- 0}

  if (BLRMrule==1){

    ## Estimate toxicity for current, next, next intermediate, prior, prior intermediate dose

    action <- c(rep("Not Selected", 5),rep(NA,length(nextdose)-1))

    currptox <- stats::predict(model, data.frame(dosevec = currdose), type = "response")
    currvec <- round(cbind(currdose, currptox),3)

    nextptox <- stats::predict(model, data.frame(dosevec = nextdose), type = "response")
    nextvec <- round(cbind(nextdose, nextptox),3)

    nextdose0 = (currdose + nextdose[1])/2
    nextptox0 = stats::predict(model, data.frame(dosevec = nextdose0), type = "response")

    intervec <- round(cbind(nextdose0, nextptox0),3)

    priordose = dose[max(which(dose==currdose)-1,1)]
    priorptox = stats::predict(model, data.frame(dosevec = priordose), type = "response")
    priorvec = round(cbind(priordose, priorptox),3)

    priordose0 = (currdose + priordose)/2
    priorptox0 = stats::predict(model, data.frame(dosevec = priordose0), type = "response")
    priorvec0 = round(cbind(priordose0, priorptox0),3)

    nextptox11 <- stats::predict(model, data.frame(dosevec = 2*currdose), type = "response")

    if ((min(nextptox[1],nextptox11) > 0.95) & priorptox < 0.05 ){BLRMrule <- 0}


    ## mTPI suggests escalate
    if (decision1 == "Escalate"){
      if (currptox <= target & nextptox[1] <= target){
        decision2 <- decision1
        action[5] <- "Selected"
        priorvec0 <- cbind("NA","NA")
        intervec <- cbind("NA","NA")
      }
      else if (currptox <= target & nextptox[1] > target){
        if (nextptox0 > target){
          decision2 <- "Stay"
          action[1] <- "Selected"
          priorvec0 <- cbind("NA","NA")
        }
        else {
          decision2 <- decision1
          action[4] <- "Selected"
          priorvec0 <- cbind("NA","NA")
        }
      }
      else if (currptox > target){
        decision2 <- "Stay"
        action[1] <- "Selected"
        priorvec0 <- cbind("NA","NA")
        intervec <- cbind("NA","NA")
      }
    }

    ## mTPI suggests stay
    else if (decision1 == "Stay"){
      if (currptox <= target){
        decision2 <- decision1
        action[1] <- "Selected"
        priorvec0 <- cbind("NA","NA")
        intervec <- cbind("NA","NA")
      }
      else if (currptox > target){
        if (priorptox0 <= target){
          action[2] <- "Selected"
          intervec <- cbind("NA","NA")
        }
        else {
          action[3] <- "Selected"
          intervec <- cbind("NA","NA")
        }
        decision2 <- "Deescalate"

      }
    }

    ## mPTI suggests deescalate
    else if (decision1  == "Deescalate" || decision1 == "Eliminate"){
      if (priorptox0 <= target){
        action[2] <- "Selected"
        intervec <- cbind("NA","NA")
      }
      else {
        action[3] <- "Selected"
        intervec <- cbind("NA","NA")
      }
      decision2 <- decision1
    }

    ## Output table
    dosetab <- rbind(currvec, priorvec0, priorvec, intervec, nextvec)
    nextdoselb <- rep("Next Dose",nrow(nextvec))
    doselevel <- c("Current Dose","Prior Intermediate Dose", "Prior Dose", "Next Intermediate Dose", nextdoselb)
    dosetab <- cbind(doselevel, dosetab, action)
    colnames(dosetab) <- c("Dose Level","Dose", "Estimated Toxicity","Action")
    rownames(dosetab) <- seq(1:nrow(dosetab))
  }

  if (BLRMrule == 0){
    decision2 <- decision1
  }

  out=list()
  if (BLRMrule==0){
    out=list(
      mTPIboundary = mTPIbound,
      decision_table = decision_table,
      decision = decision2
    )
  }

  else if (BLRMrule==1){
    out=list(
      mTPIboundary = mTPIbound,
      decision_table = decision_table,
      decision = decision2,
      dose_toxicity_and_action_table = dosetab
    )
  }
  return(out);
}
