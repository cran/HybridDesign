#'
#' Generate operating characteristics for single-agent dose-finding studies using the Hybrid design
#'
#' Obtain the operating characteristics of the Hybrid design for single-agent dose-finding studies by simulation.
#'
#' @usage get_oc_hybrid(trueMTD, trueDLTvec, dose, target, ncohort, cohortsize,
#'                      eps1=0.05, eps2=0.05,a=1, b=1, cutoff.eli=0.95,
#'                      tox.control=TRUE, cut.tox=0.8, esc.control=FALSE, cut.esc=0.5,
#'                      ntrial, seednum=10000)
#' @param trueMTD the dosage of true maximum tolerated dose (MTD)
#' @param trueDLTvec a vector of true dose-limiting toxicity (DLT) rates at each dose level
#' @param dose a vector containing the numerical dosage of each dose level
#' @param target target toxicity rate
#' @param ncohort the total number of cohorts
#' @param cohortsize the cohort size
#' @param eps1 mTPI design parameter epsilon1. Default: 0.05
#' @param eps2 mTPI design parameter epsilon2. Default: 0.05
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
#' @param ntrial the total number of trials to be simulated
#' @param seednum the random seed for simulation

#' @return This function returns the operating characteristics of the Hybrid design as a list,
#'        including:
#'        (1) Percentage of correct selection of the true MTD in all simulated trials,
#'        (2) Percentage of selecting a dose above MTD in all simulated trials,
#'        (3) Percentage of selecting a dose below MTD in all simulated trials,
#'        (4) Average number of patients treated at MTD in all simulated trials.
#'
#' @examples
#' \donttest{
#' get_oc_hybrid(trueMTD=12, trueDLTvec=c(0.15,0.20,0.25,0.30,0.35), dose=c(3, 6, 12, 18, 24),
#'               target=0.25, ncohort=10, cohortsize=3, eps1=0.05, eps2=0.05, a=1, b=1,
#'               cutoff.eli=0.95, tox.control=TRUE, cut.tox=0.8, esc.control=FALSE, cut.esc=0.5,
#'               ntrial=100, seednum=10000)
#' }
#'
#' @export

get_oc_hybrid <- function(trueMTD, trueDLTvec, dose, target, ncohort, cohortsize, eps1=0.05, eps2=0.05, a=1, b=1, cutoff.eli=0.95,
                          tox.control=TRUE, cut.tox=0.8, esc.control=FALSE, cut.esc=0.5, ntrial, seednum=10000){

  dosevec1 <- dose
  nmax <- ncohort*cohortsize
  dselect <- rep(0, ntrial)

  ncorrect <- 0;
  overdose <- 0;
  underdose <- 0;
  MTDpts <- rep(0, ntrial)

  if(target<0.05) {stop("Error: the target is too low! \n");}
  if(target>0.6)  {stop("Error: the target is too high! \n");}

  set.seed(seednum)

  mTPIbound <- get_boundary_mtpi(target, ncohort, cohortsize, eps1, eps2, a, b, cutoff.eli, tox.control, cut.tox, esc.control, cut.esc)$full_boundary_tab

  mTPIbound[4,-which(mTPIbound[4,]<10000)] <- 9999

  for (trial in 1:ntrial){

    dosevec <- dose
    trueDLT <- trueDLTvec
    dosenum <- length(trueDLT)
    DLTvec <- rep(0,dosenum)   ### Cumulative number of DLTs  per dose level
    ptsvec <- rep(0,dosenum)    ### Cumulative number of subjects per dose level
    ptsvec0 <- NULL
    DLTvec0 <- NULL
    dosevec0 <- NULL
    elimdose <- 0

    i=1;imax <- dosenum

    repeat{

      DLTvec[i] <- stats::rbinom(1,cohortsize,trueDLT[i])+DLTvec[i]
      ptsvec[i] <- ptsvec[i]+cohortsize

      ptsvec0 <- ptsvec[ptsvec!=0]
      DLTvec0 <- DLTvec[ptsvec!=0]
      dosevec0 <- dosevec[ptsvec!=0]


      {
        if (DLTvec[i] >= mTPIbound[4, ptsvec[i]]){
          imax <- i-1
          if (imax == 0){break}

        }
        j <- min((i+1),imax)
      }

      rule <- ifelse(trueDLT[min((i+1), imax)]>target, TRUE, FALSE)


      temp <- hybrid(dosevec0, DLTvec0, ptsvec0, dosevec[i], dosevec[j], target, ncohort, cohortsize,  eps1, eps2, a, b, cutoff.eli,
                     tox.control, regrule=rule, cut.tox, esc.control, cut.esc)

      nextstep <- temp$dose_toxicity_and_action_table[,4]=="Selected"
      action <- temp$dose_toxicity_and_action_table[,1][nextstep][1]


      decision <- temp$decision

      if (is.null(action)==TRUE){
        if (decision=="Escalate"){i <- min((i+1), imax)}
        else if (decision=="Stay"){i <- i}
        else if (decision=="Deescalate"){i <- max((i-1), 1)}
        else if (decision=="Eliminate"){
          imax <- i-1
          elimdose <- dosevec[i]
          if (imax == 0){break}
          else {i <- i-1}
        }
      }

      else if (is.null(action)==FALSE){

        if (action=="Prior Dose"||action=="Prior Intermediate Dose" ){
          i <- max((i-1), 1)
          if (decision=="Eliminate"){
            imax <- i-1
            elimdose <- dosevec[i]
            if (imax == 0){break}
          }
        }

        else if (action == "Current Dose"){i <- i}

        else if (action == "Next Intermediate Dose"){
          lower <- dosevec1[dosevec1<as.numeric(temp$dose_toxicity_and_action_table[,2][temp$dose_toxicity_and_action_table[,1]==action][1])][which.max(dosevec1[dosevec1<as.numeric(temp$dose_toxicity_and_action_table[,2][temp$dose_toxicity_and_action_table[,1]==action][1])]-as.numeric(temp$dose_toxicity_and_action_table[,2][temp$dose_toxicity_and_action_table[,1]==action][1]))]
          upper <- dosevec1[dosevec1>as.numeric(temp$dose_toxicity_and_action_table[,2][temp$dose_toxicity_and_action_table[,1]==action][1])][which.min(dosevec1[dosevec1>as.numeric(temp$dose_toxicity_and_action_table[,2][temp$dose_toxicity_and_action_table[,1]==action][1])]-as.numeric(temp$dose_toxicity_and_action_table[,2][temp$dose_toxicity_and_action_table[,1]==action][1]))]

          if (which(dosevec==upper)-which(dosevec==lower)==1){
            dosevec <- c(dosevec[1:i], (dosevec[i]+dosevec[i+1])/2, dosevec[(i+1):dosenum])
            ptsvec <- c(ptsvec[1:i], 0, ptsvec[(i+1):dosenum])
            DLTvec <- c(DLTvec[1:i], 0, DLTvec[(i+1):dosenum])
            trueDLT <-c(trueDLT[1:i], (trueDLT[i]+trueDLT[i+1])/2, trueDLT[(i+1):dosenum])
            dosenum <- dosenum + 1
            imax <- imax + 1
            i <- min((i+1), imax)
          }
          else if (which(dosevec==upper)-which(dosevec==lower)>1){i <- i}
        }

        else if (action == "Next Dose"){i <- min((i+1), imax)}
      }

      if (sum(ptsvec)==nmax){break}  ## Stopping rule
    }


    if (imax==0){dselect[trial] <- -99}

    else {

      dselect[trial] <- hybrid_MTD_selection(target, dosevec0, ptsvec0, DLTvec0,elimdose)$MTD

    }

    MTDpts[trial] <- ptsvec[dosevec==trueMTD]
    if (dselect[trial]==trueMTD){ncorrect <- ncorrect+1}
    else if (dselect[trial]>trueMTD){overdose <- overdose+1}
    else if (dselect[trial]<trueMTD){underdose <- underdose+1}

  }

  avgMTDpts <- mean(MTDpts)
  pcorrect <- ncorrect/(ntrial)
  poverd <- overdose/(ntrial)
  punderd <- underdose/(ntrial)

  out = list(pcorrect=pcorrect, poverd=poverd, punderd=punderd, avgMTDpts=avgMTDpts)

  return(out);
}
