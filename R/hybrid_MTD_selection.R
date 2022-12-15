#'
#' Select the maximum tolerated dose (MTD) for single-agent dose-finding studies
#'
#' Select the maximum tolerated dose (MTD) when the single-agent dose-finding study is completed
#'
#' @usage hybrid_MTD_selection(target, dose, npts, nDLT, elimdose)
#'
#' @param target the target toxicity rate
#' @param dose a vector containing the numerical dosage of each dose level
#' @param npts a vector containing the number of patients treated at each dose level
#' @param nDLT a vector containing the number of patients who experienced dose-limiting
#'              toxicity at each dose level
#' @param elimdose the dosage at the dose level which is excluded due to excessive toxicity
#'
#' @details \code{hybrid.MTD.selection()} selects the MTD based on isotonic estimates of toxicity
#'          probabilities. The isotonic estimates are obtained by the pooled-adjacent-violators
#'          algorithm (PAVA) (Barlow, 1972 <doi: 10.1080/01621459.1972.10481216>).
#'
#' @note  The dose levels above \code{elim} are all excluded for MTD selection.
#'
#' @return The selected dosage as MTD
#'
#' @examples
#' hybrid_MTD_selection(target=0.3, dose=c(2,4,8,16,22,28,40), npts=c(2,4,8,16,22,28,40),
#'                      nDLT=c(0,0,0,0,1,0,2), elimdose=28)
#'
#' @export

hybrid_MTD_selection <- function(target, dose, npts, nDLT, elimdose){
  pava <- function(x, wt = rep(1, length(x))) {
    n <- length(x)
    if (n <= 1)
      return(x)
    if (any(is.na(x)) || any(is.na(wt))) {
      stop("Missing values in 'x' or 'wt' not allowed")
    }
    lvlsets <- (1:n)
    repeat {
      viol <- (as.vector(diff(x)) < 0)
      if (!(any(viol)))
        break
      i <- min((1:(n - 1))[viol])
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
      lvlsets[ilvl] <- lvl1
    }
    x
  }

  elimvec <- rep(0, length(dose))

  if (elimdose != 0){
    elipos <- which(dose==elimdose)

    for (i in elipos:length(dose)){
      elimvec[i] <- 1
    }
  }


  if (elimvec[1] == 1) {
    selectdose = -99
  }

  else {
    adm.set = (elimvec == 0)
    adm.index = which(adm.set == T)
    nDLT.adm = nDLT[adm.set]
    npts.adm = npts[adm.set]
    phat = (nDLT.adm + 0.05)/(npts.adm + 0.1)
    phat.var = (nDLT.adm + 0.05) * (npts.adm - nDLT.adm + 0.05)/((npts.adm +
                                                                    0.1)^2 * (npts.adm + 0.1 + 1))
    phat = pava(phat, wt = 1/phat.var)
    phat = phat + (1:length(phat)) * 1e-10

    selectd = sort(abs(phat - target), index.return = T)$ix[1]
    selectdose = adm.index[selectd]

  }
  out = list(MTD = dose[selectdose])
  return(out)

}
