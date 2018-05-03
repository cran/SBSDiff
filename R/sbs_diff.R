#' @title Satorra-Bentler Scaled Chi-Squared Difference Test (Based on Chi-Squared Values)
#' @description Takes chi-squared values from nested models estimated using maximum likelihood
#' with robust standard errors, model degrees of freedom, scaling correlation factors and returns:
#' (1) change in model chi-squared (2) change in model degrees of freedom and
#' (3) the probability of rejecting the null.
#' @param chi0 chi-squared value for the more restrictive model
#' @param chi1 chi-squared value for the less restrictive model
#' @param df0 degrees of freedom for the more restrictive model (with more degrees of freedom)
#' @param df1 degrees of freedom for the less restrictive model (with fewer degrees of freedom)
#' @param c0 scaling correction factor for the more restrictive model
#' @param c1 scaling correction factor for the less restrictive model
#' @return Change in model chi-squared, change in model degrees of freedom and the probability of rejecting the null
#' @examples
#' chi0 <- 50
#'
#' chi1 <- 40
#'
#' df0 <- 10
#'
#' df1 <- 9
#'
#' c0 <- 1
#'
#' c1 <- 1
#'
#' sbs.chi(chi0,chi1,df0,df1,c0,c1)
#' @export

sbs.chi<-function(chi0,chi1,df0,df1,c0,c1){
  cd<-(df0*c0-df1*c1)/(df0-df1)
  TRd<-(chi0*c0-chi1*c1)/cd
  delta.df<-df0-df1
  p<-1-pchisq(TRd,delta.df)
  results1<-list(TRd,delta.df,p)
  results1<-as.data.frame(results1)
  colnames(results1) <- c("SB-Scaled Diff","Delta DF","p")
  rownames(results1)<-""
  print(results1)}

#' @title Satorra-Bentler Scaled Chi-Squared Difference Test (Based on Loglikelihood Values)
#' @description Takes loglikelihood values from nested models estimated using maximum likelihood
#' with robust standard errors, number of free parameters, scaling correlation factors and returns:
#' (1) Satorra-Bentler scaled change in model chi-squared (2) change in model degrees of freedom and
#' (3) the probability of rejecting the null.
#' @param L0 loglikelihood value for the more restrictive model (should be a negatige value)
#' @param L1 loglikelihood value for the less restrictive model (should be a negatige value)
#' @param p0 number of free parameters for the more restrictive model (with fewer freely estimated parameters)
#' @param p1 number of free parametersfor the less restrictive model (with more freely estimated parameters)
#' @param c0 scaling correction factor for the more restrictive model
#' @param c1 scaling correction factor for the less restrictive model
#' @return Change in model chi-squared, change in model degrees of freedom and the probability of rejecting the null
#' @examples
#' L0 <- -50
#'
#' L1 <- -45
#'
#' p0 <- 9
#'
#' p1 <- 10
#'
#' c0 <- 1
#'
#' c1 <- 1
#'
#' sbs.log(L0,L1,p0,p1,c0,c1)
#' @export

sbs.log<-function(L0,L1,p0,p1,c0,c1){
  cd<-((p0 * c0) - (p1 * c1))/(p0-p1)
  TRd<- (-2*(L0 - L1))/cd
  delta.df<-p1-p0
  p<-1-pchisq(TRd,delta.df)
  results2<-list(TRd,delta.df,p)
  results2<-as.data.frame(results2)
  colnames(results2) <- c("SB-Scaled Diff","delta DF","p")
  rownames(results2)<-""
  print(results2)}


