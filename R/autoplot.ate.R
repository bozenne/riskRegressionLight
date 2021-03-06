### autoplot.ate.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: apr 28 2017 (14:19) 
## Version: 
## last-updated: Mar  3 2019 (19:45) 
##           By: Thomas Alexander Gerds
##     Update #: 58
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * autoplot.ate (documentation)
#' @title Plot Average Risks
#' @description Plot average risks.
#' @name autoplot.ate
#' 
#' @param object Object obtained with the function \code{ate}.
#' @param ci [logical] If \code{TRUE} display the confidence intervals for the average risks.
#' @param band [logical] If \code{TRUE} display the confidence bands for the average risks.
#' @param plot [logical] Should the graphic be plotted.
#' @param digits [integer, >0] Number of decimal places.
#' @param alpha [numeric, 0-1] Transparency of the confidence bands. Argument passed to \code{ggplot2::geom_ribbon}.
#' @param ... not used. Only for compatibility with the plot method.
#' 
#' @seealso
#' \code{\link{ate}} to compute average risks.

## * autoplot.ate (examples)
#' @rdname autoplot.ate
#' @examples
#' \dontrun{
#' library(survival)
#' library(rms)
#' library(ggplot2)
#' #### simulate data ####
#' n <- 1e2
#' set.seed(10)
#' dtS <- sampleData(n,outcome="survival")
#' 
#' 
#' #### Cox model ####
#' fit <- cph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
#'
#' #### Average treatment effect ####
#' seqTimes <- sort(unique(fit$y[,1]))
#' seqTimes5 <- seqTimes[seqTimes>5 & seqTimes<10]
#' ateFit <- ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
#'               times = seqTimes, B = 0, band = TRUE, nsim.band = 500, y = TRUE,
#'               mc.cores=1)
#'
#' #### display #### 
#' autoplot(ateFit)
#' 
#' outGG <- autoplot(ateFit, band = TRUE, ci = TRUE, alpha = 0.1)
#' dd <- as.data.frame(outGG$data[Treatment == 0])
#' outGG$plot + facet_wrap(~Treatment, labeller = label_both)
#' }

## * autoplot.ate (code)
#' @rdname autoplot.ate
#' @method autoplot ate
#' @export
autoplot.ate <- function(object,
                         ci = FALSE,
                         band = FALSE,
                         plot = TRUE,
                         digits = 2,
                         alpha = NA, ...){

    ## for CRAN check
    Treatment <- NULL
    
    ## initialize and check          
    if(ci[[1]]==TRUE && (object$se[[1]]==FALSE || is.null(object$conf.level))){
        stop("argument \'ci\' cannot be TRUE when no standard error have been computed \n",
             "set arguments \'se\' and \'confint\' to TRUE when calling ate \n")
    }
    if(band[[1]] && (object$band[[1]]==FALSE  || is.null(object$conf.level))){
        stop("argument \'band\' cannot be TRUE when the quantiles for the confidence bands have not been computed \n",
             "set arguments \'band\' and \'confint\' to TRUE when calling ate \n")
    }
    
    dots <- list(...)
    if(length(dots)>0){
        txt <- names(dots)
        txt.s <- if(length(txt)>1){"s"}else{""}
        stop("unknown argument",txt.s,": \"",paste0(txt,collapse="\" \""),"\" \n")
    }

    ## display
    dataL <- copy(object$meanRisk)
    dataL[,row := as.numeric(as.factor(Treatment))]
    if(ci){
        data.table::setnames(dataL,
                             old = c("meanRisk.lower","meanRisk.upper"),
                             new = c("lowerCI","upperCI"))
    }
    if(band){
        data.table::setnames(dataL,
                             old = c("meanRisk.lowerBand","meanRisk.upperBand"),
                             new = c("lowerBand","upperBand"))
    }

    gg.res <- predict2plot(dataL = dataL,
                           name.outcome = "meanRisk", # must not contain space to avoid error in ggplot2
                           ci = ci, band = band,
                           group.by = "Treatment",
                           conf.level = object$conf.level,
                           alpha = alpha,
                           ylab = "Average absolute risk")
  
    if(plot){
        print(gg.res$plot)
    }
  
  return(invisible(gg.res))
}



#----------------------------------------------------------------------
### autoplot.ate.R ends here
