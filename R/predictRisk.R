### predictRisk.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun  6 2016 (09:02) 
## Version: 
## last-updated: maj  6 2019 (15:43) 
##           By: Brice Ozenne
##     Update #: 152
#----------------------------------------------------------------------
## 
### Commentary:
#
# methods for extracting predictions from various objects
# 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
# --------------------------------------------------------------------
#' Extrating predicting risks from regression models 
#' 
#' Extract event probabilities from fitted regression models and machine learning objects.
#' 
#' The function predictRisk is a generic function, meaning that it invokes
#' specifically designed functions depending on the 'class' of the first
#' argument.
#' 
#' See \code{\link{predictRisk}}.
#' 
#' @aliases predictRisk predictRisk.CauseSpecificCox
#' predictRisk.riskRegression predictRisk.FGR
#' predictRisk.prodlim predictRisk.rfsrc predictRisk.aalen
#' predictRisk.riskRegression predictRisk.cox.aalen
#' predictRisk.coxph predictRisk.cph predictRisk.default
#' predictRisk.matrix predictRisk.pecCtree
#' predictRisk.pecCforest predictRisk.prodlim predictRisk.psm
#' predictRisk.selectCox predictRisk.survfit predictRisk.randomForest
#' predictRisk.lrm predictRisk.glm
#' predictRisk.rpart
#' @usage
#' \method{predictRisk}{glm}(object,newdata,...)
#' \method{predictRisk}{cph}(object,newdata,times,...)
#' \method{predictRisk}{coxph}(object,newdata,times,...)
#' \method{predictRisk}{matrix}(object,newdata,times,cause,...)
#' \method{predictRisk}{survfit}(object,newdata,times,...)
#' \method{predictRisk}{prodlim}(object,newdata,times,cause,...)
#' \method{predictRisk}{CauseSpecificCox}(object,newdata,times,cause,...)
#' @param object A fitted model from which to extract predicted event
#' probabilities
#' @param newdata A data frame containing predictor variable combinations for
#' which to compute predicted event probabilities.
#' @param times A vector of times in the range of the response variable, for
#' which the cumulative incidences event probabilities are computed.
#' @param cause Identifies the cause of interest among the competing events.
#' @param \dots Additional arguments that are passed on to the current method.
#' @return For binary outcome a vector with predicted risks. For survival outcome with and without
#' competing risks
#' a matrix with as many rows as \code{NROW(newdata)} and as many
#' columns as \code{length(times)}. Each entry is a probability and in
#' rows the values should be increasing.
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
#' @details
#' In uncensored binary outcome data there is no need to choose a time point.
#'
#' When operating on models for survival analysis (without competing risks) the function still
#' predicts the risk, as 1 - S(t|X) where S(t|X) is survival chance of a subject characterized
#' by X.
#'
#' When there are competing risks (and the data are right censored) one needs
#' to specify both the time horizon for prediction (can be a vector) and the
#' cause of the event. The function then extracts the absolute risks F_c(t|X)
#' aka the cumulative incidence of an event of type/cause c until time t for a
#' subject characterized by X. Depending on the model it may or not be possible
#' to predict the risk of all causes in a competing risks setting. For example. a
#' cause-specific Cox (CSC) object allows to predict both cases whereas a Fine-Gray regression
#' model (FGR) is specific to one of the causes. 
#' 
#' @keywords survival
#' @export 
predictRisk <- function(object,newdata,...){
  UseMethod("predictRisk",object)
}

##' @export 
predictRisk.default <- function(object,newdata,times,cause,...){
    stop(paste0("No method available for evaluating predicted probabilities from objects in class: ",class(object),". But, you can write it yourself or ask the package manager."),call.=FALSE)
}

##' @export
predictRisk.double <- function(object,newdata,times,cause,...){
    stopifnot(NROW(object)==NROW(newdata))
    object
}

##' @export
predictRisk.integer <- function(object,newdata,times,cause,...){
    stopifnot(NROW(object)==NROW(newdata))
    object
}

##' @export
predictRisk.factor <- function(object,newdata,times,cause,...){
    stopifnot(NROW(object)==NROW(newdata))
    object
}

##' @export
predictRisk.numeric <- function(object,newdata,times,cause,...){
    stopifnot(NROW(object)==NROW(newdata))
    object
}

##' @export
predictRisk.glm <- function(object,newdata,...){
    if (object$family$family=="binomial")
        return(as.numeric(stats::predict(object,newdata=newdata,type="response")))
    else{ stop("Currently only the binomial family is implemented for predicting a status from a glm object.")
      }
}

##' @export
predictRisk.formula <- function(object,newdata,...){
    ff <- update.formula(object,"NULL~.")
    if (length(all.vars(ff))==1){
        p <- stats::model.frame(ff,newdata)[[1]]
        p
    } else{
        fit <- glm(object,data=newdata,family="binomial")
        predictRisk(fit,newdata=newdata,...)
    }
}

##' @export
predictRisk.matrix <- function(object,newdata,times,cause,...){
    if (NROW(object) != NROW(newdata) || NCOL(object) != length(times)){
        stop(paste("Prediction matrix has wrong dimensions: ",
                   NROW(object),
                   " rows and ",
                   NCOL(object),
                   " columns.\n But requested are predicted probabilities for ",
                   NROW(newdata),
                   " subjects (rows) in newdata and ",
                   length(times),
                   " time points (columns)",
                   sep=""))
    }
    object
}
   
##' @export
predictRisk.coxph <- function(object,newdata,times,...){
    p <- predictCox(object=object,
                    newdata=newdata,
                    times=times,
                    se = FALSE,
                    iid = FALSE,
                    keep.times=FALSE,
                    type="survival")$survival

    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)){
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    }
    return(1-p)
}

##' @export
predictRisk.coxphTD <- function(object,newdata,times,landmark,...){
    stopifnot(attr(object$y,"type")=="counting")
    bh <- survival::basehaz(object,centered=TRUE)
    Lambda0 <- bh[,1]
    etimes <- bh[,2]
    lp <- predict(object,newdata=newdata,type="lp")
    p <- do.call("cbind",lapply(times,function(ttt){
        index <- prodlim::sindex(eval.times=c(landmark,landmark+ttt),jump.times=etimes)
        Lambda0.diff <- c(0,Lambda0)[1+index[2]] - c(0,Lambda0)[1+index[1]]
        1-exp(-Lambda0.diff * exp(lp))
    }))
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)){
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    }
    return(p)
}

##' @export 
predictRisk.cph <- function(object,newdata,times,...){
    ## if (!match("surv",names(object),nomatch=0)) stop("Argument missing: set surv=TRUE in the call to cph!")
    ## p <- rms::survest(object,times=times,newdata=newdata,se.fit=FALSE,what="survival")$surv
    ## if (is.null(dim(p))) p <- matrix(p,nrow=NROW(newdata))
    p <- predictCox(object=object,
                    newdata=newdata,
                    times=times,
                    se = FALSE,
                    iid = FALSE,
                    keep.times=FALSE,
                    type="survival")$survival
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    return(1-p)
}

##' @export 
predictRisk.prodlim <- function(object,newdata,times,cause,...){
    ## require(prodlim)
    if (object$model[[1]]=="competing.risks" && missing(cause))
        stop(paste0("Cause is missing. Should be one of the following values: ",paste(attr(object$model.response,"states"),collapse=", ")))
    p <- predict(object=object,
                 cause=cause,
                 type="cuminc",
                 newdata=newdata,
                 times=times,
                 mode="matrix",
                 level.chaos=1)
    ## if the model has no covariates
    ## then all cases get the same prediction
    ## in this exceptional case we return a vector
    if (NROW(p)==1 && NROW(newdata)>=1)
        p <- as.vector(p)
    ## p[is.na(p)] <- 0
    if (is.null(dim(p))){
        if (length(p)!=length(times))
            stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    }
    else{
        if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
            stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    }
    p
}

predict.survfit <- function(object,newdata,times,bytimes=TRUE,type="cuminc",fill="last",...){
    if (length(class(object))!=1 || class(object)[[1]]!="survfit" || object$type[[1]] !="right")
        stop("Predictions only available \nfor class 'survfit', possibly stratified Kaplan-Meier fits.\n For class 'cph' Cox models see survest.cph.")
    if (missing(newdata))
        npat <- 1
    else
        if (is.data.frame(newdata))
            npat <- nrow(newdata)
        else stop("If argument `newdata' is supplied it must be a dataframe." )
    ntimes <- length(times)
    sfit <- summary(object,times=times)
    if (is.na(fill))
        Fill <- function(x,len){x[1:len]}
    else if (fill=="last")
        Fill <- function(x,len){
            y <- x[1:len]
            y[is.na(y)] <- x[length(x)]
            y}
    else stop("Argument fill must be the string 'last' or NA.")
    if (is.null(object$strata)){
        pp <- Fill(sfit$surv,ntimes)
        p <- matrix(rep(pp,npat),
                    ncol=ifelse(bytimes,ntimes,npat),
                    nrow=ifelse(bytimes,npat,ntimes),
                    byrow=bytimes)
    }
    else{
        covars <- attr(terms(eval.parent(object$call$formula)),"term.labels")
        if (!all(match(covars,names(newdata),nomatch=FALSE)))
            stop("Not all strata defining variables occur in newdata.")
        ## FIXME there are different ways to build strata levels
        ## how can we test which one was used???
        stratdat <- newdata[,covars,drop=FALSE]
        names(stratdat) <- covars
        NewStratVerb <- survival::strata(stratdat)
        NewStrat <- interaction(stratdat,sep=" ")
        levs <- levels(sfit$strata)
        #    print(levs)
        #    print(levels(NewStrat))
        #    print(levels(NewStratVerb))
        if (!all(choose <- match(NewStratVerb,levs,nomatch=F))
            &&
            !all(choose <- match(NewStrat,levs,nomatch=F)))
            stop("Not all strata levels in newdata occur in fit.")
        survlist <- split(sfit$surv,sfit$strata)
        pp <- lapply(survlist[choose],Fill,ntimes)
        p <- matrix(unlist(pp,use.names=FALSE),
                    ncol=ifelse(bytimes,ntimes,npat),
                    nrow=ifelse(bytimes,npat,ntimes),
                    byrow=bytimes)
    }
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    if (type=="cuminc") 1-p else p
}

##' @export
predictRisk.survfit <- function(object,newdata,times,...){
    p <- predict.survfit(object,newdata=newdata,times=times,type="cuminc",bytimes=TRUE,fill="last")
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

##' @export 
predictRisk.CauseSpecificCox <- function (object, newdata, times, cause, ...) { 
    p <- predict(object=object,
                 newdata=newdata,
                 times=times,
                 cause=cause,
                 keep.strata=FALSE,
                 se = FALSE,
                 iid = FALSE)$absRisk
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}





#----------------------------------------------------------------------
### predictRisk.R ends here
