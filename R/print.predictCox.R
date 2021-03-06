### print.predictCox.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 15 2017 (17:36) 
## Version: 
## last-updated: jul  5 2018 (15:52) 
##           By: Brice Ozenne
##     Update #: 156
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * print.predictCox (documentation)
#' @title Print Predictions From a Cox Model
#' @description Print predictions from a Cox model.
#' @name print.predictCox
#' 
#' @param x object obtained with the function \code{predictCox}.
#' @param digits [integer, >0] indicating the number of decimal places.
#' @param ... Passed to print.
#' 
#' @details to display confidence intervals/bands,
#' the \code{confint} method needs to be applied on the object.
#'
#' @seealso
#' \code{\link{confint.predictCox}} to compute confidence intervals/bands.
#' \code{\link{predictCox}} to compute the predicted cumulative hazard/survival.


## * print.predictCox (code)
#' @rdname print.predictCox
#' @method print predictCox
#' @export
print.predictCox <- function(x, digits = 3,...){
    out <- as.data.table(x)
    print(out,digits=digits,...)
    invisible(out)
}



#----------------------------------------------------------------------
### print.predictCox.R ends here
