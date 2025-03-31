#' Constructor for linear2ph Objects
#'
#' Creates an object of class \code{linear2ph}.
#' @param input the internal list input
#' @return An object of class \code{linear2ph}.
#' @export
linear2ph_class <- function(input) {
  structure(
    list(
      call = input$call,
      coefficients = input$coefficients,
      sigma = input$sigma,
      covariance = input$covariance,
      converge = input$converge,
      converge_cov = input$converge_cov
    ),
    class = "linear2ph"
  )
}

#' Print Method for linear2ph Objects
#'
#' Prints the details of a \code{linear2ph} object.
#' @param object An object of class \code{linear2ph}.
#' @export
print.linear2ph <- function(object) {
  # print function call
  cat("Call:\n")
  print(match.call())
  cat("Call:\n")

  if (object$converge) {
    cat("The parameter estimation has converged.\n")
    cat("Coefficients:\n")
    print(object$coefficients[,1])
    if(!object$converge_cov){
      cat("The variance estimation is either not requested to be estimated or did not converge.")
    }
  } else {
    cat("This model did not converge.\n")
  }
}

#' Coefficient Method for linear2ph Objects
#'
#' Prints the coefficientss of a \code{linear2ph} object.
#' @param object An object of class \code{linear2ph}.
#' @export
coef.linear2ph <- function(object) {
  object$coefficients[,1]
}

#' Coefficient Method for linear2ph Objects
#'
#' Prints the coefficientss of a \code{linear2ph} object.
#' @param object An object of class \code{linear2ph}.
#' @export
coefficients.linear2ph <- function(object) {
  object$coefficients[,1]
}

#' Summary Method for linear2ph Objects
#'
#' Summarizes the details of a \code{linear2ph} object.
#' @param object An object of class \code{linear2ph}.
#' @return An object of class \code{summary.linear2ph}, containing the call, coefficients, and covariance.
#' @export
#' @method summary linear2ph
summary.linear2ph <- function(object) {
  if (!object$converge) {
    warning("This model did not converge. No summary available.")
    return(invisible(NULL))
  }

  # Construct summary object similar to summary.lm
  summary_obj <- list(
    call = object$call,  # Store original model call
    coefficients = object$coefficients,
    covariance = object$covariance
  )

  class(summary_obj) <- "summary.linear2ph"  # Assign S3 class
  return(summary_obj)  # Unlike print, summary should return a usable object
}

#' Print Method for summary.linear2ph Objects
#'
#' Prints a structured summary of a \code{linear2ph} model.
#' @param x An object of class \code{summary.linear2ph}.
#' @return Invisibly returns \code{x}.
#' @export
#' @method print summary.linear2ph
print.summary.linear2ph <- function(x) {
  if (!inherits(x, "summary.linear2ph")) {
    stop("print.summary.linear2ph() called on a non-summary object")
  }

  # Print function call like summary.lm
  cat("Call:\n")
  print(x$call)
  cat("\n")

  # Print model coefficients
  cat("Coefficients:\n")
  print(x$coefficients)

  invisible(x)  # Prevent automatic printing when assigned
}


#' Constructor for logistic2ph Objects
#'
#' Creates an object of class \code{logistic2ph}.
#' @param input the internal list input
#' @return An object of class \code{logistic2ph}.
#' @export
logistic2ph_class <- function(input) {
  structure(
    list(
      call = input$call,
      coefficients = input$coefficients,
      covariance = input$covariance,
      converge = input$converge,
      converge_cov = input$converge_cov
    ),
    class = "logistic2ph"
  )
}

#' Print Method for logistic2ph Objects
#'
#' Prints the details of a \code{logistic2ph} object.
#' @param object An object of class \code{logistic2ph}.
#' @export
#' @method print logistic2ph
print.logistic2ph <- function(object) {
  # print function call
  cat("Call:\n")
  print(match.call())
  cat("Call:\n")

  if (object$converge) {
    cat("The parameter estimation has converged.\n")
    cat("Coefficients:\n")
    print(object$coefficients[,1])
    if(!object$converge_cov){
      cat("The variance estimation is either not requested to be estimated or did not converge.")
    }
  } else {
    cat("This model did not converge.\n")
  }
}

#' Coefficient Method for logistic2ph Objects
#'
#' Prints the coefficients of a \code{logistic2ph} object.
#' @param object An object of class \code{logistic2ph}.
#' @export
#' @method coef logistic2ph
coef.logistic2ph <- function(object) {
  object$coefficients[,1]
}

#' Coefficient Method for logistic2ph Objects
#'
#' Prints the coefficients of a \code{logistic2ph} object.
#' @param object An object of class \code{logistic2ph}.
#' @export
#' @method coef logistic2ph
coefficients.logistic2ph <- function(object) {
  object$coefficients[,1]
}

#' Summary Method for logistic2ph Objects
#'
#' Summarizes the details of a \code{logistic2ph} object.
#' @param object An object of class \code{logistic2ph}.
#' @return An object of class \code{summary.logistic2ph}, containing the call, coefficients, and covariance.
#' @export
#' @method summary logistic2ph
summary.logistic2ph <- function(object) {
  if (!object$converge) {
    warning("This model did not converge. No summary available.")
    return(invisible(NULL))
  }

  # Construct summary object similar to summary.lm
  summary_obj <- list(
    call = object$call,  # Store original model call
    coefficients = object$coefficients,
    covariance = object$covariance
  )

  class(summary_obj) <- "summary.logistic2ph"  # Assign S3 class
  return(summary_obj)  # Unlike print, summary should return a usable object
}

#' Print Method for summary.logistic2ph Objects
#'
#' Prints a structured summary of a \code{logistic2ph} model.
#' @param x An object of class \code{summary.logistic2ph}.
#' @return Invisibly returns \code{x}.
#' @export
#' @method print summary.logistic2ph
print.summary.logistic2ph <- function(x) {
  if (!inherits(x, "summary.logistic2ph")) {
    stop("print.summary.logistic2ph() called on a non-summary object")
  }

  # Print function call like summary.lm
  cat("Call:\n")
  print(x$call)
  cat("\n")

  # Print model coefficients
  cat("Coefficients:\n")
  print(x$coefficients)

  invisible(x)  # Prevent automatic printing when assigned
}
