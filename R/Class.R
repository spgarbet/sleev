#' Constructor for linear2ph Objects
#'
#' Creates an object of class \code{linear2ph}.
#' @param input the internal list input
#' @return An object of class \code{linear2ph}.
#' @export
linear2ph_class <- function(input) {
  structure(
    list(
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
  object$coefficients
}

#' Summary Method for linear2ph Objects
#'
#' Summarizes the details of a \code{linear2ph} object.
#' @param object An object of class \code{linear2ph}.
#' @return coefficient matrix and covariance matrix.
#' @export
#' @method summary linear2ph
summary.linear2ph <- function(object) {
  if (!object$converge) {
    warning("This model did not converge. No summary available.")
    return(invisible(NULL))
  }

  # Construct summary object
  summary_obj <- list(
    coefficients = object$coefficients,
    covariance = object$covariance
  )

  class(summary_obj) <- "summary.linear2ph"  # Assign class
  return(invisible(summary_obj))  # Return invisibly
}

#' Print Method for summary.linear2ph Objects
#'
#' Prints a structured summary of a \code{linear2ph} model.
#'
#' @param x An object of class \code{summary.linear2ph}.
#' @return Invisibly returns \code{x}.
#' @export
#' @method print summary.linear2ph
print.summary.linear2ph <- function(x) {
  cat("Model Summary:\n")
  cat("Coefficients:\n")
  print(x$coefficients)  # Explicitly prints coefficients

  cat("\nCovariance Matrix:\n")
  print(x$covariance)  # Explicitly prints covariance

  invisible(x)  # Return invisibly, avoiding auto-printing
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
  object$coefficients
}

#' Summary Method for logistic2ph Objects
#'
#' Summarizes the details of a \code{logistic2ph} object.
#' @param object An object of class \code{logistic2ph}.
#' @return coefficient matrix and covariance matrix.
#' @export
#' @method summary logistic2ph
summary.logistic2ph <- function(object) {
  if (!object$converge) {
    warning("This model did not converge. No summary available.")
    return(invisible(NULL))
  }

  # Construct summary object
  summary_obj <- list(
    coefficients = object$coefficients,
    covariance = object$covariance
  )

  class(summary_obj) <- "summary.logistic2ph"  # Assign class
  return(invisible(summary_obj))  # Return invisibly
}

#' Print Method for summary.logistic2ph Objects
#'
#' Prints a structured summary of a \code{logistic2ph} model.
#'
#' @param x An object of class \code{summary.logistic2ph}.
#' @return Invisibly returns \code{x}.
#' @export
#' @method print summary.logistic2ph
print.summary.logistic2ph <- function(x) {
  cat("Model Summary:\n")
  cat("Coefficients:\n")
  print(x$coefficients)  # Explicitly prints coefficients

  cat("\nCovariance Matrix:\n")
  print(x$covariance)  # Explicitly prints covariance

  invisible(x)  # Return invisibly, avoiding auto-printing
}
