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
    cat("This model has converged.\n")
    cat("Coefficients:\n")
    print(object$coefficients[,1])
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

#' Print Method for linear2ph Objects
#'
#' Prints the details of a \code{linear2ph} object.
#' @param object An object of class \code{linear2ph}.
#' @return None. Prints a summary of the object.
#' @export
#' @method summary linear2ph
summary.linear2ph <- function(object) {
  if (object$converge) {
    cat("Model Summary:\n")
    cat("Coefficients:\n")
    print(object$coefficients)
    cat("\nCovariance Matrix:\n")
    print(object$covariance)
    # Return the summary as a list (if desired)
    result <- list(
      coefficients = object$coefficients,
      covariance = object$covariance
    )
    return(invisible(result))  # Return the summary as a list
  } else {
    cat("This model did not converge. No summary available.\n")
  }
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
    cat("This model has converged.\n")
    cat("Coefficients:\n")
    print(object$coefficients[,1])
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
#' Prints the details of a \code{logistic2ph} object.
#' @param object An object of class \code{logistic2ph}.
#' @return list of coefficients and covariance
#' @export
#' @method summary logistic2ph
summary.logistic2ph <- function(object) {
  if (object$converge) {
    # Print the summary
    cat("Model Summary:\n")
    cat("Coefficients:\n")
    print(object$coefficients)  # Print coefficients
    cat("\nCovariance Matrix:\n")
    print(object$covariance)  # Print covariance matrix

    # Return the summary as a list (if desired)
    result <- list(
      coefficients = object$coefficients,
      covariance = object$covariance
    )
    return(invisible(result))  # Return the summary as a list
  } else {
    cat("This model did not converge. No summary available.\n")
    return(NULL)  # Return NULL if model did not converge
  }
}
