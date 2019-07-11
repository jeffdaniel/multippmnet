#' Fit a regularized inhomogeneous multivariate Poisson point process model
#'
#' Fit an inhomogeneous multivariate Poisson point process model via penalized
#' conditional logistic regression. An adaptive sparse group lasso penalty is
#' employed in order to exploit grouped structure in the model coefficients.
#'
#' @param Q A quadrature scheme (of class \code{quad}) containing a multivariate
#'        point pattern.
#' @param data A list of pixel images (of class \code{imlist}) representing
#'        spatial covariates.
#' @param ... Additional arguments passed to \code{asgl} to control the fitting
#'        procedure.
#' @param grp_weights Optional. A vector of adaptive weights applied to the
#'        groups of coefficients associated with each covariate in \code{data}.
#' @param ind_weights Optional. A vector of adaptive weights applied to the
#'        individual model coefficients.
#'
#' @return An object of class \code{multippmnet}.
#' @import spatstat
#' @importFrom stats as.formula
#' @importFrom stats model.matrix
#' @export
multippmnet <- function(Q, data, ...,
                        grp_weights = NULL, ind_weights = NULL) {

  # Validate point pattern data
  if (!inherits(Q, "quad")) {
    stop("The argument 'Q' must be a quadrature scheme.", call. = FALSE)
  } else {
    if (!is.multitype(Q)) {
      stop("The argument 'Q' must contain a multitype point pattern.",
           call. = FALSE)
    }
  }

  # Validate covariate data
  if (!inherits(data, "imlist")) {
    stop("The argument 'data' must be a list of pixel images.",
         call. = FALSE)
  } else {
    if (length(data) < 2) {
      stop("The argument 'data' must be a list of 2 or more pixel ",
           "images.", call. = FALSE)
    }
  }

  # Construct matrix of covariate values at quadrature points
  U <- union.quad(Q)
  x <- lapply(data, lookup.im, U$x, U$y, naok = TRUE, strict = FALSE)
  x <- do.call(cbind, x)

  fmla <- as.formula(paste0("~ marks / (",
                           paste0(names(data), collapse = " + "), ")"))
  x <- as.data.frame(x); x <- cbind(marks = marks(U), x)
  x <- model.matrix(fmla, x)[, -1]

  # Construct vectors of "responses" and offset
  b <- rep.int(Q$param$rho, nrow(x))
  y <- as.numeric(is.data(Q))

  # Ignore any quadrature points with NA covariate values
  subset <- as.vector(!(rowSums(is.na(x)) > 0))

  # Index
  index <- sort(rep(rep(1:(length(data) + 1)),
                    length(levels(marks(U)))))
  index <- index[-1]

  # Adaptive weights
  if (is.null(grp_weights)) {
    grp_weights <- rep(1, length(unique(index)))
    grp_weights[1] <- 0
  }
  if (is.null(ind_weights)) {
    ind_weights <- rep(1, ncol(x))
    ind_weights[1:(length(levels(marks(U))) - 1)] <- 0
  }

  # Fit the regularization path
  fit <- asgl::asgl(x[subset, ], y[subset], index, family = "binomial",
                    offset = -log(b[subset]), standardize = FALSE,
                    grp_weights = grp_weights, ind_weights = ind_weights,
                    ...)

  fit$Q <- Q
  fit$x <- x
  fit$y <- y
  fit$b <- b
  fit$fmla <- fmla
  fit$index <- index
  fit$subset <- subset
  fit$grp_weights <- grp_weights
  fit$ind_weights <- ind_weights
  class(fit) <- c("multippmnet", class(fit))
  fit
}
