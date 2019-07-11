#' Make predictions from a 'multippmnet' object
#'
#' Make predictions from a '\code{multippmnet}' object.
#'
#' @param object A fitted \code{multippmnet} object.
#' @param data A list of pixel images (of class \code{imlist}) containing the
#'         spatial covariates used to fit the model.
#' @param window Optional. An observation window (of class \code{owin}) defining
#'        the region within which predictions are to be made. Default is the
#'        window of the original data used to fit the model.
#' @param eps Optional. The height and width of pixels in the prediction
#'        image(s). A numeric value or numeric vector of length 2 specifying
#'        pixel dimensions in the x and y directions. Incompatible with
#'        \code{dimyx}.
#' @param dimyx Optional. The resolution of the prediction image(s). A numeric
#'        value or numeric vector of length 2 specifying the number of pixels
#'        in the y and x directions. Incompatible with \code{eps}.
#' @param s Value(s) of the penalty tuning parameter at which predictions are to
#'        be made. Default is the entire sequence used to fit the regularization
#'        path.
#' @param ... Ignored.
#' @export
predict.multippmnet <- function(object, data, window = NULL,
                                eps = NULL, dimyx = NULL, s = NULL, ...) {

  # Get window for predictions
  if (is.null(window)) {
    masque <- as.mask(object$Q$data$window, eps, dimyx)
  } else {
    masque <- as.mask(window, eps, dimyx)
  }

  # Get prediction points
  rxy <- rasterxy.mask(masque, drop = TRUE)
  xpredict <- rxy$x
  ypredict <- rxy$y

  # Make a copy of each prediction point for each type of point
  types <- levels(marks(object$Q))
  nt <- length(types)
  np <- length(xpredict)
  xpredict <- rep.int(xpredict, nt)
  ypredict <- rep.int(ypredict, nt)
  mpredict <- rep.int(types, rep.int(np, nt))

  # Construct matrix of covariate values at prediction points
  newx <- lapply(data, lookup.im, xpredict, ypredict,
                 naok = TRUE, strict = FALSE)
  newx <- do.call(cbind, newx)
  newx <- as.data.frame(newx)
  newx <- data.frame(marks = mpredict, newx)

  fmla  <- object$fmla
  newx <- model.matrix(fmla, data = newx)[,-1]

  preds <- asgl:::predict.asgl(object, newx = newx, s = s, type = "link")
  preds <- exp(preds)

  output <- as.list(seq(ncol(preds)))
  for (i in 1:length(output)) {
    predlist <- as.list(seq(nt))
    for (k in 1:length(predlist)) {
      out <- as.im(masque)
      out[] <- preds[mpredict == types[k], i]
      predlist[[k]] <- out
    }
    names(predlist) <- types
    output[[i]] <- as.imlist(predlist)
    if (length(output) == 1) {
      output <- output[[1]]
    }
  }
  return(output)

}

#' Fitted intensity from a 'multippmnet' object
#'
#' Computes the fitted intensity for a regularized inhomogeneous Poisson point
#' process model at the quadrature points used to fit the model.
#'
#' @param object A fitted \code{multippmnet} object.
#' @param s Value(s) of the penalty tuning parameter at which predictions are
#'        to be made. Default is the entire sequence used to fit the
#'        regularization path.
#' @param drop Logical value. If \code{TRUE}, quadrature points that were not
#'        used to fit the model are deleted.
#' @param ... Ignored
#'
#' @export
fitted.multippmnet <- function(object, s = NULL, drop = FALSE, ...) {
  x <- object$x
  if (drop) {
    x <- x[object$subset, ]
  }
  result <- asgl:::predict.asgl(object, newx = x, s = s, type = "link")
  result <- exp(result)
  return(result)
}
