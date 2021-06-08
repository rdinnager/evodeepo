#' Convert an object with genomic data to a `torch_tensor`
#'
#' @param x Object to be converted
#' @param ... Other arguments passed to or from other methods
#' @inheritParams adegenet::tab
#' @inheritParams torch::torch_tensor
#'
#' @return A `torch_tensor`
#' @export
#'
#' @examples
as_tensor <- function(x, dtype = torch::torch_float(),
                      device = NULL,
                      requires_grad = FALSE,
                      NA.method = c("asis", "mean", "zero"),
                      ...) {
  UseMethod("as_tensor", x)
}

#' Convert a `genlight` object from the `adegenet` package to a `torch_tensor`
#'
#' @param x `genlight` object to convert
#' @param by_pop Should the object be split by population?
#' @param ... Other arguments passed to or from other methods
#' @inheritParams adegenet::tab
#' @inheritParams torch::torch_tensor
#'
#' @return A `torch_tensor`
#' @export
#'
#' @examples
as_tensor.genlight <- function(x, dtype = torch::torch_float(),
                               device = NULL,
                               requires_grad = FALSE,
                               NA.method = c("asis", "mean", "zero"),
                               by_pop = FALSE, ...) {

  NA.method <- match.arg(NA.method)

  if(by_pop) {
    pops <- adegenet::pop(x)
  }

  x <- adegenet::tab(x, NA.method = NA.method)

  if(NA.method == "asis") {
    x[is.na(x)] <- NaN
  }

  if(by_pop & !is.null(pops) & length(pops) > 1) {
    x <- split.data.frame(x, pops) %>%
      purrr::map(~torch_tensor(.x, dtype = dtype))
  } else {
    torch::torch_tensor(as.matrix(x),
                        dtype = dtype,
                        device = device,
                        requires_grad = requires_grad)
  }


}
