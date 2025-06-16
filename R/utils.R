
#` @importFrom phutil as_persistence
#` @export
phutil::as_persistence

#` @importFrom phutil get_pairs
#` @export
phutil::get_pairs

# copied from {ggplot2}
`%||%` <- function (a, b) if (! is.null(a)) a else b

# pre-process power
ensure_p <- function(p) {
  # only allow positive integer powers
  if (p < 1 || (p != Inf && p %% 1 != 0))
    stop("`p` must be a positive integer or infinity.")
  p
}


