#' @rdname g-and-k
#' @export
rgk = function(n, A, B, g, k, c=0.8){
    z2gk(stats::rnorm(n), A, B, g, k, c)
}
