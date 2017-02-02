#' @rdname g-and-h
#' @export
rgh = function(n, A, B, g, h, c=0.8){
    z2gh(stats::rnorm(n), A, B, g, h, c)
}
