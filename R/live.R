#' @rdname forceatlas
#' @param ... Passed to \code{forceatlas}
#' @export
forceatlas_live <- function(
    m, iter=100, step=10, dim=2,
    init=matrix( runif( nrow( m ) * dim, -1, 1 ) * 1000, nrow( m ), dim ), ...
) {
    pos <- init
    for( e in 1:floor( iter / step ) ) {
        pos <- forceatlas( m, ..., iter=step, init=pos, dim=dim )
        plot( pos, xlab=NA, ylab=NA, axes=FALSE )
    }
}

