#' @export
live_fa2 <- function( m, iter=100, step=10, ... ) {
    pos
    for( e in 1:floor( iter / step ) ) {
        pos <- forceatlas2( m, ..., iter=step, init_pos=pos )
        plot( pos, xlab=NA, ylab=NA, axes=FALSE )
    }
}

