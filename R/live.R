# Rforceatlas: Rcpp implementation of the ForceAtlas2 algorithm
# Copyright (C) 2017 José Tomás Atria jtatria at nomoi dot org
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

#' @rdname forceatlas
#' @param ... Passed to \code{forceatlas}
#' @export
forceatlas_live <- function(
    m, iter=100, step=10, dim=2,
    init=matrix( runif( nrow( m ) * dim, -1, 1 ) * 1000, nrow( m ), dim ), ...
) {
    pos <- init
    for( e in 1:floor( iter / step ) ) {
        last <- pos
        pos <- forceatlas( m, iter=step, init=pos, dim=2, ... )
        plot_step( m, pos, e*step, iter, edges=edges, xlab=NA, ylab=NA, axes=axes )
        if( !is.null( last ) && max( abs( last - pos ) ) == 0 ) {
            message( "Convergence reached!" )
            break
        }
    }
    plot_step( m, pos, e*step, iter, edges=TRUE, xlab=NA, ylab=NA, axes=axes )
    return( pos )
}

#' @importFrom Matrix which
plot_step <- function( m, pos, iter, total, edges=FALSE, lwd=.1, ecol='grey', ... ) {
    plot( NA, ..., pty='s',
        xlim=c( min( pos[,1] ), max( pos[,1] ) ), ylim=c( min( pos[,2] ), max( pos[,2] ) )
    )
    if( edges && any( m != 0 ) ) {
        e <- which( m != 0, arr.ind=TRUE )
        e0 <- pos[e[,1],]
        e1 <- pos[e[,2],]
        segments( e0[,1], e0[,2], e1[,1], e1[,2], lwd=lwd, col=ecol )
        # TODO extra info?
    }
    points( pos )
    title( main=sprintf( "ForceAtlas iteration %d out of %d", iter, total ) )
}
