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
        pos <- forceatlas( m, ..., iter=step, init=pos, dim=dim )
        plot( pos, xlab=NA, ylab=NA, axes=FALSE )
    }
    return( pos )
}

