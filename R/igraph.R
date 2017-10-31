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
#' @param live Logical. If TRUE, attempt to plot intermediate convergence steps.
#' @export
layout_with_fa2 <- function( g, attr='weight', live=FALSE, ... ) {
    attr <- if( is.null( igraph::edge_attr( g, attr ) ) ) NULL else attr
    m <- igraph::as_adjacency_matrix( g, type='both', attr=attr, sparse=FALSE )
    res <- if( live ) {
        forceatlas_live( m, ... )
    } else {
        forceatlas( m, ... )
    }
    return( res )
}

#' @rdname forceatlas
#' @param ... Passed to \code{layout_with_fa2}
#' @export
with_fa2 <- function( ... ) igraph_layout_spec( layout_with_fa2, ... )

#' Deprecated layout functions.
#'
#' This is here only for compatibility with igraph's version < 0.8.0 layout API. Please use
#' \code{layout_with_fa2}.
#'
#' @param ...    Passed to \code{layout_with_f2}.
#' @param params Passed to \code{layout_with_fa2} as arguments.
#'
#' @export
layout.forceatlas2 <- function( ..., params=list() ) {
    igraph_do_call( layout_with_fa2, .args=c( list( ... ), params ) )
}

# TODO: ask for igraph to export utility functions and C headers!
igraph_layout_spec <- function( fun, ... ) {
    my_call <- match.call(sys.function(1), sys.call(1))
    my_call[[1]] <- substitute(fun)
    structure(
        list(
            fun = fun,
            call_str = sub("(", "(<graph>, ", deparse(my_call), fixed = TRUE),
            args = list(...)
        ), class = "igraph_layout_spec"
    )
}

# TODO: ask for igraph to export utility functions and C headers!
igraph_do_call <- function( f, ..., .args = list(), .env = parent.frame() ) {
    f <- substitute( f )
    call <- make_call( f, ..., .args )
    eval( call, .env )
}
