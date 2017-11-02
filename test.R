require( Rforceatlas )
require( igraph )
require( igraphdata )

data( "UKfaculty" )
g <- UKfaculty
m <- igraph::as_adjacency_matrix( g, type='both', sparse=FALSE )
init <- matrix( runif( igraph::vcount( g )*2, -1, 1 ) * 1000, igraph::vcount( g ) )
r <- forceatlas_live( m, init=init, device=cairoDevice::Cairo, step=10, iter=10000, axes=TRUE, G=0.1 )

