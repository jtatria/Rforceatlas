require( Rforceatlas )
require( igraph )
require( igraphdata )

data( "UKfaculty" )
#g <- UKfaculty
g <- sample_pa( 100 )
m <- igraph::as_adjacency_matrix( g, type='both' )
set.seed( 13 )
init <- matrix( runif( igraph::vcount( g )*2, -1, 1 ) * 1000, igraph::vcount( g ) )
r <- forceatlas_live( m, init=init, device=cairoDevice::Cairo, step=1, edges=TRUE, iter=10000, axes=TRUE )
