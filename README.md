# Rforceatlas

R/Rcpp implementation of Gephi's ForceAtlas2 graph layout algorithm.

Detailed discussion of the procedure and its parameters can be found in [Jacomy M, Venturini T, Heymann S, Bastian M (2014) ForceAtlas2, a Continuous Graph Layout Algorithm for Handy Network Visualization Designed for the Gephi Software. PLoS ONE 9(6): e98679](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0098679).

This implementation has been written from scratch in C++ using Eigen for all numerical computations.

All functions are implemented serially, support for parallel computation is planned but not implemented yet.

## Installation

```R
#install.packages( 'devtools' )
install_github( "jtatria/Rforceatlas" )
```

## Usage

The package is loaded as usual with

```R
require( Rforceatlas )
```

This will export the following functions:
* `forceatlas` is the actual algorithm implementation in C++.
* `forceatlas_live` is simple R wrapper allowing for plotting intermediate convergence steps (mimicking Gephi's original implementation).
* `layout_with_fa2`, `with_fa2`, `layout.forceatlas2` igraph compatibility functions that take igraph objects as inputs.

A brief review of the main algorithm parameters can be found in the expoted function's documentation.

Both `forceatlas` and `forceatlas_live` take an adjacency matrix as input and return a matrix of vertex positions in the requested dimensions.

igraph functions take igraph objects as inputs and extract the adjacency matrix before calling the main work function.

Typical usage:

```R
require( Rforceatlas )
m <- matrix( ( runif( 100*100 ) > .3 ), 100, 100 )
lo <- forceatlas( m )
plot( lo )

require( igraph )
g <- sample_grg( 100, .3 )
lo <- layout_with_fa2( g )
plot( g, layout=lo )
```

# Author

José Tomás Atria (jtatria@gmail.com)

# License

GPL v3 or later.
