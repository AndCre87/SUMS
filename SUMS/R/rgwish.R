## This function and the ones called are modifications of the originals presented in (Mohammadi et al. 2015)

# Sample from G-Wishart using Cpp routing rgwish.cpp
rgwish = function( n = 1, nu = 3, Psi = NULL, G = NULL, threshold = 1e-8 ){
  if( nu <= 2 )         stop( "For G-Wishart distribution parameter 'nu' must be more than 2" )
  if( is.null( G ) ) stop( "Adjacency matrix must be determined" )
  
  if( !is.matrix( G ) ) stop("The adjacency matrix of the graph must be in matrix form")
  
  if( ( sum( G == 0 ) + sum( G == 1 ) ) != ( nrow( G ) ^ 2 ) ) stop( "Element of matrix 'adj' must be 0 or 1" )
  
  G <- as.matrix( G )
  diag( G ) <- 0
  
  if( !isSymmetric( G ) )
  {
    G[ lower.tri( G ) ] <- 0
    G                   <- G + t( G )
  }
  
  p <- nrow( G )
  if( p < 1 ) stop( "'p' must be more than or equal to 1" )
  
  if( is.null( Psi )      ) Psi <- diag( p )
  if( !isSymmetric( Psi ) ) stop( "Matrix 'D' must be positive definite matrix." )
  if( nrow( Psi ) != p    ) stop( "Dimension of matrix G and D must be the same." )
  
  Ti = chol( solve( Psi ) )
  Omega  = matrix( 0, p, p )
  
  if( n > 1 ){
    samples = array( 0, dim = c( p, p, n ) )
    
    for( i in 1 : n ){
      if( sum( G ) == ( p * ( p - 1 ) ) ){
        Omega = rwish(nu, Psi)
      }else{
        Omega = rgwish_c(nu, Ti, G, threshold)
      }
      samples[ , , i ] = Omega		
    }
  }else{
    if( sum( G ) == ( p * ( p - 1 ) ) ){
      Omega = rwish(nu, Psi)
    }else{
      Omega = rgwish_c(nu, Ti, G, threshold)
    }
    samples = Omega	
  }
  
  return( samples )   
}

