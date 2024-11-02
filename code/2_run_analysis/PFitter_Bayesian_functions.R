    # This is PFitter's days function, modified for taking clambda instead of lambda.  We use the Bayesian logic that the unknown Poisson rate has a Gamma-distributed posterior, and with no observed failures, we just call the prior gamma(alpha=1,beta=1), so the posterior is gamma( alpha + 0, beta + #nobs ) -- but note that the original days(..) function expects the intersequence HD estimate of lambda, which should be about twice the estimate of the lambda that we use.  #nobs here is the number of _sequences_ -- the fact that each observation involves multiple positions is accounted for within the function by dividing l by nb.
    days.per.generation <- (1.5*2*((phi)/(1+phi))); # the coefficient of lambda/(nucleotide.substitution.rate.per.generation*num.positions) in the PFitter days(..) function.
    generations.per.day <- (1/days.per.generation); # inverse of 1.5*2*((phi)/(1+phi)), which turns out to be 0.5515512.
    days0 <- function( clambda, nb,epsilon) { days( clambda/generations.per.day, nb, epsilon ) }

      # Note about Jeffries prior from https://en.wikipedia.org/wiki/Jeffreys_prior#Poisson_distribution_with_rate_parameter:
    # for the poisson parameter lambda the Jeffries prior is sqrt( 1/lambda ).
    # Equivalently, the Jeffreys prior for {\textstyle {\sqrt {\lambda }}=\int d\lambda /{\sqrt {\lambda }}}{\textstyle {\sqrt {\lambda }}=\int d\lambda /{\sqrt {\lambda }}} is the unnormalized uniform distribution on the non-negative real line.
    
    ## Poisson posterior for the elapsed time since infection (measured in days), given fixed mutation rate given a gamma prior with the given parameters (used as pseudocounts; these can be 0).
    create.poisson.posterior.fn.using.gamma.prior <- function (
        k, # this is the total number of nucleotides in the alignment that did mutate; to center the distribution at x days, set k = x * n * mutation.rate
        n, # note this n is the total number of nucs in the alignment matrix that could have mutated
        mutation.rate = phi*epsilon*generations.per.day/1.5, # Note 1.5 also appears in PF's days(..) fn.. not sure where it comes from!
        prior.gamma.shape = 0,
        prior.gamma.rate = 0,
        min.gamma.shape = base::max( .Machine$double.eps * 100, ( ( n * mutation.rate ) * ( 1 + prior.gamma.rate ) ) ), # centers at 1 day at minimum
        log.results = FALSE,
        be.verbose = TRUE,
        ...
      ) {
        # Basically, if k = 0 (very recent infection?), and if the prior is also 0 (the default), we don't want to risk the degenerate point mass at 0 that the dgamma(..) fn gives if shape == 0. Note that we are not concerned regarding the rate because it is 1 + prior.gamma.rate.
        stopifnot( n > 0 );
        stopifnot( k >= 0 );
        stopifnot( mutation.rate > 0 );
        stopifnot( min.gamma.shape > 0 );
        if( ( k + prior.gamma.shape ) < min.gamma.shape ) {
            if( be.verbose ) {
                print( paste( "WARNING from PFitter's Bayesian code: observed k + prior.gamma.shape is too low in create.poisson.posterior.fn.using.gamma.prior; setting k + prior.gamma.shape to ", min.gamma.shape, " (so that the expected value of the days distribution is ", ( ( min.gamma.shape / ( 1 + prior.gamma.rate ) ) / ( n * mutation.rate ) ), " days)", sep = "" ) );
            }
            # this is necessarily non-negative due to the condition above (assuming k >= 0)
            k <- min.gamma.shape - prior.gamma.shape;
        }
        # Note that this accepts a vector arg days.
        .fn <- function( days ) {
            p <- dgamma( days*n*mutation.rate, shape = k + prior.gamma.shape, rate = 1 + prior.gamma.rate, log = log.results );
            return( p );
        }
        return( .fn );
    } # create.poisson.posterior.fn.using.gamma.prior (..)

    ## Results computation functions

    ## Reduce x until .fn( x ) is finite, nonzero, and stable (insofar as .fn( x ) == .fn( .mutate.x.fn( .x ) )).
    ## Also ensure that there is stability, that is that .fn( x ) == .fn( .mutate.x.fn( .x ) ).
    ensure.finite.and.nonzero <- function ( x, .fn, .mutate.x.fn = function( .x ) { .x * 0.9 }, be.verbose = FALSE, tol = 1E-10 ) {
        .y <- .fn( x );
        .next.x  <- .mutate.x.fn( x );
        .next.y  <- .fn( .next.x );
              if( !is.na( as.logical( all.equal( .y, 0, tol = tol ) ) ) || !is.finite( .y ) || !is.na( as.logical( all.equal( .y, .next.y, tol = tol ) ) ) ) {
                  while( !is.na( as.logical( all.equal( .y, 0, tol = tol ) ) ) || !is.finite( .y ) || !is.na( as.logical( all.equal( .y, .next.y, tol = tol ) ) ) ) {
                      if( be.verbose ) {
                          if( !is.na( as.logical( all.equal( .y, 0 ) ) ) ) {
                              print( paste( x, "does not yield nonzero .fn( x )"  ) );
                          } else if( !is.finite( .y ) ) {
                              print( paste( x, "does not yield finite .fn( x )"  ) );
                          } else {
                              print( paste( x, "does not satisfy .fn( x ) == .fn( .mutate.x.fn( .x ) )"  ) );
                          }
                      } # End if be.verbose
                      if( .next.x == 0 ) {
                          stop( "unable to find an x that yields finite and nonzero .fn( x )" );
                      }
                      x  <- .next.x;
                      .y <- .next.y;
                      .next.x  <- .mutate.x.fn( x );
                      .next.y  <- .fn( .next.x );
                  }
              }
      return( x );
    } # ensure.finite.and.nonzero ( x, .fn )


  compute.posterior.maximum <- 
      function( k, n, posterior.fn.creator.fn, prior.lower.bound = 0, prior.upper.bound = 365*10, be.verbose = FALSE, sample.name = "unnamed sample", tolerance = 0.0001, ... ) {
          if( be.verbose ) {
              cat( paste( "computing posterior maximum for", sample.name ), fill = TRUE );
          }
          .posterior.fn <-
              posterior.fn.creator.fn( k, n, ..., log.results = TRUE, be.verbose = be.verbose );
          .upper.bound <- ensure.finite.and.nonzero( prior.upper.bound, .posterior.fn, be.verbose = be.verbose );
          optimize( .posterior.fn, c( prior.lower.bound, .upper.bound ), maximum = TRUE, tol = tolerance )$maximum
      } # compute.posterior.maximum (..)

  ## We use these to scale things before unlogging them to ensure sufficient numerical precision for the computation.
  get.maximum.log.unnormalized.density <- function ( k, n, posterior.maximum, posterior.fn.creator.fn, be.verbose = FALSE, ... ) {
      .fn <- posterior.fn.creator.fn( k, n, log.results = TRUE, be.verbose = be.verbose, ... );
      return( .fn( posterior.maximum ) );
  } # get.maximum.log.unnormalized.density (..)
  
  create.scaled.unlogged.posterior.fn <- function ( k, n, posterior.maximum, posterior.fn.creator.fn, be.verbose = FALSE, ... ) {
      .fn <- posterior.fn.creator.fn( k, n, log.results = TRUE, be.verbose = be.verbose, ... );
    .scalar <- .fn( posterior.maximum );
    function ( x ) { exp( .fn( x ) - .scalar ) }
  }
    
  ## Now compute the quantiles.  First we must normalize the function.
  compute.total.mass <- 
      function( k, n, posterior.fn.creator.fn, posterior.maximum = posterior.maximum, prior.lower.bound = 0, prior.upper.bound = 365*10, be.verbose = FALSE, sample.name = "unnamed sample", tolerance = 0.0001, ... ) {
          if( be.verbose ) {
              cat( paste( "Integrating posterior density for", sample.name ), fill = TRUE );
          }
        .fn <- create.scaled.unlogged.posterior.fn( k = k, n = n, posterior.maximum = posterior.maximum, posterior.fn.creator.fn = posterior.fn.creator.fn, be.verbose = be.verbose, ... );
         # Sometimes the precision is so poor you need to look around a bit for a good place to compute the total mass.
        #.powers.of.10 <- 10^(0:floor( log10( prior.upper.bound )));
        #max( sapply( .powers.of.10, function( .power.of.10 ) { as.numeric( integrate( .fn, prior.lower.bound, prior.upper.bound / .power.of.10, rel.tol = tolerance )$value ) } ) );
        as.numeric( integrate( .fn, prior.lower.bound, ensure.finite.and.nonzero( prior.upper.bound, .fn, be.verbose = be.verbose ), rel.tol = tolerance )$value )
      } # compute.total.mass (..)
  
  create.normalized.posterior.density <-
    function ( k, n, posterior.fn.creator.fn, posterior.maximum, ... ) {
        .scaled.unlogged.posterior.fn <- create.scaled.unlogged.posterior.fn( k = k, n = n, posterior.maximum = posterior.maximum, posterior.fn.creator.fn = posterior.fn.creator.fn, ... );
      total.mass <-
          compute.total.mass( k, n, posterior.maximum = posterior.maximum, posterior.fn.creator.fn = posterior.fn.creator.fn, ... );
      stopifnot( total.mass > 0 );
      function( x ) { .scaled.unlogged.posterior.fn( x ) / total.mass }
    } # create.normalized.posterior.density (..)

  ## Now make the versions specific to each analysis type
  create.poisson.normalized.posterior.density <-
      function ( k, n, posterior.fn.creator.fn = create.poisson.posterior.fn.using.gamma.prior, prior.lower.bound = 0, prior.upper.bound = 365*10, be.verbose = FALSE, sample.name = "unnamed sample", tolerance = 0.0001, ... ) {
          posterior.maximum <-
              compute.posterior.maximum( k, n, posterior.fn.creator.fn = create.poisson.posterior.fn.using.gamma.prior, prior.lower.bound = prior.lower.bound, prior.upper.bound = prior.upper.bound, be.verbose = be.verbose, sample.name = sample.name, tolerance = tolerance, ... );
          create.normalized.posterior.density( k, n, posterior.fn.creator.fn = create.poisson.posterior.fn.using.gamma.prior, posterior.maximum = posterior.maximum, prior.lower.bound = prior.lower.bound, prior.upper.bound = prior.upper.bound, be.verbose = be.verbose, sample.name = sample.name, tolerance = tolerance, ... )
    } # create.poisson.normalized.posterior.density (..)

  ## Now make the versions specific to each analysis type
     # Ok, so now we need to determine specific quantiles.  I could do this using the memoizing integration table I developed for the MBS code, that hacks it right into integrate(), but instead I'll do the following.
     ## Given a list.of.desired.quantiles, I can start by dividing the entire range from prior.lower.bound to prior.upper.bound into K (eg, 10) partitions.  I keep a list of length K+1 with a tuple at each position: (x, density.at.x, cdf.at.x, desired.quantiles.in.bin ), where desired.quantiles.in.bin are the results, named by the quantile, corresponding to those of the all.desired.quantiles that are below the "cdf.at.x" at the _previous_ position, but that are _not_ below the "cdf.at.x" at this position.  Thus eg c( "0.5" = 10.234, "0.9" = 103.201 ).  "x" marks the top of the bin, and the first x to be evaluated is 0, then 1/K, .. up to 1.  For x=prior.lower.bound the tuple has no desired.quantiles.in.bin unless the function has non-neglible density at this point; then it is all of the desired quantiles that fall below that density, which is also the cdf value there.  All but the first bin can be further broken down by repeating the process for any bin for which desired.quantiles.in.bin is non-empty. We declare to have found an individual quantile whenever abs( cdf.at.x - desired.quantile ) < tolerance, and we continue the recursion until all quantiles are found.

     ## To speed this up, if the lower limit is below 0.5, we just call it zero and be done quickly.  A lot of time used to be spent carefully resolving the lower limit and this is wasted since we effectively round at the end anyway.
     compute.desired.quantiles.results.recursively <-
         function ( k, n, density.fn.creator.fn = create.poisson.normalized.posterior.density, num.partitions = 10, all.desired.quantiles = c( 0.025, 0.5, 0.975 ), prior.lower.bound = 0, prior.upper.bound = 365*10, be.verbose = FALSE, sample.name = "unnamed sample", tolerance = 0.0001, ... ) {
        #print( num.partitions );
        if( be.verbose ) {
            cat( "computing quantiles results for", sample.name, fill = TRUE );
        }
        recursively.compute.desired.quantiles <- function ( density.fn, lower = prior.lower.bound, upper = ensure.finite.and.nonzero( prior.upper.bound, density.fn ), cdf.at.lower.bound = with( list2env( list( ...fn = density.fn, ...lower = lower ) ), as.numeric( integrate( ...fn, lower=0, upper=...lower, rel.tol = tolerance )$value ) ), desired.quantiles = all.desired.quantiles, include.lower.point.partition = ( lower == prior.lower.bound ) ) {
            if( be.verbose ) {
                cat( "recursively.compute.desired.quantiles( lower=", lower, ", upper=", upper, ", cdf.at.lower.bound=", cdf.at.lower.bound, ", desired.quantiles=c(", paste( desired.quantiles, collapse = "," ), "), include.lower.point.partition=", include.lower.point.partition, sep = "", fill = TRUE )
            }
            # Start at zero if include.lower.point.partition; or at 1 otherwise.
            partition.number <- as.numeric( !include.lower.point.partition );
            return.quantiles <- c();
            prev.x <- lower;
            cdf.at.prev.x <- cdf.at.lower.bound;
            while( ( partition.number <= num.partitions ) && ( length( desired.quantiles ) > 0 ) ) {
              recursively.compute.desired.quantiles.helper.results.list <-
                recursively.compute.desired.quantiles.helper(
                    density.fn,
                    partition.number,
                    lower, upper, prev.x, cdf.at.prev.x, desired.quantiles
                );
              prev.x <-
                  recursively.compute.desired.quantiles.helper.results.list$x;
              cdf.at.prev.x <-
                  recursively.compute.desired.quantiles.helper.results.list$cdf.at.x;
              quantile.results.in.bin <-
                  recursively.compute.desired.quantiles.helper.results.list$quantile.results.in.bin;
            if( be.verbose ) {
                cat( "after calling recursively.compute.desired.quantiles.helper, got: ( prev.x=", prev.x, ", cdf.at.prev.x=", cdf.at.prev.x, ", quantile.results.in.bin=c(", paste( quantile.results.in.bin, collapse = "," ), "), names( quantile.results.in.bin )=c(", paste( names( quantile.results.in.bin ), collapse = "," ), ")", sep = "", fill = TRUE )
            }
              if( length( quantile.results.in.bin ) > 0 ) {
                  stopifnot( all( as.numeric( names( quantile.results.in.bin ) ) %in% desired.quantiles ) );
                  return.quantiles <-
                      c( return.quantiles, quantile.results.in.bin );
                  desired.quantiles <-
                      setdiff( desired.quantiles, names( quantile.results.in.bin ) );
              }
              partition.number <- partition.number + 1;
            } # End while scanning partitions
            if( length( desired.quantiles ) > 0 ) {
                warning( paste( "Could not find quantiles: ", paste( desired.quantiles, collapse = ", " ), " in range ", lower, " to ", upper, ". Using ", upper, ".", sep = "" ) );
                .new.quantile.results <- rep( upper, length( desired.quantiles ) );
                names( .new.quantile.results ) <- desired.quantiles;
                quantile.results.in.bin <- c( quantile.results.in.bin, .new.quantile.results );
                return.quantiles <- c( return.quantiles, .new.quantile.results );
                desired.quantiles <- c();
            }
            if( be.verbose ) {
                cat( "after calling recursively.compute.desired.quantiles, got: return.quantiles=c(", paste( return.quantiles, collapse = "," ), "), names( return.quantiles )=c(", paste( names( return.quantiles ), collapse = "," ), ")", sep = "", fill = TRUE )
            }
            return( return.quantiles );
        } # recursively.compute.desired.quantiles (..)
        recursively.compute.desired.quantiles.helper <-
          function ( density.fn, partition.number, lower, upper, prev.x, cdf.at.prev.x, desired.quantiles ) {
            if( be.verbose ) {
                cat( "recursively.compute.desired.quantiles.helper( partition.number=", partition.number, ", lower=", lower, ", upper=", upper, ", prev.x=", prev.x, ", cdf.at.prev.x=", cdf.at.prev.x, ", desired.quantiles=c(", paste( desired.quantiles, collapse = "," ), ") )", sep = "", fill = TRUE )
            }
            if( ( partition.number > 0 ) && any( desired.quantiles <= cdf.at.prev.x ) ) {
                stop( paste( "Could not find quantiles: ", paste( desired.quantiles[ desired.quantiles <= cdf.at.prev.x ], collapse = ", " ), " in range ", lower, " to ", upper, ".", sep = "" ) );
            }
            
            x <- lower + partition.number * ( ( upper - lower ) / num.partitions );
            if( partition.number == 0 ) {
                # Special case: partition number is zero.  This is for the very first partition evaluated, which is just needed in case there is mass at the prior.lower.bound and there are desired.quantiles falling in that region.
                stopifnot( lower == prev.x );
                cdf.within.bin <- 0;
                # This can happen when the integral is a poor approximation at the upper bound so that it's actually larger using a smaller upper bound.  So don't "stopifnot"
                #stopifnot( cdf.within.bin <= 1.0 + 1.01 );
                cdf.at.x <- cdf.at.prev.x; # Don't add the density at x because it's already added.
            } else {
                cdf.within.bin <- as.numeric( integrate( density.fn, prev.x, x, rel.tol = tolerance )$value );
                # This can happen when the integral is a poor approximation at the upper bound so that it's actually larger using a smaller upper bound.  So don't "stopifnot"
                # if( cdf.within.bin > 1.0+0.01 ) {
                #     ## TODO: REMOVE
                #     print( paste("CDF WITHIN BIN:",cdf.within.bin));
                # }
                # stopifnot( cdf.within.bin <= 1.0+0.01 );
                cdf.at.x <- cdf.at.prev.x + cdf.within.bin;
            }
            desired.quantiles.in.bin <-
                desired.quantiles[ ( desired.quantiles <= ( cdf.at.x + tolerance ) ) ];
            density.at.x <- density.fn( x );
            if( length( desired.quantiles.in.bin ) == 0 ) {
              return( list( x = x, density.at.x = density.at.x, cdf.at.x = cdf.at.x, quantile.results.in.bin = c() ) );
            }
            if( partition.number == 0 ) {
                # In the zeroth partition we do not recurse, we just report them at x.
                quantile.results.in.bin <-
                    rep( x, length( desired.quantiles.in.bin ) );
                names( quantile.results.in.bin ) <- desired.quantiles.in.bin;
            } else if( round( prev.x ) == round( x ) ) {
              # Then don't bother getting more precise, we just need it to the nearest day.
              quantile.results.in.bin <- rep( round( x ), length( desired.quantiles.in.bin ) );
              names( quantile.results.in.bin ) <- desired.quantiles.in.bin;
            } else {
                # RECURSE, build up quantile.results.in.bin.
                # First for the base case: check the desired quantiles -- are any close enough?
                desired.quantiles.within.tolerance <-
                    desired.quantiles.in.bin[ abs( desired.quantiles.in.bin - cdf.at.x ) < tolerance ];
                if( length( desired.quantiles.within.tolerance ) > 0 ) {
                    # Then we've found them so they can be removed from the recursion.
                    quantile.results.in.bin <-
                        rep( x, length( desired.quantiles.within.tolerance ) );
                    names( quantile.results.in.bin ) <- desired.quantiles.within.tolerance;
                    desired.quantiles.in.bin <-
                        setdiff( desired.quantiles.in.bin, desired.quantiles.within.tolerance );
                } else {
                    quantile.results.in.bin <- c();
                }
                if( length( desired.quantiles.in.bin ) > 0 ) {
                      .rv <-
                          recursively.compute.desired.quantiles(
                              density.fn = density.fn,
                              lower = prev.x,
                              upper = x,
                              cdf.at.lower.bound = cdf.at.prev.x,
                              desired.quantiles = desired.quantiles.in.bin,
                              include.lower.point.partition = FALSE
                          );
                      # At this point desired.quantiles.in.bin is the percentiles like "0.5"; these will become the _names_ of the returned value of desired.quantiles.in.bin.
                      quantile.results.in.bin <-
                          c( quantile.results.in.bin,
                             .rv );
                  }
              }
            return( list( x = x, density.at.x = density.at.x, cdf.at.x = cdf.at.x, quantile.results.in.bin = quantile.results.in.bin ) );
          } # recursively.compute.desired.quantiles.helper (..)
        ## NOTE That we do not pass the prior bounds or tolerance to the underlying density fn.
        density.fn <-
            density.fn.creator.fn( k, n, prior.lower.bound = 0, prior.upper.bound = 365*10, be.verbose = FALSE, sample.name = "unnamed sample", tolerance = 0.0001, ... );
        desired.quantiles <-
            recursively.compute.desired.quantiles( density.fn );
        return( desired.quantiles );
    } # compute.desired.quantiles.results.recursively (..)
compute.poisson.desired.quantiles.results.recursively <-
    compute.desired.quantiles.results.recursively;

compute.posterior.probability.by.day <- 
    function ( k, n, from = 0, to = 365*2, density.fn.creator.fn = create.poisson.normalized.posterior.density, prior.lower.bound = 0, prior.upper.bound = 365*10, be.verbose = FALSE, sample.name = "unnamed sample", tolerance = 0.0001, ... ) {
        #print( num.partitions );
l        if( be.verbose ) {
            cat( "computing cumulative posterior probability by day results for", sample.name, fill = TRUE );
        }
        density.fn <-
            density.fn.creator.fn( k, n, prior.lower.bound = prior.lower.bound, prior.upper.bound = prior.upper.bound, be.verbose = be.verbose, sample.name = sample.name, tolerance = tolerance, ... );
        sapply( from:to, function ( day ) {
            with( list2env( list( ...fn = density.fn, ...lower = day - 0.5, ...upper = ( day + 0.5 ) ) ), as.numeric( integrate( ...fn, lower=...lower, upper=...upper, rel.tol = tolerance )$value ) )
        } );
} # compute.posterior.probability.by.day (..)

