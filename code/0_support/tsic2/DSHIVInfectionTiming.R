### R code from vignette source '/Users/pedlefsen/src/from-git/projects/pedlefse/ds/DSHIVInfectionTiming.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: DSHIVInfectionTiming.Rnw:84-88
###################################################
## R packages needed

# Setup for prettier Sweave output.
old.continue.option <- options( continue = " " )


###################################################
### code chunk number 2: DSHIVInfectionTiming.Rnw:368-373
###################################################


ten.million.samples.from.eclipse.phase.weibull <- 4.8 + rweibull( 10000000, shape = 1.35, scale = 9 );
# sampled density curve plotted
hist( ten.million.samples.from.eclipse.phase.weibull, breaks = 500, freq = FALSE, main = "10M draws from weibull( location = 4.8, shape = 1.35, scale = 9 )", xlab = "Days", xlim = c( 0, 50 ) )


###################################################
### code chunk number 3: DSHIVInfectionTiming.Rnw:380-381
###################################################
plot( c( 4.8 + ( 1:50000 )/500 ), c( dweibull( c( ( 1:50000 )/500 ), shape = 1.35, scale = 9 ) ), xlim = c( 0, 50 ), pch = 20 , main = "Density of weibull( location = 4.8, shape = 1.35, scale = 9 )", xlab = "Day", ylab = "Density" )


###################################################
### code chunk number 4: DSHIVInfectionTiming.Rnw:389-390
###################################################
plot( c( ((0:48)/10), 4.8 + qweibull( c( ( 1:5000 )/5000 ), shape = 1.35, scale = 9 ) ), c( rep( 0, 49 ), ( 1:5000 )/5000 ), xlim = c( 0, 50 ), ylim = c( 0, 1.2 ), pch = 20, main = "CDF of weibull( location = 4.8, shape = 1.35, scale = 9 )", xlab = "Day", ylab = "Cumulative Probability" )


###################################################
### code chunk number 5: DSHIVInfectionTiming.Rnw:405-514
###################################################
## "Dynamics" data type was defined by Phillip Labuschagne for "Assay Dynamics" in his implementation of the IDT methodology, and I have used it to try to keep modification requirements down.

## "additional.location.shift" is convenient when the assays within a class have just constant delays in the time to positivity, as we are assuming here, see below. It is also convenient for integrating the function over a range of values, as we do. Applying density.fn.returned.with.additional.location.shift( time ) is equivalent to density.fn.returned.with.no.additional.location.shift( ( time - additional.location.shift ) ).
get.dynamics.density.function <- function ( the.dynamics, additional.location.shift = 0, log.results = FALSE ) {
    if( the.dynamics$form != "weib3" ) {
        stop( paste( "Unrecognized form of dynamics:", the.dynamics$form ) );
    }
    location.shift <- the.dynamics$params$location + additional.location.shift;
    return( function( time ) {
        ifelse( time <= location.shift, 0, dweibull( time - location.shift, shape = the.dynamics$params$shape, scale = the.dynamics$params$scale, log = log.results ) )
    } );
} # get.dynamics.density.function (..)

get.dynamics.quantile.function <- function ( the.dynamics, additional.location.shift = 0 ) {
    if( the.dynamics$form != "weib3" ) {
        stop( paste( "Unrecognized form of dynamics:", the.dynamics$form ) );
    }
    location.shift <- the.dynamics$params$location + additional.location.shift;
    return( function( q ) {
        location.shift + qweibull( q, shape = the.dynamics$params$shape, scale = the.dynamics$params$scale )
    } );
} # get.dynamics.quantile.function (..)

get.dynamics.random.draw.function <- function ( the.dynamics, additional.location.shift = 0 ) {
    if( the.dynamics$form != "weib3" ) {
        stop( paste( "Unrecognized form of dynamics:", the.dynamics$form ) );
    }
    location.shift <- the.dynamics$params$location + additional.location.shift;
    return( function( n = 1 ) {
        location.shift + rweibull( n, shape = the.dynamics$params$shape, scale = the.dynamics$params$scale )
    } );
} # get.dynamics.random.draw.function (..)

get.dynamics.cumulative.distribution.function <- function ( the.dynamics, additional.location.shift = 0, log.results = FALSE ) {
    if( the.dynamics$form != "weib3" ) {
        stop( paste( "Unrecognized form of dynamics:", the.dynamics$form ) );
    }
    location.shift <- the.dynamics$params$location + additional.location.shift;
    return( function( time ) {
        ifelse( time <= location.shift, 0, pweibull( time - location.shift, shape = the.dynamics$params$shape, scale = the.dynamics$params$scale, log = log.results ) )
    } );
} # get.dynamics.cumulative.distribution.function (..)

get.dynamics.minimum <- function ( the.dynamics, additional.location.shift = 0, log.results = FALSE ) {
    if( the.dynamics$form != "weib3" ) {
        stop( paste( "Unrecognized form of dynamics:", the.dynamics$form ) );
    }
    location.shift <- the.dynamics$params$location + additional.location.shift;
    return( location.shift );
} # get.dynamics.minimum (..)

# Figure out where the cdf is essentially 1. Arbitrarily, double the value where the quantile function becomes the closest value indistinguishable from 1 on my mac (64 bit) laptop. I fdigured out that 1-E-7 resulves to 1-1E-7 = 0.9999999, but 1-1E-8 is 1 on my machine, so I'm using 1E-7 as the smallest value subtractable from 1 and still 1 - then I look that up in the quantile function, and return twice that value.
get.dynamics.maximum <- function ( the.dynamics, additional.location.shift = 0, log.results = FALSE, arbitrary.scale.factor = 2, largest.number.less.than.1 = 1-1E-7 ) {
    arbitrary.scale.factor * get.dynamics.quantile.function( the.dynamics )( largest.number.less.than.1 )
} # get.dynamics.maximum (..)

## These are the functions that we called Tau in the text, describing the eclipse phase distribution and sequential ITRIs.
## These are made up for now, until we can get better data from eg Kevin Delaney et al.
biomarker_dynamics <- list(
    "RNA" = list(
      class = 'RNA',
      full_assayname = 'RNA viral load at least 1 copy per ml',
      short_assayname = 'RNA.1',
      form = 'weib3',
      source = "delaney_2017",
      fun = 'weib3_assay_dynamics',
      params = list(location = 4.8, shape = 1.350, scale = 9) ),
   "Ag/Ab" = list(
      class = 'Ag/Ab',
      full_assayname = 'p24 at least 1 unit per ml',
      short_assayname = 'p24.1',
      form = 'weib3',
      source = "delaney_2017",
      fun = 'weib3_assay_dynamics',
      params = list(location = -7, shape = 4, scale = 14) ), # This works, approximately
    "IgG_Supp" = list(
      class = 'IgG_Supp',
      full_assayname = 'IgG at least 1 unit per ml',
      short_assayname = 'Ab.1',
      form = 'weib3',
      source = "delaney_2017",
      fun = 'weib3_assay_dynamics',
      params = list(location = 10, shape = 1.733, scale = 6) ) # This works, approximately
); # biomarker_dynamics list

## within a class measuring a common biomarker eg RNA, the various assays are here represented as essentially the same, just perhaps shifted using a linear offset (modifying the "location" parameter of the weib3 model) to ensure the median matches. For now we are assuming that all of the "\alpha" functions above are simple constant shifts to reflect only constant linear delay within a class. This could be modified in future, but see above where we take advantage of it in simplifying the code.
get.time.from.biomarker.crossing.threshold.to.assay.positivity <- function ( the.assay.dynamics, the.biomarker.dynamics = biomarker_dynamics[[ the.assay.dynamics$class ]] ) {
    # Shift by the difference in medians; but the biomarker median is the median sum over each biomarker class in between.
    assay.crossing.threshold.median.time <-
        get.dynamics.quantile.function( the.assay.dynamics )( .5 ); # this is a window period distribution
    biomarker.crossing.threshold.median.time.since.preceding.event <-
        get.dynamics.quantile.function( the.biomarker.dynamics )( .5 ); # this is an ITRI distribution
    ## Get all of the preceding events, add their medians too.
    if( the.biomarker.dynamics$class == names( biomarker_dynamics )[ 1 ] ) {
        # If this is the fastest class, don't try looking for faster.
        biomarker.crossing.threshold.median.time <- 
            biomarker.crossing.threshold.median.time.since.preceding.event;
    } else {
        biomarker.crossing.threshold.median.time <- 
            biomarker.crossing.threshold.median.time.since.preceding.event +
                sum(
                    sapply( 1:( which( the.biomarker.dynamics$class == names( biomarker_dynamics ) ) - 1 ), function( .which.biomarker ) {
                        get.dynamics.quantile.function( biomarker_dynamics[[ .which.biomarker ]] )( .5 ) # this is an ITRI distribution
                    } )
            );
    } # End if this is the fastest class, don't try looking for faster .. else ..
    return( assay.crossing.threshold.median.time - biomarker.crossing.threshold.median.time );
} # get.time.from.biomarker.crossing.threshold.to.assay.positivity (..)



###################################################
### code chunk number 6: DSHIVInfectionTiming.Rnw:517-700
###################################################

## "additional.location.shift" is convenient for integrating the function over a range of values. Applying density.fn.returned.with.additional.location.shift( time ) is equivalent to density.fn.returned.with.no.additional.location.shift( ( time - additional.location.shift ) ).
density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold <- function ( time, additional.location.shift = 0, log.results = FALSE, biomarker.dynamics.class = "RNA", the.biomarker.dynamics = biomarker_dynamics[[ biomarker.dynamics.class ]] ) {
    get.dynamics.density.function( the.biomarker.dynamics, additional.location.shift = additional.location.shift, log.results = log.results )( time )
} # density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold (..)

## "additional.location.shift" is convenient for integrating the function over a range of values. Applying density.fn.returned.with.additional.location.shift( time ) is equivalent to density.fn.returned.with.no.additional.location.shift( ( time - additional.location.shift ) ).
integrated.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold <- function ( time.low, time.high, additional.location.shift = 0, log.results = FALSE, biomarker.dynamics.class = "RNA", the.biomarker.dynamics = biomarker_dynamics[[ biomarker.dynamics.class ]] ) {
    the.cdf.fn <- get.dynamics.cumulative.distribution.function( the.biomarker.dynamics, additional.location.shift = additional.location.shift, log.results = log.results )
    if( log.results ) {
        log( exp( the.cdf.fn( time.high ) ) - exp( the.cdf.fn( time.low ) ) )
    } else {
        the.cdf.fn( time.high ) - the.cdf.fn( time.low )
    }
} # integrated.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold (..)

## We have approximated the dynamics fairly well.
## dynamics for geenious:
.w <- rweibull( 100000, shape = 1.733, scale = 14.483) 
summary( get.dynamics.random.draw.function( biomarker_dynamics[[ 3 ]] )( 10000 ) + get.dynamics.random.draw.function( biomarker_dynamics[[ 2 ]] )( 10000 ) + get.dynamics.random.draw.function( biomarker_dynamics[[ 1 ]] )( 10000 ) )
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   14.34   28.46   33.17   34.06   38.63   87.05 
summary(21.151+.w)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   21.16   28.17   32.86   34.05   38.64   83.34 
var(21.151+.w)
# [1] 59.12597
var( get.dynamics.random.draw.function( biomarker_dynamics[[ 3 ]] )( 10000 ) + get.dynamics.random.draw.function( biomarker_dynamics[[ 2 ]] )( 10000 ) + get.dynamics.random.draw.function( biomarker_dynamics[[ 1 ]] )( 10000 ) )
# [1] 59.59556

## dynamics for Ag/Ab combo:
.z <- rweibull( 100000, shape = 1.725, scale = 13.185) 
summary( get.dynamics.random.draw.function( biomarker_dynamics[[ 2 ]] )( 10000 ) + get.dynamics.random.draw.function( biomarker_dynamics[[ 1 ]] )( 10000 ) )
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.09114 13.75140 17.82177 18.72593 22.68128 60.80533 
summary( 7.209+.z )
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   7.223  13.596  17.852  18.944  23.133  59.620 
var( get.dynamics.random.draw.function( biomarker_dynamics[[ 2 ]] )( 10000 ) + get.dynamics.random.draw.function( biomarker_dynamics[[ 1 ]] )( 10000 ) )
# [1] 50.29131
var( 7.209+.z )
# [1] 49.39667

#get.time.from.biomarker.crossing.threshold.to.assay.positivity( all_assay_dynamics[[1]] )
#[1] 0.8282798
#> get.time.from.biomarker.crossing.threshold.to.assay.positivity( all_assay_dynamics[[2]] )
#[1] 0.435775
#> get.time.from.biomarker.crossing.threshold.to.assay.positivity( all_assay_dynamics[[3]] )
#[1] 0.5825461

# pdf( "RNA_positive_curve.pdf" )
# plot( c( 0, 4.8, 4.8 + ( 1:50000 )/500 ), density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold( c( 0, 4.8, c( ( 1:50000 )/500 ) ) ), xlim = c( 0, 50 ), lwd = 10, type = "l", main = "Time until RNA crosses positivity threshold", xlab = "Day since infection-causing exposure", ylab = "Density" )
# dev.off() 
# 
# pdf( "RNA_positive_curve_reflected.pdf" )
# plot( c( 0, -4.8, -4.8 - ( 1:50000 )/500 ), density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold( c( 0, 4.8, c( ( 1:50000 )/500 ) ) ), lwd = 10, type = "l", xlim = c( -50, 0 ), main = "Posterior density of time of infection prior to RNA crossing day given a uniform prior", xlab = "Day", ylab = "Posterior Density" )
# dev.off() 
# 
# pdf( "RNA_positive_curve_reflected_log.pdf" )
# plot( c( 0, -4.8, -4.8 - ( 1:50000 )/500 ), density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold( c( 0, 4.8, c( ( 1:50000 )/500 ) ), log.results = TRUE ), lwd = 10, type = "l", xlim = c( -50, 0 ), main = "Log posterior density of time of infection prior to RNA crossing day given a uniform prior", xlab = "Day", ylab = "Log Posterior Density" )
# dev.off() 

## The gap between negative and positive tests corresponds to "don't know" uncertainty about the time when the threshold was crossed. With a uniform or other prior, let's say uniform, it becomes diluted; that is the marginal posterior is an even mixture of the posteriors if each moment were when the threshold were crossed, so each moment gets marginal posterior proportional to sum over window of 1 day to X days (the days between the adjacent tests).

# Here we consider that if there is a range of uncertainty about the biomarker crossing time, the posterior is a mixture over the prior of the possible times of that crossing, which effectively shift the distribution so eg if the threshold was crossed 2 days before "time", we need to evaluate the time since infection function at time - 2 days. You specify this by calling density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold(..) with the additional argument additional.location.shift = 2 (note the additional shift is positive in this case). NOTE: if log.results is TRUE, integration happens on natural scale, result is logged.
unnormalized.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold <- function ( time, shift.low = 0, shift.high = 1, biomarker.dynamics.class = "RNA", the.biomarker.dynamics = biomarker_dynamics[[ biomarker.dynamics.class ]], min.time = get.dynamics.minimum( the.biomarker.dynamics ), max.time = get.dynamics.maximum( the.biomarker.dynamics ), prior.over.shift = function ( shift ) { 1 }, log.results = FALSE, abs.tol = 0.0001, rel.tol = 0.0001, subdivisions = 100L ) {
    if( ( time <= min.time ) || ( time >= max.time ) ) { 
        return( ifelse( log.results, -Inf, 0 ) );
    }
    # when ( time - shift ) < min.time, prob is 0.
    if( ( shift.low > 0 ) && ( ( time - shift.low ) <= min.time ) ) {
        return( ifelse( log.results, -Inf, 0 ) );
    }
    # replace if necessary
    shift.high.max <- base::max( 0, time - min.time );
    if( shift.high > shift.high.max ) {
        shift.high <- shift.high.max;
    }
    if( is.null( prior.over.shift ) ) {
        if( TRUE || be.verbose ) { ## TODO: PUT BACK / REMOVE TRUE ||
            print( "Using a uniform prior" );
        }
        if( shift.low == shift.high ) {
            # Special case: if this is zero-width, return the density (not zero, since this is used for the posterior density).
            the.density.fn.as.fn.of.shift <- 
                function( shift ) { density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold( time, additional.location.shift = shift, log.results = FALSE, the.biomarker.dynamics = the.biomarker.dynamics ) };
            .x <- the.density.fn.as.fn.of.shift( shift.low );
        } else {
            .x <- integrated.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold( time - shift.high, time - shift.low, additional.location.shift = 0, log.results = FALSE, the.biomarker.dynamics = the.biomarker.dynamics );
        }
    } else {
        .fntointegrate <- function( shift ) { prior.over.shift( shift ) * density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold( time, additional.location.shift = shift, log.results = FALSE, the.biomarker.dynamics = the.biomarker.dynamics ) }
        if( shift.low == shift.high ) {
            # Special case: if this is zero-width, return the density (not zero, since this is used for the posterior density).
            .x <- .fntointegrate( shift.low );
        } else {
            .x <- integrate( .fntointegrate, lower = shift.low, upper = shift.high, rel.tol = rel.tol, abs.tol = abs.tol, subdivisions = subdivisions )[[1]];
        }
    } # End if is.null( prior.over.shift ) .. else ..
    return( ifelse( log.results, log( .x ), .x ) );
} # unnormalized.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold (..)

# unnormalized.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold <- function ( time, shift.low = 0, shift.high = 1, biomarker.dynamics.class = "RNA", the.biomarker.dynamics = biomarker_dynamics[[ biomarker.dynamics.class ]], min.time = get.dynamics.minimum( the.biomarker.dynamics ), max.time = get.dynamics.maximum( the.biomarker.dynamics ), prior.over.shift = function ( shift ) { 1 }, log.results = FALSE, abs.tol = 0.0001, rel.tol = 0.0001, subdivisions = 100L ) {
#     if( ( time <= min.time ) || ( time >= max.time ) ) { 
#         return( ifelse( log.results, -Inf, 0 ) );
#     }
#     # when ( time - shift ) < min.time, prob is 0.
#     if( ( shift.low > 0 ) && ( ( time - shift.low ) <= min.time ) ) {
#         return( ifelse( log.results, -Inf, 0 ) );
#     }
#     # replace if necessary
#     shift.high.max <- base::max( 0, time - min.time );
#      if( shift.high > shift.high.max ) {
#          shift.high <- shift.high.max;
#      }
#     .fntointegrate <- function( shift ) { prior.over.shift( shift ) * density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold( time, additional.location.shift = shift, log.results = FALSE, the.biomarker.dynamics = the.biomarker.dynamics ) }
#         # Usually 0.
#         .lower <- shift.low;
#         # Shift can't exceed time
#         .upper <- shift.high;
#     if( .lower == .upper ) {
#         # Special case: if this is zero-width, return the density (not zero, since this is used for the posterior density).
#         .x <- .fntointegrate( .lower );
#     } else {
#         .x <- integrate( .fntointegrate, lower = .lower, upper = .upper, abs.tol = abs.tol )[[1]];
#     }
#     return( ifelse( log.results, log( .x ), .x ) );
# } # unnormalized.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold (..)


# This is the density relative to the time at the right end of the uncertainty interval about when the biomarker crossed the threshold.
create.unnormalized.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold.from.uncertainty.interval <- function( low.boundary, high.boundary, biomarker.dynamics.class = "RNA", the.biomarker.dynamics = biomarker_dynamics[[ biomarker.dynamics.class ]], prior.lower.bound = get.dynamics.minimum( the.biomarker.dynamics ), prior.upper.bound = get.dynamics.maximum( the.biomarker.dynamics ), prior.over.biomarker.crossing.threshold.time = function ( time ) { 1 }, log.results = FALSE, be.verbose = FALSE, sample.name = "unnamed sample", abs.tol = 0.0001, rel.tol = 0.0001, subdivisions = 100L ) {
    if( is.null( prior.over.biomarker.crossing.threshold.time ) ) {
        prior.over.shift <- NULL;
    } else {
        prior.over.shift <- function( shift ) { ifelse( shift < 0, 0, ifelse( shift > ( high.boundary - low.boundary ), 0, ifelse( high.boundary == low.boundary, 1, prior.over.biomarker.crossing.threshold.time( low.boundary + shift ) ) ) ) };
    }
    Vectorize( function( time ) {
        ifelse( time >= high.boundary, 0,
            unnormalized.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold(
                   high.boundary - time,
                   shift.low = 0,
                   shift.high = ( high.boundary - low.boundary ),
                   the.biomarker.dynamics = the.biomarker.dynamics,
                   min.time = low.boundary - prior.upper.bound,
                   max.time = high.boundary - prior.lower.bound,
                   prior.over.shift = prior.over.shift, 
                   log.results = log.results,
                   rel.tol = rel.tol, abs.tol = abs.tol, subdivisions = subdivisions
               )
        )
    } )
} # create.unnormalized.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold.from.uncertainty.interval (..)

# ## This is what our posterior would be if we condition on the time that a particular assay becomes positive.  We don't know this precisely; see below.
# create.unnormalized.posterior.on.time.of.preceding.event.given.diagnostic.assay.positivity.time <- function ( the.assay.dynamics, rna.biomarker.dynamics = biomarker_dynamics[[ the.assay.dynamics$class ]], prior.lower.bound = get.dynamics.minimum( the.assay.dynamics ), prior.upper.bound = get.dynamics.maximum( the.assay.dynamics ), log.results = FALSE, be.verbose = FALSE, abs.tol = 0.0001, rel.tol = 0.0001, subdivisions = 100L ) {
#     additional.location.shift <-
#         get.time.from.biomarker.crossing.threshold.to.assay.positivity( the.assay.dynamics, the.biomarker.dynamics = the.biomarker.dynamics );
#     return( function( time ) { density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold( time, additional.location.shift = additional.location.shift, log.results = log.results, the.biomarker.dynamics = the.biomarker.dynamics ) } );
# } # create.unnormalized.posterior.on.time.of.preceding.event.given.diagnostic.assay.positivity.time (..)


## In case the upper and lower ends of the the biomarker curves use the same assay, use this version of the below two functions.
## This is what our posterior is given an uncertain range of potential times that an assay crossed the threshold. Here we make use of the assumed form of the "\alpha" functions to skip a step. NOTE That if we allowed those distributions to be anything other than a constant shift, we would have to make modifications here.
# create.unnormalized.posterior.on.time.of.preceding.event.from.biomarker.diagnostic.interval <- function( last.negative.time, first.positive.time, the.assay.dynamics, the.biomarker.dynamics = biomarker_dynamics[[ the.assay.dynamics$class ]], prior.lower.bound = get.dynamics.minimum( the.assay.dynamics ), prior.upper.bound = get.dynamics.maximum( the.assay.dynamics ), prior.over.biomarker.crossing.threshold.time = function ( time ) { 1 }, log.results = FALSE, be.verbose = FALSE, sample.name = "unnamed sample", abs.tol = 0.0001, rel.tol = 0.0001, subdivisions = 100L ) {
#     ## This is "alpha", here assumed to be a constant for the assay:
#     constant.delay.for.assay.from.biomarker.class.canonical.event <-
#         get.time.from.biomarker.crossing.threshold.to.assay.positivity( the.assay.dynamics, the.biomarker.dynamics = the.biomarker.dynamics );
#     create.unnormalized.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold.from.uncertainty.interval( low.boundary = last.negative.time - constant.delay.for.assay.from.biomarker.class.canonical.event, high.boundary = first.positive.time - constant.delay.for.assay.from.biomarker.class.canonical.event, the.biomarker.dynamics = the.biomarker.dynamics, prior.lower.bound = prior.lower.bound, prior.upper.bound = prior.upper.bound, prior.over.biomarker.crossing.threshold.time = prior.over.biomarker.crossing.threshold.time, log.results = log.results, be.verbose = be.verbose, sample.name = sample.name, rel.tol = rel.tol, abs.tol = abs.tol, subdivisions = subdivisions )
# } # create.unnormalized.posterior.on.time.of.preceding.event.from.biomarker.diagnostic.interval (..)
# 
# ## In case the upper and lower ends of the the biomarker curves use different assays, use this version of the above function.
# ## This is what our posterior is given an uncertain range of potential times that an assay crossed the threshold. Here we make use of the assumed form of the "\alpha" functions to skip a step. NOTE That if we allowed those distributions to be anything other than a constant shift, we would have to make modifications here.
# create.unnormalized.posterior.on.time.of.preceding.event.from.biomarker.diagnostic.interval.twotests <- function( last.negative.time, last.negative.assay.dynamics, first.positive.time, first.positive.assay.dynamics, biomarker.dynamics.class = last.negative.assay.dynamics$class, the.biomarker.dynamics = biomarker_dynamics[[ biomarker.dynamics.class ]], prior.lower.bound = get.dynamics.minimum( the.assay.dynamics ), prior.upper.bound = get.dynamics.maximum( last.negative.assay.dynamics ), prior.over.biomarker.crossing.threshold.time = function ( time ) { 1 }, log.results = FALSE, be.verbose = FALSE, sample.name = "unnamed sample", abs.tol = 0.0001, rel.tol = 0.0001, subdivisions = 100L ) {
#     ## This is "alpha", here assumed to be a constant but separeate for each assay; they should be in the same biomarker class.
#     stopifnot( biomarker.dynamics.class == first.positive.assay.dynamics$class );
#     stopifnot( biomarker.dynamics.class == last.negative.assay.dynamics$class );
#     last.negative.constant.delay.for.assay.from.biomarker.class.canonical.event <-
#         get.time.from.biomarker.crossing.threshold.to.assay.positivity( last.negative.assay.dynamics, the.biomarker.dynamics = the.biomarker.dynamics );
#     first.positive.constant.delay.for.assay.from.biomarker.class.canonical.event <-
#         get.time.from.biomarker.crossing.threshold.to.assay.positivity( first.positive.assay.dynamics, the.biomarker.dynamics = the.biomarker.dynamics );
#     create.unnormalized.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold.from.uncertainty.interval( low.boundary = last.negative.time - last.negative.constant.delay.for.assay.from.biomarker.class.canonical.event, high.boundary = first.positive.time - first.positive.constant.delay.for.assay.from.biomarker.class.canonical.event, the.biomarker.dynamics = the.biomarker.dynamics, prior.lower.bound = prior.lower.bound, prior.upper.bound = prior.upper.bound, prior.over.biomarker.crossing.threshold.time = prior.over.biomarker.crossing.threshold.time, log.results = log.results, be.verbose = be.verbose, sample.name = sample.name, rel.tol = rel.tol, abs.tol = abs.tol, subdivisions = subdivisions )
# } # create.unnormalized.posterior.on.time.of.preceding.event.from.biomarker.diagnostic.interval.twotests (..)


###################################################
### code chunk number 7: DSHIVInfectionTiming.Rnw:703-851
###################################################

## Some general-purpose functions.
## Results computation functions
ensure.finite.and.nonzero <- function ( x, .fn, .mutate.x.fn = function( .x ) { .x + 1 }, be.verbose = FALSE, .selfcount = 0, .MAX.SELFCOUNT = 1000 ) {
    stopifnot( is.finite( x ) );
    .y <- .fn( x );
    if( ( .y == 0 ) || !is.finite( .y ) ) {
        while( ( .y == 0 ) || !is.finite( .y ) ) {
            .selfcount <- .selfcount + 1;
            stopifnot( .selfcount < .MAX.SELFCOUNT );

            if( !is.finite( x ) || ( x == 0 ) ) {
                stop( "unable to find an x that yields finite and nonzero .fn( x )" );
            }
            if( be.verbose ) {
                if( .y == 0 ) {
                    print( paste( x, "does not yield nonzero .fn( x )"  ) );
                } else {
                    print( paste( x, "does not yield finite .fn( x )"  ) );
                }
            } # End if be.verbose
            x <- .mutate.x.fn( x );
            .y <- .fn( x );
        } # End while( ( .y == 0 ) || !is.finite( .y ) )
    }
    return( x );
} # ensure.finite.and.nonzero ( x, .fn )
ensure.finite.and.nonzero.from.below <- function( x, .fn, ... ) { ensure.finite.and.nonzero( x, .fn, .mutate.x.fn = function( .x ) { .x + 1 }, ... ) }
ensure.finite.and.nonzero.from.above <- function( x, .fn, ... ) { ensure.finite.and.nonzero( x, .fn, .mutate.x.fn = function( .x ) { .x - 1 }, ... ) }

## NOT PRESENTLY USED DIRECTLY, BUT COULD BE..
compute.density.maximum <- function  ( density.fn, lower, upper, be.verbose = FALSE, sample.name = "unnamed sample", abs.tol = 0.0001, rel.tol = 0.0001, subdivisions = 100L, ... ) {
    if( be.verbose ) {
        cat( paste( "Computing density maximum for", sample.name ), fill = TRUE );
    }

    ## Make sure the priors are set to just outside the range, if they are way outside it.
    if( !( ( density.fn( lower ) > 0 ) && ( density.fn( upper ) > 0 ) ) ) {
        if( be.verbose ) {
            cat( paste( "Computing non-zero range for density for", sample.name ), fill = TRUE );
        }
        lower <-
            ensure.finite.and.nonzero.from.below( lower, density.fn, be.verbose = be.verbose );
        upper <-
            ensure.finite.and.nonzero.from.above( upper, density.fn, be.verbose = be.verbose );
    }

    density.maximum <-
        optimize( density.fn, c( lower, upper ), maximum = TRUE, tol = abs.tol )$maximum;
    return( density.maximum );
} # compute.density.maximum( .. )

is.uniform.over.range <- function  ( density.fn, lower, upper, density.maximum = NULL, be.verbose = FALSE, sample.name = "unnamed sample", abs.tol = 0.0001, rel.tol = 0.0001, subdivisions = 100L, ... ) {
    ## Make sure the priors are set to just outside the range, if they are way outside it.
    density.fn.result.lower <- density.fn( lower );
    density.fn.result.upper <- density.fn( upper );
    if( !( ( density.fn.result.lower > 0 ) && ( density.fn.result.upper > 0 ) ) ) {
        if( be.verbose ) {
            cat( paste( "Computing non-zero range for density for", sample.name ), fill = TRUE );
        }
        new.lower <-
            ensure.finite.and.nonzero.from.below( lower, density.fn, be.verbose = be.verbose );
        if( new.lower != lower ) {
            lower <- new.lower;
            density.fn.result.lower <- density.fn( lower );
        }
        new.upper <-
            ensure.finite.and.nonzero.from.above( upper, density.fn, be.verbose = be.verbose );
        if( new.upper != upper ) {
            upper <- new.upper;
            density.fn.result.upper <- density.fn( upper );
        }
    } # End if we need to recompute the bounds.

    if( !isTRUE( all.equal( density.fn.result.lower, density.fn.result.upper, tol = rel.tol ) ) ) {
        return( FALSE );
    }
    
    if( is.null( density.maximum ) || is.na( density.maximum ) ) {
        density.maximum.result <-
            optimize( density.fn, c( lower, upper ), maximum = TRUE, tol = abs.tol );
        density.maximum <- density.maximum.result$maximum;
        density.fn.result.maximum <- density.maximum.result$objective;
    }
    return( isTRUE( all.equal( density.fn.result.lower, density.fn.result.maximum, tol = rel.tol ) ) && isTRUE( all.equal( density.fn.result.upper, density.fn.result.maximum, tol = rel.tol ) ) );
} # is.uniform.over.range( .. )

# set posterior.maximum to NULL to compute it.
create.normalized.posterior.density <- function ( posterior.fn.creator.fn, prior.lower.bound, prior.upper.bound, posterior.maximum = NULL, be.verbose = FALSE, sample.name = "unnamed sample", abs.tol = 0.0001, rel.tol = 0.0001, subdivisions = 100L, add.numeric.stability.at.the.sacrifice.of.speed = TRUE, total.mass.rel.tol = 100*rel.tol, ... ) {
    ## Make sure the priors are set to just outside the range, if they are way outside it.
    .posterior.fn <-
        posterior.fn.creator.fn( log.results = FALSE, be.verbose = be.verbose, rel.tol = rel.tol, abs.tol = abs.tol, subdivisions = subdivisions, prior.lower.bound = prior.lower.bound, prior.upper.bound = prior.upper.bound, ... );
    .lower <-
        ensure.finite.and.nonzero.from.below( prior.lower.bound, .posterior.fn, be.verbose = be.verbose );
    .upper <-
        ensure.finite.and.nonzero.from.above( prior.upper.bound, .posterior.fn, be.verbose = be.verbose );
    prior.lower.bound <- .lower - 1;
    prior.upper.bound <- .upper + 1;

    .posterior.fn <-
        posterior.fn.creator.fn( log.results = FALSE, be.verbose = be.verbose, rel.tol = rel.tol, abs.tol = abs.tol, subdivisions = subdivisions, prior.lower.bound = prior.lower.bound, prior.upper.bound = prior.upper.bound, ... );
    
    if( add.numeric.stability.at.the.sacrifice.of.speed ) {
        if( be.verbose ) {
            cat( paste( "Computing log posterior maximum for", sample.name ), fill = TRUE );
        }
        if( is.null( posterior.maximum ) ) {
            posterior.maximum <-
                optimize( .posterior.fn, c( .lower, .upper ), maximum = TRUE, tol = abs.tol )$maximum;
        }
        .log.posterior.fn <-
            posterior.fn.creator.fn( log.results = TRUE, be.verbose = be.verbose, rel.tol = rel.tol, abs.tol = abs.tol, subdivisions = subdivisions, prior.lower.bound = prior.lower.bound, prior.upper.bound = prior.upper.bound, ... );
    
        .scalar <- .log.posterior.fn( posterior.maximum );
        if( be.verbose ) {
            cat( paste( "Done computing log posterior maximum for", sample.name, "- got", .scalar, "= .log.posterior.fn( posterior.maximum =", posterior.maximum, ")" ), fill = TRUE );
        }
        .scaled.unlogged.posterior.fn <- function ( x ) {
            ifelse( ( ( x <= prior.lower.bound ) | ( x >= prior.upper.bound ) ), 0, exp( .log.posterior.fn( x ) - .scalar ) )
        }
    } else {
        .scaled.unlogged.posterior.fn <- .posterior.fn;
    } # End if add.numeric.stability.at.the.sacrifice.of.speed .. else ..

    if( be.verbose ) {
        cat( paste( "Integrating posterior density for", sample.name, "over the range", .lower, "to", .upper, "with rel.tol =", total.mass.rel.tol ), fill = TRUE );
    }

    # Note we find you don't need the tolerance to be as precise here, and making total.mass.rel.tol on the order of 100 times less precise than rel.tol seems ok.
    total.mass <- as.numeric( integrate( .scaled.unlogged.posterior.fn, lower = .lower, upper = .upper, abs.tol = abs.tol, rel.tol = total.mass.rel.tol )$value );
    stopifnot( total.mass > 0 );

    if( be.verbose ) {
        cat( paste( "Done integrating posterior density for", sample.name, " - got total mass", total.mass ), fill = TRUE );
    }

    ## TODO: DEHACKIFY. FOR NOW be.verbose ONLY IF sample.name STARTS WITH "RNA"
    if( be.verbose && ( substr( sample.name, 1, 3 ) == "RNA" ) ) {
        function( x ) {
            .rv <- .scaled.unlogged.posterior.fn( x ) / total.mass;
            cat( paste( "Calculating normalized posterior density for", sample.name, "got", paste( x, .rv, sep = "=", collapse = "," ) ), fill = TRUE );
            return( .rv )
        }
    } else {
        function( x ) { .scaled.unlogged.posterior.fn( x ) / total.mass }
    }
} # create.normalized.posterior.density (..)



###################################################
### code chunk number 8: DSHIVInfectionTiming.Rnw:857-872
###################################################

## this function returns the posterior density for the prior event (by default if biomarker.dynamics.class == "RNA", it's the infection-causing-exposure time ), given boundaries on the time when the biomarker crosses the threshold, and a prior distribution on that time. 
create.biomarker.normalized.posterior.density <- function ( low.boundary = 0, high.boundary = 1, posterior.maximum = NULL, biomarker.dynamics.class = "RNA", the.biomarker.dynamics = biomarker_dynamics[[ biomarker.dynamics.class ]], prior.lower.bound = low.boundary - get.dynamics.maximum( the.biomarker.dynamics ), prior.upper.bound = high.boundary - get.dynamics.minimum( the.biomarker.dynamics ), be.verbose = FALSE, sample.name = "unnamed sample", abs.tol = 0.0001, rel.tol = 0.0001, subdivisions = 100L, ... ) {
    create.normalized.posterior.density( posterior.fn.creator.fn = create.unnormalized.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold.from.uncertainty.interval, low.boundary = low.boundary, high.boundary = high.boundary, posterior.maximum = posterior.maximum, the.biomarker.dynamics = the.biomarker.dynamics, prior.lower.bound = prior.lower.bound, prior.upper.bound = prior.upper.bound, be.verbose = be.verbose, sample.name = sample.name, rel.tol = rel.tol, abs.tol = abs.tol, subdivisions = subdivisions, ... )
} # create.biomarker.normalized.posterior.density (..)

# create.biomarker.assay.normalized.posterior.density <- function ( last.negative.time, first.positive.time, the.assay.dynamics, posterior.maximum = NULL, biomarker.dynamics.class = the.assay.dynamics$class, the.biomarker.dynamics = biomarker_dynamics[[ biomarker.dynamics.class ]], prior.lower.bound = last.negative.time - get.dynamics.maximum( the.biomarker.dynamics ), prior.upper.bound = last.negative.time - get.dynamics.minimum( the.biomarker.dynamics ), prior.over.biomarker.crossing.threshold.time = function ( time ) { 1 }, be.verbose = FALSE, sample.name = "unnamed sample", abs.tol = 0.0001, rel.tol = 0.0001, subdivisions = 100L, ... ) {
#     create.normalized.posterior.density( posterior.fn.creator.fn = create.unnormalized.posterior.on.time.of.preceding.event.from.biomarker.diagnostic.interval, last.negative.time = last.negative.time, first.positive.time = first.positive.time, the.assay.dynamics = the.assay.dynamics, the.biomarker.dynamics = the.biomarker.dynamics, posterior.maximum = posterior.maximum, prior.lower.bound = prior.lower.bound, prior.upper.bound = prior.upper.bound, prior.over.biomarker.crossing.threshold.time = prior.over.biomarker.crossing.threshold.time, be.verbose = be.verbose, sample.name = sample.name, rel.tol = rel.tol, abs.tol = abs.tol, subdivisions = subdivisions, ... )
# } # create.biomarker.assay.normalized.posterior.density (..)
# 
# ## For use when the assays differ for the two diagnostic interval endpoints (they should be for the same class!)
# create.biomarker.assay.normalized.posterior.density.twotests <- function ( last.negative.time, last.negative.assay.dynamics, first.positive.time, first.positive.assay.dynamics, posterior.maximum = NULL, biomarker.dynamics.class = last.negative.assay.dynamics$class, the.biomarker.dynamics = biomarker_dynamics[[ biomarker.dynamics.class ]], prior.lower.bound = last.negative.time - get.dynamics.maximum( the.biomarker.dynamics ), prior.upper.bound = last.negative.time - get.dynamics.minimum( the.biomarker.dynamics ), prior.over.biomarker.crossing.threshold.time = function ( time ) { 1 }, be.verbose = FALSE, sample.name = "unnamed sample", abs.tol = 0.0001, rel.tol = 0.0001, subdivisions = 100L, ... ) {
#     create.normalized.posterior.density( posterior.fn.creator.fn = create.unnormalized.posterior.on.time.of.preceding.event.from.biomarker.diagnostic.interval.twotests, last.negative.time = last.negative.time, last.negative.assay.dynamics = last.negative.assay.dynamics, first.positive.time = first.positive.time, first.positive.assay.dynamics = first.positive.assay.dynamics, the.biomarker.dynamics = the.biomarker.dynamics, posterior.maximum = posterior.maximum, prior.lower.bound = prior.lower.bound, prior.upper.bound = prior.upper.bound, prior.over.biomarker.crossing.threshold.time = prior.over.biomarker.crossing.threshold.time, be.verbose = be.verbose, sample.name = sample.name, rel.tol = rel.tol, abs.tol = abs.tol, subdivisions = subdivisions, ... )
# } # create.biomarker.assay.normalized.posterior.density.twotests (..)



###################################################
### code chunk number 9: DSHIVInfectionTiming.Rnw:875-1157
###################################################

## This is Phillip's, which I have copied here to help me with testing.
# ihist <- data.frame(
#   ptid = 'p0',
#   sample_date = c(17000.5, 17030.5, 17030.5, 17060.5),
#   test = c('abbott_real_time_weib3_delaney_and_manufacturer', 'architect_weib3_delaney', 'architect_weib3_delaney', 'geenius_fr_weib3_delaney'),
#   result = c('-', '+', '-', '+'),
#  stringsAsFactors = FALSE )

## This is Phillip's, but without the confusing instantaneous transition for the combo assay.
# ihist <- data.frame(
#   ptid = 'p0',
#   sample_date = c(17000.5, 17031.5, 17030.5, 17060.5),
#   test = c('abbott_real_time_weib3_delaney_and_manufacturer', 'architect_weib3_delaney', 'architect_weib3_delaney', 'geenius_fr_weib3_delaney'),
#   result = c('-', '+', '-', '+'),
#  stringsAsFactors = FALSE )

# These have just RNA results.
ihist.RNAonly.oneday1 <- data.frame(
  ptid = 'RNA-only',
  sample_date = c(17000.5, 17001.5),
  test = c('abbott_real_time_weib3_delaney_and_manufacturer', 'abbott_real_time_weib3_delaney_and_manufacturer' ),
  result = c('-', '+'),
  stringsAsFactors = FALSE )

ihist.RNAonly.oneday31 <- data.frame(
  ptid = 'RNA-only',
  #sample_date = c(17000.5, 17031.5),
  sample_date = c(17030.5, 17031.5),
  test = c('abbott_real_time_weib3_delaney_and_manufacturer', 'abbott_real_time_weib3_delaney_and_manufacturer' ),
  result = c('-', '+'),
  stringsAsFactors = FALSE )

ihist.RNAonly.31days <- data.frame(
  ptid = 'RNA-only',
  sample_date = c(17000.5, 17031.5),
  test = c('abbott_real_time_weib3_delaney_and_manufacturer', 'abbott_real_time_weib3_delaney_and_manufacturer' ),
  result = c('-', '+'),
  stringsAsFactors = FALSE )

## These have just combo
ihist.Comboonly.oneday31 <- data.frame(
  ptid = 'Combo-only',
  sample_date = c(17030.5, 17031.5),
  test = c('architect_weib3_delaney', 'architect_weib3_delaney' ),
  result = c('-', '+'),
  stringsAsFactors = FALSE )

ihist.Comboonly.31days <- data.frame(
  ptid = 'Combo-only',
  sample_date = c(17000.5, 17031.5),
  test = c('architect_weib3_delaney', 'architect_weib3_delaney' ),
  result = c('-', '+'),
  stringsAsFactors = FALSE )

ihist.Comboonly.oneday31.withunhelpfulRNA <- data.frame(
  ptid = 'Combo-only',
  sample_date = c(17031.5, 17030.5, 17031.5),
  test = c('abbott_real_time_weib3_delaney_and_manufacturer', 'architect_weib3_delaney', 'architect_weib3_delaney' ),
  result = c('+', '-', '+'),
  stringsAsFactors = FALSE )

## These have combo and RNA
ihist.RNACombo <- data.frame(
  ptid = 'Combo-RNA',
  sample_date = c(17000.5, 17001.5, 17030.5, 17031.5),
  test = c('abbott_real_time_weib3_delaney_and_manufacturer', 'abbott_real_time_weib3_delaney_and_manufacturer' ),
  result = c('-', '+'),
  stringsAsFactors = FALSE )


# The relevant bits of the all_assay_dynamics list are:
all_assay_dynamics <- list()
all_assay_dynamics[[tolower('abbott_real_time_weib3_delaney_and_manufacturer')]] <- list(
  class = 'RNA',
  full_assayname = 'Abbott Real Time HIV01 v1.0 m2000sp/m2000rt',
  short_assayname = 'Abbott Real Time',
  form = 'weib3',
  source = "manufacturer's details and delaney_2017",
  fun = 'weib3_assay_dynamics',
  params = list(location = 5.1252, shape = 1.350, scale = 9.660)
)

all_assay_dynamics[[tolower('architect_weib3_delaney')]] <- list(
  class = 'Ag/Ab',
  full_assayname = 'Abbott Architect HIV Ag/Ab Combo',
  short_assayname = 'Architect',
  form = 'weib3',
  source = 'delaney_2017',
  fun = 'weib3_assay_dynamics',
  params = list(location = 7.209, shape = 1.725, scale = 13.185)
)

all_assay_dynamics[[tolower('geenius_fr_weib3_delaney')]] <- list(
  class = 'IgG_Supp',
  full_assayname = 'BioRad Geenius Fully Reactive',
  short_assayname = 'Geenius_FR',
  form = 'weib3',
  source = 'delaney_2017',
  fun = 'weib3_assay_dynamics',
  params = list(location = 21.151, shape = 1.733, scale = 14.483)
)

## We go from latest to earliest. The adjusted.date values bound the time at which each biomarker crossed its threshold. So we start from the latest, integrating over that to get the prior for the next one, etc.
## The argument is a list data frames, one per biomarker in biomarker_dynamics and in that same order (which is the real-world order of those biomarkers crossing their respective thresholds: HIV-1 viral RNA, then p24 antigen, and then antibody), or fewer but only by removing them from the end of such a list (that is, all biomarkers should be present in order starting from RNA, but possibly not continuing all the way through to the last biomarker). It calls itself recursively, removing the elements of this list from the end and incorporating the evidence into the next prior.
create.normalized.posterior.density.of.infection.causing.exposure.time <- function ( remaining.biomarker.bounds.to.process.last.to.first, prior.over.last.biomarker.crossing.threshold.time = function ( time ) { 1 }, sample.name = "unnamed sample", ... ) {
    # Take the bottom one off, use it to construct a posterior, used as a prior; tail recurse. Base case: RNA.
    if( length( remaining.biomarker.bounds.to.process.last.to.first ) == 1 ) {
        the.biomarker <- names( remaining.biomarker.bounds.to.process.last.to.first )[ 1 ];
        ## Our process requires all biomarkers to be present, and RNA is the first one, so when only one is left it should be RNA.
        stopifnot( the.biomarker == "RNA" );
    } else {
        the.biomarker <- names( remaining.biomarker.bounds.to.process.last.to.first )[ length( remaining.biomarker.bounds.to.process.last.to.first ) ];
    }

    ## ERE I AM, got to the narrowest SNOTE debugging place, have to sleep now tho.
    # if( is.null( prior.over.last.biomarker.crossing.threshold.time ) ) {
    #     prior.over.last.biomarker.crossing.threshold.time.is.uniform <- TRUE;
    # } else {
    #     prior.over.last.biomarker.crossing.threshold.time.is.uniform <- 
    #         is.uniform.over.range( prior.over.last.biomarker.crossing.threshold.time, lower, upper );
    #     if( prior.over.last.biomarker.crossing.threshold.time.is.uniform ) {
    #         if( be.verbose ) {
    #             cat( "Replacing the given prior with a uniform prior for efficiency as it is effectively uniform" );
    #         }
    #     }
    # }
    # if( prior.over.last.biomarker.crossing.threshold.time.is.uniform ) {
    #    prior.over.last.biomarker.crossing.threshold.time <- NULL;
    # }
    
    ## boundaries on the biomarker threshold-crossing time are given in the "adjusted.date" column.
    low.boundary <- remaining.biomarker.bounds.to.process.last.to.first[[ the.biomarker ]][ 1, "adjusted.date" ];
    high.boundary <- remaining.biomarker.bounds.to.process.last.to.first[[ the.biomarker ]][ 2, "adjusted.date" ];
    
    if( length( remaining.biomarker.bounds.to.process.last.to.first ) == 1 ) {
        return( create.biomarker.normalized.posterior.density( low.boundary = low.boundary, high.boundary = high.boundary, posterior.maximum = NULL, biomarker.dynamics.class = "RNA", prior.over.biomarker.crossing.threshold.time = prior.over.last.biomarker.crossing.threshold.time, sample.name = paste( the.biomarker, remaining.biomarker.bounds.to.process.last.to.first[[ the.biomarker ]][ 1, "ptid" ] ), ... ) );
    }

    ## We process these from latest to earliest. The adjusted.date values bound the time at which each biomarker crossed its threshold. So we start from the latest, integrating over that to get the prior for the next one, etc. That's what this function does:

    posterior.for.preceding.event <-
        create.biomarker.normalized.posterior.density( low.boundary = low.boundary, high.boundary = high.boundary, posterior.maximum = NULL, biomarker.dynamics.class = the.biomarker, prior.over.biomarker.crossing.threshold.time = prior.over.last.biomarker.crossing.threshold.time, sample.name = paste( the.biomarker, remaining.biomarker.bounds.to.process.last.to.first[[ the.biomarker ]][ 1, "ptid" ] ), ... );
    remaining.biomarker.bounds.to.process.last.to.first.excluding.this.one <-
        remaining.biomarker.bounds.to.process.last.to.first[ -length( remaining.biomarker.bounds.to.process.last.to.first ) ];
    return( create.normalized.posterior.density.of.infection.causing.exposure.time( remaining.biomarker.bounds.to.process.last.to.first = remaining.biomarker.bounds.to.process.last.to.first.excluding.this.one, prior.over.last.biomarker.crossing.threshold.time = posterior.for.preceding.event, ... ) );
} # create.normalized.posterior.density.of.infection.causing.exposure.time ( .. )

# Prepare the ihist such that it has all biomarkers (from RNA through the last biomarker with any observed diagnostic data in the given ihist) represented with exactly one negative and exactly one positive result (the most constraining, and defaults or borrowing from other biomarker categories where data are missing).
prepare.diagnostic.data.for.analysis <- function ( ihist, earliest.sample.date = 0, latest.sample.date = as.numeric( Sys.Date() ), be.verbose = FALSE ) {
    stopifnot( class( ihist ) == 'data.frame' );
    stopifnot( all( names( ihist ) == c( 'ptid', 'sample_date', 'test', 'result' ) ) );
    stopifnot( nrow( ihist ) >= 1 );
    
    ihist.with.biomarker <- ihist;

    # Add the biomarker as a fifth column (character vector)
    ihist.biomarker.classes <- sapply( all_assay_dynamics[ as.character( ihist$test ) ], function( the.dynamics ) { as.character( the.dynamics$class ) } );
    ihist.with.biomarker$biomarker <- ihist.biomarker.classes;
    
    ## Add the within-class delay of the assay as a sixth column.
    ihist.assay.delay <- apply( ihist.with.biomarker, 1, function( .row ) { get.time.from.biomarker.crossing.threshold.to.assay.positivity( all_assay_dynamics[[ .row[ "test" ] ]], the.biomarker.dynamics = biomarker_dynamics[[ .row[ "biomarker" ] ]] ) } );
    ihist.with.biomarker$delay <- ihist.assay.delay;
    ## Add the adjusted date as a seventh column. This is the sample date minus the assay delay, so it reflects the constraints on the tie at which the biomarker crossed the canonical "1 copy per ml" level that are implied by the assay results.
    ihist.with.biomarker$adjusted.date <- ihist$sample_date - ihist.assay.delay;

    ihist.with.biomarker$result <- as.character( ihist.with.biomarker$result );
    ihist.with.biomarker$test <- as.character( ihist.with.biomarker$test );
    
    ## Order the results so that faster biomarkers are at the top, and within a biomarker class we have last negative before first positive, and within those the earliest (adjusted) dates appear first. We can just ignore all biomarkers past the last one we have available.
   
    the.biomarkers.in.increasing.order <- names( biomarker_dynamics );
    if( length( setdiff( the.biomarkers.in.increasing.order, ihist.with.biomarker$biomarker ) ) > 0 ) {
        .last.biomarker.i <- length( the.biomarkers.in.increasing.order );
        while( ( .last.biomarker.i > 0 ) && !( the.biomarkers.in.increasing.order[ .last.biomarker.i ] %in% ihist.with.biomarker$biomarker ) ) { .last.biomarker.i <- .last.biomarker.i - 1; }
        if( .last.biomarker.i != length( the.biomarkers.in.increasing.order ) ) {
            if( be.verbose ) { cat( paste( "NOTE: No diagnostics available for biomarkers at/after", the.biomarkers.in.increasing.order[ .last.biomarker.i + 1 ] ), fill = TRUE ) }
            the.biomarkers.in.increasing.order <-
                the.biomarkers.in.increasing.order[ 1:.last.biomarker.i ];
        }
    } # End if there might be an opportunity to shrink trailing biomarkers out of consideration.
  
    ## Since we assume that within each class there's only a linear shift, then there's a single best last negative and best first positive.
    the.ptid <- ihist$ptid[ 1 ];
    stopifnot( all( ihist$ptid == the.ptid ) );
    .last.negative.result <- ihist.with.biomarker[ 1, , drop = FALSE ];
    .last.negative.result[ "sample_date" ] <- earliest.sample.date;
    .last.negative.result[ "test" ] <- "default";
    .last.negative.result[ "result" ] <- "-";
    .last.negative.result[ "biomarker" ] <- NA;
    .last.negative.result[ "adjusted.date" ] <- earliest.sample.date;
    .first.positive.result <- .last.negative.result;
    .first.positive.result[ "result" ] <- "+";
    .first.positive.result[ "sample_date" ] <- latest.sample.date;
    .first.positive.result[ "adjusted.date" ] <- latest.sample.date;
    ihist.reordered.by.biomarker <- lapply( 1:length( the.biomarkers.in.increasing.order ), function( .biomarker.i ) {
        the.biomarker <- the.biomarkers.in.increasing.order[ .biomarker.i ];
         .ihist.negative.results.for.biomarker <-
            ihist.with.biomarker[ ( ihist.with.biomarker$biomarker == the.biomarker ) & ( ihist.with.biomarker$result == "-" ), , drop = FALSE ];
        if( nrow( .ihist.negative.results.for.biomarker ) == 0 ) {
            .ihist.negative.results.for.biomarker <- .last.negative.result;
            .ihist.negative.results.for.biomarker[ "biomarker" ] <- the.biomarker;
        } else if( nrow( .ihist.negative.results.for.biomarker ) > 1 ) {
            .ihist.negative.results.for.biomarker <-
                .ihist.negative.results.for.biomarker[ order( .ihist.negative.results.for.biomarker$adjusted_date ), , drop = FALSE ];
        }
        # As we go, update the last negative. Note we'll have to sweep back through to fix the first positives, but we might as well do this as we go here:
        .last.negative.result <<-
            .ihist.negative.results.for.biomarker[ nrow( .ihist.negative.results.for.biomarker ), , drop = FALSE ];
  
        .ihist.positive.results.for.biomarker <-
            ihist.with.biomarker[ ( ihist.with.biomarker$biomarker == the.biomarker ) & ( ihist.with.biomarker$result == "+" ), , drop = FALSE ];
        if( nrow( .ihist.positive.results.for.biomarker ) == 0 ) {
            .ihist.positive.results.for.biomarker <- .first.positive.result;
            .ihist.positive.results.for.biomarker[ "biomarker" ] <- the.biomarker;
        } else if( nrow( .ihist.positive.results.for.biomarker ) > 1 ) {
            .ihist.positive.results.for.biomarker <-
                .ihist.positive.results.for.biomarker[ order( .ihist.positive.results.for.biomarker$adjusted_date ), , drop = FALSE ];
        }
        # As we go, we update the last negative. But we'll have to sweep back through to fix the first positives, since the dependence is the other way.
        # .first.positive.result <<-
        #     .ihist.positive.results.for.biomarker[ 1, , drop = FALSE ];
  
        return( rbind( .ihist.negative.results.for.biomarker, .ihist.positive.results.for.biomarker ) );
    } );
    names( ihist.reordered.by.biomarker ) <- the.biomarkers.in.increasing.order;
    
    ihist.reordered.by.biomarker.fixed <- lapply( 1:length( the.biomarkers.in.increasing.order ), function( .biomarker.i ) {
        .rv <- ihist.reordered.by.biomarker[[ .biomarker.i ]];
        if( .biomarker.i == length( the.biomarkers.in.increasing.order ) ) {
            # nothing to do.
            return( .rv );
        }
        if( "default" %in% ihist.reordered.by.biomarker[[ .biomarker.i ]][ ihist.reordered.by.biomarker[[ .biomarker.i ]][ , "result" ] == "+", "test" ] ) {
            .first.positive.i <-
                which( ( ihist.reordered.by.biomarker[[ .biomarker.i ]][ , "test" ] == "default" ) &
                      ( ihist.reordered.by.biomarker[[ .biomarker.i ]][ , "result" ] == "+" ) );
            stopifnot( length( .first.positive.i ) == 1 );
            .first.positive.result <-
                ihist.reordered.by.biomarker[[ .biomarker.i + 1 ]][ nrow( ihist.reordered.by.biomarker[[ .biomarker.i + 1 ]] ), , drop = FALSE ];
            .rv[ .first.positive.i, ] <- .first.positive.result;
        }
        return( .rv );
    } );
    names( ihist.reordered.by.biomarker.fixed ) <- the.biomarkers.in.increasing.order;
    
    ## The most informative are determined by the adjusted.date (latest negative, earliest positive).
    ihist.most.informative.by.biomarker <- lapply( 1:length( the.biomarkers.in.increasing.order ), function( .biomarker.i ) {
        .first.positive.i <- base::min( which( ihist.reordered.by.biomarker.fixed[[ .biomarker.i ]][ , "result" ] == "+" ) );
        .last.negative.first.positive.for.biomarker <- ihist.reordered.by.biomarker.fixed[[ .biomarker.i ]][ .first.positive.i + -1:0, , drop = FALSE ];
        return( .last.negative.first.positive.for.biomarker );
    } );
    names( ihist.most.informative.by.biomarker ) <- the.biomarkers.in.increasing.order;
  
    return( ihist.most.informative.by.biomarker );
} # prepare.diagnostic.data.for.analysis ( ihist )

construct.aggregate.interpreter <- function ( ihist, prior.lower.bound = -Inf, prior.upper.bound = Inf, lower.bound.offset.from.first.positive = -100, minimum.bound.width = 150, be.verbose = FALSE, ... ) {
    # The prior bounds can be brought in a bit if the data have anything to say.
    if( !is.finite( prior.lower.bound ) ) {
        prior.lower.bound <- ( base::min( ihist[ ihist[ , "result" ] == "+", "sample_date" ] ) + lower.bound.offset.from.first.positive );
    }
    if( !is.finite( prior.upper.bound ) ) {
        # Have it be the later of: a) the last positive test, or b) the last negative test + minimum.bound.width.
        prior.upper.bound <-
            ( base::max( ihist[ ihist[ , "result" ] == "+", "sample_date" ] ) );
        # Also make sure it's at least 50 days after min, for examples like that.
        prior.upper.bound <- base::max( prior.upper.bound, prior.lower.bound + minimum.bound.width );
    }

    if( be.verbose ) {
        cat( "preparing diagnostic data for aggregation: sorting and selecting most informative by biomarker..", fill = TRUE );
    }
    ihist.most.informative.by.biomarker <- prepare.diagnostic.data.for.analysis( ihist, be.verbose = be.verbose, earliest.sample.date = prior.lower.bound, latest.sample.date = prior.upper.bound );
    if( be.verbose ) {
        cat( "..done.", fill = TRUE );
    }
  
    ## We process these from latest to earliest. The adjusted.date values bound the time at which each biomarker crossed its threshold. So we start from the latest, integrating over that to get the prior for the next one, etc. That's what this function does:
    return( create.normalized.posterior.density.of.infection.causing.exposure.time( ihist.most.informative.by.biomarker, prior.lower.bound = prior.lower.bound, prior.upper.bound = prior.upper.bound, be.verbose = be.verbose, ... ) );
} # construct.aggregate.interpreter ( ihist )



###################################################
### code chunk number 10: DSHIVInfectionTiming.Rnw:1160-1349
###################################################

## This is a general-purpose helper function.

# Ok, so now we need to determine specific quantiles.  I could do this using the memoizing integration table I developed for the MBS code, that hacks it right into integrate(), but instead I'll do the following.
## Given a list.of.desired.quantiles, I can start by dividing the entire range from prior.lower.bound to prior.upper.bound into K (eg, 10) partitions.  I keep a list of length K+1 with a tuple at each position: (x, density.at.x, cdf.at.x, desired.quantiles.in.bin ), where desired.quantiles.in.bin are the results, named by the quantile, corresponding to those of the all.desired.quantiles that are below the "cdf.at.x" at the _previous_ position, but that are _not_ below the "cdf.at.x" at this position.  Thus eg c( "0.5" = 10.234, "0.9" = 103.201 ).  "x" marks the top of the bin, and the first x to be evaluated is 0, then 1/K, .. up to 1.  For x=prior.lower.bound the tuple has no desired.quantiles.in.bin unless the function has non-neglible density at this point; then it is all of the desired quantiles that fall below that density, which is also the cdf value there.  All but the first bin can be further broken down by repeating the process for any bin for which desired.quantiles.in.bin is non-empty. We declare to have found an individual quantile whenever abs( cdf.at.x - desired.quantile ) < tolerance, and we continue the recursion until all quantiles are found.

## To speed this up, if the lower limit is below 0.5, we just call it zero and be done quickly.  A lot of time used to be spent carefully resolving the lower limit and this is wasted since we effectively round at the end anyway.
compute.desired.quantiles.results.recursively <- function ( ihist, density.fn = NULL, density.fn.creator.fn = construct.aggregate.interpreter, num.partitions = 10, all.desired.quantiles = c( 0.025, 0.5, 0.975 ), prior.lower.bound = -Inf, prior.upper.bound = Inf, lower.bound.offset.from.first.positive = -100, minimum.bound.width.when.no.prior.upper.bound.is.specified = 150, be.verbose = FALSE, sample.name = as.character( ihist[ 1, "ptid" ] ), abs.tol = 1E-4, rel.tol = 1E-2, subdivisions = 100L, ... ) {
    #print( num.partitions );
    if( be.verbose ) {
        cat( "computing quantiles results for", sample.name, fill = TRUE );
    }
    # The prior bounds can be brought in a bit if the data have anything to say. For now as a hack use a year before the earliest mentioned date for the lower bounds, and a year after the latest mentioned sample date for upper bounds. 
    if( !is.finite( prior.lower.bound ) ) {
        prior.lower.bound <- ( base::min( ihist[ ihist[ , "result" ] == "+", "sample_date" ] ) + lower.bound.offset.from.first.positive );
    }
    if( !is.finite( prior.upper.bound ) ) {
        # Have it be the later of: a) the last positive test, or b) the last negative test + minimum.bound.width.when.no.prior.upper.bound.is.specified.
        prior.upper.bound <-
            ( base::max( ihist[ ihist[ , "result" ] == "+", "sample_date" ] ) );
        # Also make sure it's at least minimum.bound.width.when.no.prior.upper.bound.is.specified days after min, for examples like that.
        prior.upper.bound <-
            base::max( prior.upper.bound, prior.lower.bound + minimum.bound.width.when.no.prior.upper.bound.is.specified );
    }

    recursively.compute.desired.quantiles <- function ( density.fn, lower = prior.lower.bound, upper = prior.upper.bound, cdf.at.lower.bound = with( list2env( list( ...fn = density.fn, ...lower = lower ) ), as.numeric( integrate( ...fn, lower=0, upper=...lower, rel.tol = rel.tol, abs.tol = abs.tol, subdivisions = subdivisions )$value ) ), desired.quantiles = all.desired.quantiles, include.lower.point.partition = ( lower == prior.lower.bound ) ) {
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
            cat( paste( "after calling recursively.compute.desired.quantiles.helper, got: ( prev.x=", prev.x, ", cdf.at.prev.x=", cdf.at.prev.x, ", quantile.results.in.bin=c(", paste( quantile.results.in.bin, collapse = "," ), "), names( quantile.results.in.bin )=c(", paste( names( quantile.results.in.bin ), collapse = "," ), ")", sep = "" ), fill = TRUE )
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
            cat( paste( "recursively.compute.desired.quantiles.helper( partition.number=", partition.number, ", lower=", lower, ", upper=", upper, ", prev.x=", prev.x, ", cdf.at.prev.x=", cdf.at.prev.x, ", desired.quantiles=c(", paste( desired.quantiles, collapse = "," ), ") )", sep = "" ), fill = TRUE )
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
            cdf.within.bin <- as.numeric( integrate( density.fn, prev.x, x, rel.tol = rel.tol, abs.tol = abs.tol, subdivisions = subdivisions )$value );
            # This can happen when the integral is a poor approximation at the upper bound so that it's actually larger using a smaller upper bound.  So don't "stopifnot"
            # if( cdf.within.bin > 1.0+0.01 ) {
            #     ## TODO: REMOVE
            #     print( paste("CDF WITHIN BIN:",cdf.within.bin));
            # }
            # stopifnot( cdf.within.bin <= 1.0+0.01 );
            cdf.at.x <- cdf.at.prev.x + cdf.within.bin;
        }
        desired.quantiles.in.bin <-
            desired.quantiles[ ( desired.quantiles <= ( cdf.at.x + rel.tol ) ) ];
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
                desired.quantiles.in.bin[ abs( desired.quantiles.in.bin - cdf.at.x ) < rel.tol ];
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
    if( is.null( density.fn ) ) {
        density.fn <-
            density.fn.creator.fn( ihist = ihist, be.verbose = be.verbose, prior.lower.bound = prior.lower.bound, prior.upper.bound = prior.upper.bound, sample.name = sample.name, rel.tol = rel.tol, abs.tol = abs.tol, subdivisions = subdivisions, ... );
    }
    desired.quantiles <-
        recursively.compute.desired.quantiles( density.fn );
    return( desired.quantiles );
} # compute.desired.quantiles.results.recursively (..)

# NOTE that unlike estimate_lb_med_ub(..), this does not return "max_agg", which appears to be not used anywhere -- excluding because it would be slow to compute and it is not used. Use compute.density.maximum(..) if you want it.
## Here we ASSUME THAT fun is NORMALIZED ALREADY. THAT IS, IT SUMS TO 1.
estimate.lb.med.ub.for.density.fun <- function ( fun, range_start, range_end, verbose = FALSE, label = 'unlabeled', 
                                 warn_low_AOC = NA, extra_tiles = NULL, date_splits = NULL ) {
    desired.quantiles <- compute.desired.quantiles.results.recursively( density.fn = fun, all.desired.quantiles = c( 0.025, 0.5, 0.975, extra_tiles ), prior.lower.bound = range_start, prior.upper.bound = range_end, be.verbose = verbose, sample.name = ifelse( label == 'unlabeled', as.character( ihist[ 1, "ptid" ] ), label ) );
    extra.computed.quantiles <- NULL;
    if( !is.null( extra_tiles ) ) {
        extra.computed.quantiles <- desired.quantiles[ 4:length( desired.quantiles ) ];
    }

    aoc_left_of_date <- NULL;
    if( !is.null( date_splits ) ) {
      for (indx in 1:length(date_splits)){
        stopifnot(date_splits[indx] > range_start)
        stopifnot(date_splits[indx] < range_end)
        if (date_splits[indx]%%1 != 0.5){
          warning('Some date split not at midday')
        }
      }
      last.date <- range_start;
      cum.sum <- 0;
      for (indx in 1:length(date_splits)){
          # Note that I (Paul) am using integrate rather than pracma::integral as Phillip had used.
          mass.of.interval <- integrate( fun, lower = last.date, upper = date_splits[ indx ] )[[1]];
          cum.sum <- cum.sum + mass.of.interval;
          aoc_left_of_date <- c( aoc_left_of_date, cum.sum );
          last.date <- date_splits[ indx ];
      }
    }

    return( list( lb = unname( desired.quantiles[ 1 ] ), med = unname( desired.quantiles[ 2 ] ), ub = unname( desired.quantiles[ 3 ] ), aoc = 1, extra_tiles = extra_tiles, extra_computed_tiles = extra.computed.quantiles, date_splits = date_splits, aoc_left_of_date = aoc_left_of_date ) );
} # estimate.lb.med.ub.for.density.fun (..)



###################################################
### code chunk number 11: DSHIVInfectionTiming.Rnw:1353-1357
###################################################
# With an interval width of 10 between the biomarker results, the prior on the location in that interval (eg from a Weibull posterior induced by a combo result) makes a difference:
# compute.desired.quantiles.results.recursively( prior.over.biomarker.crossing.threshold.time = function( time ) { dweibull( time + 5 - 4.8, shape = 10, scale = 10 ) }, high.boundary = 10 )
# 0.025   0.5 0.975 
#    14    20    37 


###################################################
### code chunk number 12: DSHIVInfectionTiming.Rnw:1361-1364
###################################################
# compute.desired.quantiles.results.recursively( prior.over.biomarker.crossing.threshold.time = function( time ) { dweibull( time + 15 - 4.8, shape = 10, scale = 10 ) }, high.boundary = 10 )
# 0.025   0.5 0.975 
#     6    12    29 


###################################################
### code chunk number 13: DSHIVInfectionTiming.Rnw:1368-1426
###################################################
##ENDCOPYMARK

# if( FALSE ) {
#     # the areas, by time.
#     area.under.unnormalized.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold.from.uncertainty.interval <- function( low.boundary, high.boundary, biomarker.dynamics.class = "RNA", the.biomarker.dynamics = biomarker_dynamics[[ biomarker.dynamics.class ]], prior.lower.bound = get.dynamics.minimum( the.biomarker.dynamics ), prior.upper.bound = get.dynamics.maximum( the.biomarker.dynamics ), prior.over.biomarker.crossing.threshold.time = function ( shift ) { 1 }, prior.over.time = function ( time ) { 1 }, subdivisions = 100L, log.results = FALSE, ... ) {
#         .fn <- create.unnormalized.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold.from.uncertainty.interval( low.boundary, high.boundary, prior.over.biomarker.crossing.threshold.time = prior.over.biomarker.crossing.threshold.time, prior.lower.bound = prior.lower.bound, prior.upper.bound = prior.upper.bound, ... );
#         .x <- integrate( function( time ) { prior.over.time( time ) * .fn( time ) }, lower = prior.lower.bound, upper = prior.upper.bound, subdivisions = subdivisions )[[1]]
#         ifelse( log.results, log( .x ), .x );
#     } # area.under.unnormalized.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold.from.uncertainty.interval ( .. )
#     
#     area.under.unnormalized.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold.using.uniform.prior.by.day <- sapply( 0:100, function( day ) { area.under.unnormalized.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold.from.uncertainty.interval( low.boundary = day - 0.5, high.boundary = day + 0.5, prior.upper.bound = 150 ); } );
#     
#     # The weibull is shifted by 10 days...
#     area.under.unnormalized.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold.using.weibull.prior.by.day <- sapply( 0:100, function( day ) { area.under.unnormalized.density.of.time.of.preceding.event.from.time.of.biomarker.crossing.threshold.from.uncertainty.interval( low.boundary = day - 0.5, high.boundary = day + 0.5, prior.upper.bound = 150, prior.over.biomarker.crossing.threshold.time = function( time ) { dweibull( time + 10 - 4.8, shape = 1.35, scale = 9 ) } ); } );
# } # END IF ( FALSE )

# marginalized.posterior.from.rna.pos.result.on.time.of.preceding.event <- Vectorize( function ( time.center, interval.width, prior.over.time = function ( time ) { 1 } ) {
#     unnormalized.marginalized.posterior.from.rna.pos.result.on.time.of.preceding.event( time.center, interval.width, prior.over.time = prior.over.time ) /
#     area.under.unnormalized.marginalized.posterior.from.rna.pos.result.for.interval.width( interval.width, prior.over.time = prior.over.time )
# }, vectorize.args = "time.center" ); # marginalized.posterior.from.rna.pos.result.on.time.of.preceding.event (..)
#     
# integrate( marginalized.posterior.from.rna.pos.result.on.time.of.preceding.event, interval.width = 10, lower = 0, upper = 100 )
# # 1 with absolute error < 8.6e-05
# integrate( marginalized.posterior.from.rna.pos.result.on.time.of.preceding.event, interval.width = 20, lower = 0, upper = 100 )
# # 1 with absolute error < 0.00011
# 
# .pod <- function( time ) { dweibull( time - 4.8, shape = 1.35, scale = 9 ) };
# .fn <- function( time ) { 
# integrate( marginalized.posterior.from.rna.pos.result.on.time.of.preceding.event, interval.width = 20, prior.over.time = .pod, lower = 0, upper = 100 )
# 
# ## ENDMARK
# 
# pdf( "Marginalized_posterior_on_time_since_infection_given_RNA_negative_for_interval_width_1.pdf" )
# plot( (1:5000/100), marginalized.posterior.from.rna.pos.result.on.time.of.preceding.event( (1:5000/100), interval.width = 1 ), xlim = c( 0, 50 ), pch = 20 , main = "Posterior prob of time since infection interval of 1 day", xlab = "Day", ylab = "Density" )
# dev.off()
# 
# pdf( "Marginalized_posterior_on_time_since_infection_given_RNA_negative_for_interval_width_10.pdf" )
# plot( (1:5000/100), marginalized.posterior.from.rna.pos.result.on.time.of.preceding.event( (1:5000/100), interval.width = 10 ), xlim = c( 0, 50 ), pch = 20 , main = "Posterior prob of time since infection interval of 10 days", xlab = "Day", ylab = "Density" )
# dev.off()
# 
# pdf( "Marginalized_posterior_on_time_since_infection_given_RNA_negative_for_interval_width_20.pdf" )
# plot( (1:5000/100), marginalized.posterior.from.rna.pos.result.on.time.of.preceding.event( (1:5000/100), interval.width = 20 ), xlim = c( 0, 50 ), pch = 20 , main = "Posterior prob of time since infection interval of 20 days", xlab = "Day", ylab = "Density" )
# dev.off()
# 
# ## Now since the later tests are all only evaluated based on their time since RNA positivity, ie Combo and Ab assays, and if we assume a similar form (3 parameter Weibull) then we have essentially to use its posterior as a prior over the evidence from the RNA tests. 
# 
# pdf( "Marginalized_posterior_on_time_since_infection_given_RNA_negative_with_weibull_prior_for_interval_width_1.pdf" )
# plot( (1:5000/100), marginalized.posterior.from.rna.pos.result.on.time.of.preceding.event( (1:5000/100), interval.width = 1, prior.over.days = function( days ) { dweibull( days - 4.8, shape = 1.35, scale = 9 ) } ), xlim = c( 0, 50 ), pch = 20 , main = "Posterior prob of time since infection interval of 1 day", xlab = "Day", ylab = "Density" )
# dev.off()
# 
# pdf( "Marginalized_posterior_on_time_since_infection_given_RNA_negative_with_weibull_prior_for_interval_width_10.pdf" )
# plot( (1:5000/100), marginalized.posterior.from.rna.pos.result.on.time.of.preceding.event( (1:5000/100), interval.width = 10, prior.over.days = function( days ) { dweibull( days - 4.8, shape = 1.35, scale = 9 ) } ), xlim = c( 0, 50 ), pch = 20 , main = "Posterior prob of time since infection interval of 10 days", xlab = "Day", ylab = "Density" )
# dev.off()
# 
# pdf( "Marginalized_posterior_on_time_since_infection_given_RNA_negative_with_weibull_prior_for_interval_width_20.pdf" )
# plot( (1:5000/100), marginalized.posterior.from.rna.pos.result.on.time.of.preceding.event( (1:5000/100), interval.width = 20, prior.over.days = function( days ) { dweibull( days - 4.8, shape = 1.35, scale = 9 ) } ), xlim = c( 0, 50 ), pch = 20 , main = "Posterior prob of time since infection interval of 20 days", xlab = "Day", ylab = "Density" )
# dev.off()



