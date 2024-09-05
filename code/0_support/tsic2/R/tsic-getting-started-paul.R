## To ensure you don't accidentally keep stuff in memory and get buggy output if there's an unnoticed Error, do this first:
rm( list = ls() )

# This loads and sources all necessary background stuff:
source( "tsic-getting-started-fromsource.R" )

##############
## Params
## Choose one scenario among names( scenarios ) [see below]
scenario <- "ComboRNA";
# Make the RNA negative date X.1.days.apart days before the positive date:
X.1.days.apart <- 7;
# Make the other one negative date X.2 days before the positive date:
X.2.days.apart <- 28;
# Shift them so the RNA dates are Y days earlier.
Y.days.earlier <- 0;

#### Usually don't need to change these:
# Align them so the RNA + date is on the Combo - date, keeping diagnostic uncertainty window lengths.
align.diagnosis.dates <- TRUE;
# Make the scales for the aggregate curves free, rather than forced to go to 1 - especially helpful for posterior density plotting.
scales <- "free_y"; # Note this is a change from what Phillip had shown.
##############

##############
## Data
## This is the default ihist shown in isic-getting-started.Rmd. See below where we also can use pre-defined ihists.
# ihist <- data.frame(
#   ptid = c('p0', 'p0'),
#   sample_date = c(as.numeric(as.Date('2019-04-01')), as.numeric(as.Date('2019-05-01'))),
#   test = c('architect_weib3_delaney', 'architect_weib3_delaney'),
#   result = c('-', '+'),
#   stringsAsFactors = FALSE
# )
# ihist$sample_date <- ihist$sample_date + 0.5

# This is exactly the same as ihist.Comboonly.oneday6
# ihist <- data.frame(
#   ptid = c('p0', 'p0'),
#   sample_date = c(as.numeric( as.Date( '2019-04-20' ) ), as.numeric( as.Date( '2019-04-21' ) )),
#   test = c('architect_weib3_delaney', 'architect_weib3_delaney'),
#   result = c('-', '+'),
#   stringsAsFactors = FALSE
# )
# ihist$sample_date <- ihist$sample_date + 0.5

#### Choose one scenario:
scenarios <-
    list(
         "RNA" = list(
                      biomarker.with.fastest.result = "RNA",
                      ihist = ihist.RNAonly.oneday1,
                      output.filename.base = "RNA" ),
         "Combo" = list(
                      biomarker.with.fastest.result = "Combo",
                      ihist = ihist.Comboonly.oneday6,
                      output.filename.base = "Combo" ),
         "Ab" = list(
                      biomarker.with.fastest.result = "Ab",
                      ihist = ihist.Abonly.oneday21,
                      output.filename.base = "Ab" ),
         "ComboRNA" = list(
                      biomarker.with.fastest.result = "RNA",
                      ihist = ihist.ComboRNA.oneday6.oneday1,
                      output.filename.base = "ComboRNA" )
         );

stopifnot( scenario %in% names( scenarios ) );

## Add this scenario to the environment.
cat( paste( "Using scenario", scenario ), fill = TRUE );
attach( scenarios[[ scenario ]], name = scenario );

## Make the RNA negative date X.1.days.apart days before the positive date:
output.filename.base <- paste( output.filename.base, X.1.days.apart, paste( "days", biomarker.with.fastest.result, sep = "" ), sep = "" );
ihist$sample_date[ 1 ] <- ( ihist$sample_date[ 2 ] - X.1.days.apart );

## Make the other one negative date X.2 days before the positive date:
if( nrow( ihist ) > 2 ) {
    ihist$sample_date[ 3 ] <- ( ihist$sample_date[ 4 ] - X.2.days.apart );
    output.filename.base <- paste( output.filename.base, X.2.days.apart, "daysCombo", sep = "" );
}

## Align them so the RNA + date is on the Combo - date, keeping diagnostic uncertainty window lengths.
if( align.diagnosis.dates && ( nrow( ihist ) > 2 ) ) {
     days.to.shift.by.to.align.times.across.tests <-
         ( ihist$sample_date[ 2 ] - ihist$sample_date[ 3 ] );
    ihist$sample_date[ 1:2 ] <-
        ihist$sample_date[ 1:2 ] - days.to.shift.by.to.align.times.across.tests;
     #output.filename.base <- paste( output.filename.base, "RNAPosAtComboNeg", sep = "" );
}
## Shift them so the RNA dates are Y days earlier.
if( ( Y.days.earlier != 0 ) && ( nrow( ihist ) > 2 ) ) {
    ihist$sample_date[ 1:2 ] <-
        ihist$sample_date[ 1:2 ] - Y.days.earlier;
    output.filename.base <- paste( output.filename.base, "RNAshifted", Y.days.earlier, "days", sep = "" );
}

##############

##############
## More setup after setting the ihist data input value:
range_start <- min(ihist$sample_date) - 100
range_end <- max(ihist$sample_date) + 100
##############

##############
## Functions
## USE BOTH OLD AND NEW!
## OLD:
aggregate_interpreter_constructor_fn.old <- construct_aggregate_interpreter;
estimate_lb_med_ub_fn.old <- estimate_lb_med_ub;
## NEW:
aggregate_interpreter_constructor_fn.new <- construct.aggregate.interpreter;
estimate_lb_med_ub_fn.new <- estimate.lb.med.ub.for.density.fun;

agg_interpreter.old <- aggregate_interpreter_constructor_fn.old(ihist)

## By doing it the second of these ways you can extract the list of inputs to the final call (the base case, for RNA) to get at the prior used there.n
#agg_interpreter.new <- aggregate_interpreter_constructor_fn.new(ihist)
agg_interpreter.new.in.list.with.inputs <- aggregate_interpreter_constructor_fn.new(ihist, return.list.with.inputs = TRUE )
agg_interpreter.new <- agg_interpreter.new.in.list.with.inputs[[ "normalized.posterior.density.of.infection.causing.exposure.time" ]];
##############

## and some usage examples
# agg_interpreter.old(17928) # 2019-02-01
# agg_interpreter.new(17928) # 2019-02-01
# agg_interpreter.old(17988) # 2019-04-02
# agg_interpreter.new(17988) # 2019-04-02
# agg_interpreter.old(18020) # 2019-05-04
# agg_interpreter.new(18020) # 2019-05-04

##############
## Results: quantiles
lb_med_ub.old <- estimate_lb_med_ub_fn.old(fun = agg_interpreter.old,
                                range_start = range_start,
                                range_end = range_end,
                                   verbose = FALSE)
lb_med_ub.new <- estimate_lb_med_ub_fn.new(fun = agg_interpreter.new,
                                range_start = range_start,
                                range_end = range_end,
                                           verbose = FALSE)
lb_med_ub.new[[ "max_agg" ]] <- NA;

## The result as two columns. Below we add a third for the prior.
results.table <- cbind( "CDFs" = unlist( lb_med_ub.old ), "IntegratedPDFs" = unlist( lb_med_ub.new ) )[ 1:3, , drop = FALSE ];

print( apply( results.table, 1:2, round ) )

###########################################
## The prior for the RNA threshold crossing time.
prior.over.rna.threshold.crossing.time <-
    agg_interpreter.new.in.list.with.inputs[[ "prior.over.last.biomarker.crossing.threshold.time" ]];

if( is.null( prior.over.rna.threshold.crossing.time ) ) {
    lb_med_ub.prior <- list( lb = range_start, med = ( range_start + ( ( range_end - range_start ) / 2 ) ), ub = range_end );
    cat( paste( "Prior is uniform over", range_start, "to", range_end ), fill = TRUE );
    prior.over.rna.threshold.crossing.time <- function( time ) { 1 }
} else {
    lb_med_ub.prior <- estimate_lb_med_ub_fn.new(fun = prior.over.rna.threshold.crossing.time,
                                range_start = range_start,
                                range_end = range_end,
                                           verbose = FALSE)
    cat( paste( "Prior", sapply( c( "lb", "med", "ub" ), function( stat.key ) { paste( stat.key, "=", sprintf( "%0.2f", lb_med_ub.prior[[ stat.key ]] ) ) } ), collapse = "\n" ), fill = TRUE );
}
results.table <- cbind( results.table, "IntegratedPDFs_PriorPDF" = unlist( lb_med_ub.prior )[ 1:3 ] );

##############
## Results: plot
the_plot.old <- 
plot_iihist(ihist = ihist, lb_med_ub = lb_med_ub.old, 
            range_start = as.numeric(as.Date('2019-02-01')),
            range_end = as.numeric(as.Date('2019-05-01')), 
            plot_aggregate = TRUE,
            produce_plot = FALSE,
            scales = scales,
            show_test_dates = 0.5,
                x_breaks =
                  seq( as.numeric(as.Date('2019-02-01')), as.numeric(as.Date('2019-05-01')), by = 2 ),
            aggregate_interpreter = agg_interpreter.old
            )
the_plot.old <- the_plot.old + 
  ggplot2::theme_grey(base_size = 9) + # 18 is a good number of non-vignette applications
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9, angle = 45, vjust = 1, hjust=1)) +
  ggplot2::guides(color = FALSE)

the_plot.new <- 
plot_iihist(ihist = ihist, lb_med_ub = lb_med_ub.new, 
            range_start = as.numeric(as.Date('2019-02-01')),
            range_end = as.numeric(as.Date('2019-05-01')), 
            plot_aggregate = TRUE,
            produce_plot = FALSE,
            scales = scales,
            show_test_dates = 0.5,
                x_breaks =
                  seq( as.numeric(as.Date('2019-02-01')), as.numeric(as.Date('2019-05-01')), by = 2 ),
            aggregate_interpreter = agg_interpreter.new
            )
the_plot.new <- the_plot.new + 
  ggplot2::theme_grey(base_size = 9) + # 18 is a good number of non-vignette applications
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9, angle = 45, vjust = 1, hjust=1)) +
  ggplot2::guides(color = FALSE)

######
## Prior plot
agg_interpreter.prior <- prior.over.rna.threshold.crossing.time;
ihist.RNA.tests <-
    agg_interpreter.new.in.list.with.inputs[[ "remaining.biomarker.bounds.to.process.last.to.first" ]][[ "RNA" ]]$test;
ihist.prior <-
    ihist[ !( ihist$test %in% ihist.RNA.tests ), , drop = FALSE ];
if( nrow( ihist.prior ) == 0 ) {
    ihist.prior <- NULL;
}

if( !is.null( ihist.prior ) ) {
    # Add the biomarker as a fifth column (character vector)
    ihist.prior.with.biomarker <- ihist.prior;
    ihist.prior.biomarker.classes <- sapply( all_assay_dynamics[ as.character( ihist.prior$test ) ], function( the.dynamics ) { as.character( the.dynamics$class ) } );
    ihist.prior.with.biomarker$biomarker <- ihist.prior.biomarker.classes;
        
    ihist.prior.assay.delay <- apply( ihist.prior.with.biomarker, 1, function( .row ) { get.time.from.biomarker.crossing.threshold.to.assay.positivity( all_assay_dynamics[[ .row[ "test" ] ]], the.biomarker.dynamics = biomarker_dynamics[[ .row[ "biomarker" ] ]] ) } );
    
    ## Add the adjusted date as a seventh column. This is the sample date minus the assay delay, so it reflects the constraints on the tie at which the biomarker crossed the canonical "1 copy per ml" level that are implied by the assay results.
    ihist.prior.with.biomarker$adjusted.date <- ihist.prior$sample_date - ihist.prior.assay.delay;
    
    ## Show the ITRIs rather than the window period CDFs.
    assay_interpreters.prior <-
        sapply( 1:nrow( ihist.prior ), function( i ) {
            the.biomarker <- all_assay_dynamics[[ ihist.prior[ i, "test" ] ]][[ "class" ]];
            assay.itri.dynamics <- biomarker_dynamics[[ the.biomarker ]];
            ## Apply the delay by using "adjusted.date".
            construct_assay_result_interpreter(assay_dynamics = assay.itri.dynamics,
                                               result = ihist.prior[i, 'result'],
                                               sample_date = ihist.prior.with.biomarker[i, 'adjusted.date'])
        } );
        
    the_plot.prior <- 
        plot_iihist(ihist = ihist.prior, lb_med_ub = lb_med_ub.prior,
                range_start = as.numeric(as.Date('2019-02-01')),
                range_end = as.numeric(as.Date('2019-05-01')), 
                plot_aggregate = TRUE,
                custom_aggregate_label = "Posterior distribution of\nRNA\nthreshold-crossing time",
                produce_plot = FALSE,
                scales = scales,
                show_test_dates = 0.5,
                    x_breaks =
                      seq( as.numeric(as.Date('2019-02-01')), as.numeric(as.Date('2019-05-01')), by = 2 ),
                assay_interpreters = assay_interpreters.prior,
                aggregate_interpreter = agg_interpreter.prior
                )
    the_plot.prior <- the_plot.prior + 
      ggplot2::theme_grey(base_size = 9) + # 18 is a good number of non-vignette applications
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9, angle = 45, vjust = 1, hjust=1)) +
          ggplot2::guides(color = FALSE)
} # End if !is.null( ihist.prior )

cat( paste( "Writing output to files beginning with", output.filename.base ), fill = TRUE );

write.csv( apply( results.table, 1:2, function( x ) { sprintf( "%0.2f", x ) } ), file = paste( output.filename.base, ".csv", sep = "" ) )

pdf( paste( output.filename.base, "_CDFs.pdf", sep = "" ) );
print(the_plot.old)
dev.off()

pdf( paste( output.filename.base, "_IntegratedPDFs.pdf", sep = "" ) );
print(the_plot.new)
dev.off()

if( !is.null( ihist.prior ) ) {
    pdf( paste( output.filename.base, "_IntegratedPDFs_PriorPDF.pdf", sep = "" ) );
    print(the_plot.prior)
    dev.off()
}

## Remove this scenario from the environment.
detach( scenario, character.only = TRUE );
