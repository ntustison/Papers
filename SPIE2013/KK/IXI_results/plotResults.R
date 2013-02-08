# resampler <- function( sp )
#   {
#   n <- nrow( sp.frame )
#   resample.rows <- sample( 1:n,size = n, replace = TRUE )
#   return( sp.frame[resample.rows,] )
#   }
#
# sp.spline.estimator <- function( x, y, m=300 )
#   {
#   # Fit spline to data, with cross-validation to pick lambda
#   fit <- smooth.spline( x, y, cv = TRUE );
#   # Set up a grid of m evenly-spaced points on which to evaluate the spline
#   eval.grid <- seq( from = min( x ), to = max( x ), length.out = m )
#     # Slightly inefficient to re-define the same grid every time we call this,
#     # but not a big overhead
#   # Do the prediction and return the predicted values
#   return( predict( fit, x = eval.grid)$y )  # We only want the predicted values
#   }
#
# sp.spline.cis <- function( B, alpha, m = 300 ) {
#   spline.main <- sp.spline.estimator( sp.frame, m = m )
#   # Draw B boottrap samples, fit the spline to each
#   spline.boots <- replicate( B, sp.spline.estimator( sp.resampler(), m = m ) )
#     # Result has m rows and B columns
#   cis.lower <- 2 * spline.main - apply( spline.boots, 1, quantile, probs = 1 - 0.5 * alpha )
#   cis.upper <- 2 * spline.main - apply( spline.boots, 1, quantile, probs = 0.5 * alpha )
#   return( list( main.curve = spline.main, lower.ci = cis.lower, upper.ci = cis.upper,
#               x = seq( from = min( sp.today ), to = max( sp.today ), length.out = m ) ) )
# }

calculate.spline <- function( spline.data, m = 300 )
  {
  fit <- smooth.spline( spline.data$x, spline.data$y, df = 3 );
  res <- residuals( fit );
  eval.grid <- seq( min( spline.data$x ), to = max( spline.data$x ), length.out = m );
  return( list( x = eval.grid, y = predict( fit, x = eval.grid )$y, res = res ) );
  }

resampler <- function( spline.data )
  {
  n <- length( spline.data$x );
  resample <- sample( 1:n, size = n, replace = TRUE );

  resampled.data <- list();
  resampled.data$x = spline.data$x[resample];
  resampled.data$y = spline.data$y[resample];

  return( resampled.data );
  }

calculate.spline.intervals <- function( spline.data, alpha, B, m = 300 )
  {
  spline.main <- calculate.spline( spline.data, m );
  spline.boots <- spline.data$y;

  numberofsamples <- length( spline.main$y );

  spline.boots <- matrix( 0, nrow = B, ncol = numberofsamples );

  for( i in seq( 1, B, by = 1 ) )
    {
    resampled.data <- resampler( spline.data );
    resampled.spline.main <- calculate.spline( resampled.data, m );
    for( j in seq( 1, numberofsamples, by = 1 ) )
      {
      spline.boots[i,j] <- resampled.spline.main$y[j];
      }
    }

  q_upper <- rep( 0, numberofsamples );
  q_lower <- rep( 0, numberofsamples );

  for( i in seq( 1, numberofsamples, by = 1 ) )
    {
    q_lower[i] <- as.numeric( quantile( spline.boots[,i], probs = 1 - 0.5 * alpha ) );
    q_upper[i] <- as.numeric( quantile( spline.boots[,i], probs = 0.5 * alpha ) );
    }

  return( list( lower = q_lower, upper = q_upper ) );
  }

male_color = "darkred";
female_color = "navyblue";

label_results <- read.csv( file = "labelresults.csv", header = TRUE );

male_label_results <- subset( label_results, SEX == 1 & AGE > 0 & AGE < 90 );
female_label_results <- subset( label_results, SEX == 2 & AGE > 0 & AGE < 90 );
both_label_results <- rbind( male_label_results, female_label_results );

male_label_results <- male_label_results[order( male_label_results$AGE ),];
female_label_results <- female_label_results[order( female_label_results$AGE ),];
both_label_results <- both_label_results[order( both_label_results$AGE ),];

cortical_labels <- c( "L occipital", "R occipital",
                     "L cingulate", "R cingulate",
                     "L insula",    "R insula",
                     "L temporal pole",   "R temporal pole",
                     "L superior temporal", "R superior temporal",
                     "L infero temporal", "R infero temporal",
                     "L parahippocampal", "R parahippocampal",
                     "L frontal pole",    "R frontal pole",
                     "L superior frontal","R superior frontal",
                     "L middle frontal",  "R middle frontal",
                     "L inferior",        "R inferior",
                     "L orbital frontal", "R orbital frontal",
                     "L precentral",      "R precentral",
                     "L superior parietal", "R superior parietal",
                     "L inferior parietal", "R inferior parietal",
                     "L postcentral",       "R postcentral" );

cortical_results <- read.csv( file = "corticalresults.csv", header = TRUE );

male_cortical_results <- subset( cortical_results, SEX == 1 & AGE > 0  & AGE < 90  & CORTICAL_THICKNESS > 0 );
female_cortical_results <- subset( cortical_results, SEX == 2 & AGE > 0 & AGE < 90 & CORTICAL_THICKNESS > 0  );
both_cortical_results <- rbind( male_cortical_results, female_cortical_results );

male_cortical_results <- male_cortical_results[order( male_cortical_results$AGE ),];
female_cortical_results <- female_cortical_results[order( female_cortical_results$AGE ),];
both_cortical_results <- both_cortical_results[order( both_cortical_results$AGE ),];

pdf( file = "whole_cortex_results.pdf", 8, 6 );

plot( c( 20, 90 ), c( 0.5, 5.5 ), type = "n",
  xlab = expression( paste( bold( "Age (years)" ) ) ),
  ylab = expression( paste( bold( "Thickness (mm)" ) ) ),
  frame.plot = FALSE,  );
title( main = "Cortical thickness results", col.main = "black", font.main = 2 );

male_data <- list();
male_data$x <- male_cortical_results$AGE;
male_data$y <- male_cortical_results$CORTICAL_THICKNESS;
male_data$sampled.spline <- calculate.spline( male_data, 300 );
male_data$int <- calculate.spline.intervals( male_data, 0.05, 1000, 300 );

points( male_data, pch = 22, lwd = 0.5, cex = 0.75, col = male_color );
lines( male_data$sampled.spline, lty=1, col = male_color, lwd = 5 );
lines( male_data$sampled.spline$x, male_data$int$lower, lty=5, col = male_color, lwd = 3 );
lines( male_data$sampled.spline$x, male_data$int$upper, lty=5, col = male_color, lwd = 3 );

female_data <- list();
female_data$x <- female_cortical_results$AGE;
female_data$y <- female_cortical_results$CORTICAL_THICKNESS;
female_data$sampled.spline <- calculate.spline( female_data, 300 );
female_data$int <- calculate.spline.intervals( female_data, 0.05, 1000, 300 );

points( female_data, pch = 23, lwd = 0.5, cex = 0.75, col = female_color );
lines( female_data$sampled.spline, lty=1, col = female_color, lwd = 5 );
lines( female_data$sampled.spline$x, female_data$int$lower, lty=5, col = female_color, lwd = 3 );
lines( female_data$sampled.spline$x, female_data$int$upper, lty=5, col = female_color, lwd = 3 );

legend( "topright", c( "Female", "Male" ),
   col = c( female_color, male_color ),
 pch = c( 23, 22 ),
 title = expression( paste( bold( "Gender" ) ) ),
# box.lwd = 2,
  bty = "n",
  lwd = 2 );
axis( 1, tick = TRUE, lwd = 2, font = 2 );
axis( 2, tick = TRUE, lwd = 2, font = 2 );
dev.off();

both_data <- list()
both_data$x <- both_cortical_results$AGE;
both_data$y <- both_cortical_results$CORTICAL_THICKNESS;
both_data$sampled.spline <- calculate.spline( both_data, 300 );
# both_data$int <- calculate.spline.intervals( both_data, 0.05, 1000, 300 );

# perform hypothesis testing
null.hypothesis.ss <- sum( ( both_data$sampled.spline$res )^2 );
alt.hypothesis.ss <- sum( ( male_data$sampled.spline$res )^2 ) + sum( ( female_data$sampled.spline$res )^2 );
diff.ss = ( null.hypothesis.ss - alt.hypothesis.ss );
rel.diff.ss = ( diff.ss ) / alt.hypothesis.ss;

null.hypothesis.df <- length( both_data$sampled.spline$res ) - 3;
alt.hypothesis.df <- length( male_data$sampled.spline$res ) + length( female_data$sampled.spline$res ) - 2 * 3;
diff.df <- ( null.hypothesis.df - alt.hypothesis.df );
rel.diff.df = ( diff.df ) / alt.hypothesis.df;

fstat <- rel.diff.ss / rel.diff.df;
r <- pf( fstat, diff.df, alt.hypothesis.df );
cat( "Whole cortex p-value: ",  2 * min( r, 1 - r ), "\n", sep = "" );

print( "========================================================" );
print( "           gender differences                           " );
print( "========================================================" );

for( i in seq( 1, 32, by = 1 ) )
  {
  pdf( file = paste( "label", i, "_results.pdf", sep = "" ), 8, 6 );

  plot( c( 20, 90 ), c( 0.5, 5.5 ), type = "n",
    xlab = expression( paste( bold( "Age (years)" ) ) ),
    ylab = expression( paste( bold( "Thickness (mm)" ) ) ),
    frame.plot = FALSE  );
  title( main = paste( "Cortical thickness results (", cortical_labels[i], ")", sep = "" ), col.main = "black", font.main = 2 );

  male_data <- list();
  male_data$x <- male_label_results$AGE;
  male_data$y <- male_label_results[,i+3];
  male_data$sampled.spline <- calculate.spline( male_data, 300 );
  male_data$int <- calculate.spline.intervals( male_data, 0.05, 1000, 300 );

  points( male_data, pch = 22, lwd = 0.5, cex = 0.75, col = male_color );
  lines( male_data$sampled.spline, lty=1, col = male_color, lwd = 5 );
  lines( male_data$sampled.spline$x, male_data$int$lower, lty=5, col = male_color, lwd = 3 );
  lines( male_data$sampled.spline$x, male_data$int$upper, lty=5, col = male_color, lwd = 3 );

  female_data <- list();
  female_data$x <- female_label_results$AGE;
  female_data$y <- female_label_results[,i+3];
  female_data$sampled.spline <- calculate.spline( female_data, 300 );
  female_data$int <- calculate.spline.intervals( female_data, 0.05, 1000, 300 );

  points( female_data, pch = 23, lwd = 0.5, cex = 0.75, col = female_color );
  lines( female_data$sampled.spline, lty=1, col = female_color, lwd = 5 );
  lines( female_data$sampled.spline$x, female_data$int$lower, lty=5, col = female_color, lwd = 3 );
  lines( female_data$sampled.spline$x, female_data$int$upper, lty=5, col = female_color, lwd = 3 );

  legend( "topright", c( "Female", "Male" ),
     col = c( female_color, male_color ),
   pch = c( 23, 22 ),
   title = expression( paste( bold( "Gender" ) ) ),
  # box.lwd = 2,
    bty = "n",
    lwd = 2 );
  axis( 1, tick = TRUE, lwd = 2, font = 2 );
  axis( 2, tick = TRUE, lwd = 2, font = 2 );
#
  dev.off();

  # hypothesis testing

#   both_data <- list()
#   both_data$x <- both_label_results$AGE;
#   both_data$y <- both_label_results[,i+3];
#   both_data$sampled.spline <- calculate.spline( both_data, 300 );
#   # both_data$int <- calculate.spline.intervals( both_data, 0.05, 1000, 300 );
#
#   # perform hypothesis testing
#   null.hypothesis.ss <- sum( ( both_data$sampled.spline$res )^2 );
#   alt.hypothesis.ss <- sum( ( male_data$sampled.spline$res )^2 ) + sum( ( female_data$sampled.spline$res )^2 );
#   diff.ss = ( null.hypothesis.ss - alt.hypothesis.ss );
#   rel.diff.ss = ( diff.ss ) / alt.hypothesis.ss;
#
#   null.hypothesis.df <- length( both_data$sampled.spline$res ) - 3;
#   alt.hypothesis.df <- length( male_data$sampled.spline$res ) + length( female_data$sampled.spline$res ) - 2 * 3;
#   diff.df <- ( null.hypothesis.df - alt.hypothesis.df );
#   rel.diff.df = ( diff.df ) / alt.hypothesis.df;
#
#   fstat <- rel.diff.ss / rel.diff.df;
#   r <- pf( fstat, diff.df, alt.hypothesis.df );
#   cat( cortical_labels[i], " p-value: ",  2 * min( r, 1 - r ), "\n", sep = "" );

  male.20To40 <- male_data$y[male_data$x >= 20 & male_data$x <= 40];
  male.40To60 <- male_data$y[male_data$x >= 40 & male_data$x <= 60];
  male.60plus <- male_data$y[male_data$x >= 60];

  male.20To40.mean <- mean( male.20To40 );
  male.20To40.sd <- sd( male.20To40 );
  male.40To60.mean <- mean( male.40To60 );
  male.40To60.sd <- sd( male.40To60 );
  male.60plus.mean <- mean( male.60plus );
  male.60plus.sd <- sd( male.60plus );

  female.20To40 <- female_data$y[female_data$x >= 20 & female_data$x <= 40];
  female.40To60 <- female_data$y[female_data$x >= 40 & female_data$x <= 60];
  female.60plus <- female_data$y[female_data$x >= 60];

  female.20To40.mean <- mean( female.20To40 );
  female.20To40.sd <- sd( female.20To40 );
  female.40To60.mean <- mean( female.40To60 );
  female.40To60.sd <- sd( female.40To60 );
  female.60plus.mean <- mean( female.60plus );
  female.60plus.sd <- sd( female.60plus );

  cat( cortical_labels[i], " & ",
    female.20To40.mean, ' \\pm ', 1.96 * female.20To40.sd, " & ",
    male.20To40.mean, ' \\pm ', 1.96 * male.20To40.sd, " & ",
    female.40To60.mean, ' \\pm ', 1.96 * female.40To60.sd, " & ",
    male.40To60.mean, ' \\pm ', 1.96 * male.40To60.sd, " & ",
    female.60plus.mean, ' \\pm ', 1.96 * female.60plus.sd, " & ",
    male.60plus.mean, ' \\pm ', 1.96 * male.60plus.sd, " \\\\",
    "\n", sep = "" );

  }

print( "========================================================" );
print( "           male hemispherical differences               " );
print( "========================================================" );

for( i in seq( 1, 32, by = 2 ) )
  {
  left_data <- list();
  left_data$x <- male_label_results$AGE;
  left_data$y <- male_label_results[,i+3];
  left_data$sampled.spline <- calculate.spline( left_data, 300 );
#   left_data$int <- calculate.spline.intervals( left_data, 0.05, 1000, 300 );

  right_data <- list();
  right_data$x <- male_label_results$AGE;
  right_data$y <- male_label_results[,i+1+3];
  right_data$sampled.spline <- calculate.spline( right_data, 300 );
#   right_data$int <- calculate.spline.intervals( right_data, 0.05, 1000, 300 );

  both_data <- list();
  both_data$x <- c( left_data$x, right_data$x );
  both_data$y <- c( left_data$y, right_data$y );
  both_data <- both_data[order( both_data$x )];

  both_data$sampled.spline <- calculate.spline( both_data, 300 );
  # both_data$int <- calculate.spline.intervals( both_data, 0.05, 1000, 300 );

  # perform hypothesis testing
  null.hypothesis.ss <- sum( ( both_data$sampled.spline$res )^2 );
  alt.hypothesis.ss <- sum( ( left_data$sampled.spline$res )^2 ) + sum( ( right_data$sampled.spline$res )^2 );
  diff.ss = ( null.hypothesis.ss - alt.hypothesis.ss );
  rel.diff.ss = ( diff.ss ) / alt.hypothesis.ss;

  null.hypothesis.df <- length( both_data$sampled.spline$res ) - 3;
  alt.hypothesis.df <- length( left_data$sampled.spline$res ) + length( right_data$sampled.spline$res ) - 2 * 3;
  diff.df <- ( null.hypothesis.df - alt.hypothesis.df );
  rel.diff.df = ( diff.df ) / alt.hypothesis.df;

  fstat <- rel.diff.ss / rel.diff.df;
  r <- pf( fstat, diff.df, alt.hypothesis.df );
  cat( cortical_labels[i], " p-value: ",  2 * min( r, 1 - r ), "\n", sep = "" );
  }


print( "========================================================" );
print( "         female hemispherical differences               " );
print( "========================================================" );

for( i in seq( 1, 32, by = 2 ) )
  {
  left_data <- list();
  left_data$x <- female_label_results$AGE;
  left_data$y <- female_label_results[,i+3];
  left_data$sampled.spline <- calculate.spline( left_data, 300 );
#   left_data$int <- calculate.spline.intervals( left_data, 0.05, 1000, 300 );

  right_data <- list();
  right_data$x <- female_label_results$AGE;
  right_data$y <- female_label_results[,i+1+3];
  right_data$sampled.spline <- calculate.spline( right_data, 300 );
#   right_data$int <- calculate.spline.intervals( right_data, 0.05, 1000, 300 );

#   both_data <- list();
#   both_data$x <- c( left_data$x, right_data$x );
#   both_data$y <- c( left_data$y, right_data$y );
#   both_data <- both_data[order( both_data$x )];
#
#   both_data$sampled.spline <- calculate.spline( both_data, 300 );
#   # both_data$int <- calculate.spline.intervals( both_data, 0.05, 1000, 300 );
#
#   # perform hypothesis testing
#   null.hypothesis.ss <- sum( ( both_data$sampled.spline$res )^2 );
#   alt.hypothesis.ss <- sum( ( left_data$sampled.spline$res )^2 ) + sum( ( right_data$sampled.spline$res )^2 );
#   diff.ss = ( null.hypothesis.ss - alt.hypothesis.ss );
#   rel.diff.ss = ( diff.ss ) / alt.hypothesis.ss;
#
#   null.hypothesis.df <- length( both_data$sampled.spline$res ) - 3;
#   alt.hypothesis.df <- length( left_data$sampled.spline$res ) + length( right_data$sampled.spline$res ) - 2 * 3;
#   diff.df <- ( null.hypothesis.df - alt.hypothesis.df );
#   rel.diff.df = ( diff.df ) / alt.hypothesis.df;
#
#   fstat <- rel.diff.ss / rel.diff.df;
#   r <- pf( fstat, diff.df, alt.hypothesis.df );
#   cat( cortical_labels[i], " p-value: ",  2 * min( r, 1 - r ), "\n", sep = "" );

  tt <- t.test( left_data$y, right_data$y, paired = TRUE, alternative = "two.sided", var.equal = false );
  cat( cortical_labels[i], " p-value: ", tt$p.value, "\n", sep = "" );
  }
