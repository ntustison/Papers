syn_files <- list();

# syn_files[[1]] <- list.files( "./output_syn_300_0/", full.names = TRUE, pattern = "*.csv" );
# syn_files[[2]] <- list.files( "./output_syn_300_0_downsample/", full.names = TRUE, pattern = "*.csv" );
# syn_files[[3]] <- list.files( "./output_syn_500_0/", full.names = TRUE, pattern = "*.csv" );
# syn_files[[4]] <- list.files( "./output_syn_500_0_downsample/", full.names = TRUE, pattern = "*.csv" );
# syn_files[[5]] <- list.files( "./output_syn_700_0/", full.names = TRUE, pattern = "*.csv" );
# syn_files[[6]] <- list.files( "./output_syn_700_0_downsample/", full.names = TRUE, pattern = "*.csv" );
# syn_files[[7]] <- list.files( "./output_bsplinesyn_10x12x10_order2/", full.names = TRUE, pattern = "*.csv" );
# syn_files[[8]] <- list.files( "./output_bsplinesyn_10x12x10_order2_downsample/", full.names = TRUE, pattern = "*.csv" );
# syn_files[[9]] <- list.files( "./output_bsplinesyn_10x12x10_order3/", full.names = TRUE, pattern = "*.csv" );
# syn_files[[10]] <- list.files( "./output_bsplinesyn_10x12x10_order3_downsample/", full.names = TRUE, pattern = "*.csv" );
#
# syn_files[[11]] <- list.files( "./output_gauss_300_0/", full.names = TRUE, pattern = "*.csv" );
# syn_files[[12]] <- list.files( "./output_dmffd_10x12x10_order2/", full.names = TRUE, pattern = "*.csv" );
# syn_files[[13]] <- list.files( "./output_dmffd_10x12x10_order3/", full.names = TRUE, pattern = "*.csv" );
# syn_files[[14]] <- list.files( "./output_original_SyN/", full.names = TRUE, pattern = "*.csv" );
#
#
# legendEntries <- c( "SyN[3]", "SyNd[3]", "SyN[5]", "SyNd[5]", "SyN[7]", "SyNd[7]",
#                     "BSyN[10x12x10,2]", "BSyNd[10x12x10,2]", "BSyN[10x12x10,3]", "BSyNd[10x12x10,3]",
#                     "Gauss[0.5,3]", "DMFFD[0.5,10x12x10,2]", "DMFFD[0.5,10x12x10,3]", "OSyN[3]"
#                   );
# colors <- c( "red", "orange", "yellow", "green", "blue", "purple", "brown", "black", "pink", "cyan", "magenta", "aquamarine", "azure4", "chocolate" );


syn_files[[1]] <- list.files( "./output_syn_300_0_downsample/", full.names = TRUE, pattern = "*.csv" );
syn_files[[2]] <- list.files( "./output_syn_500_0_downsample/", full.names = TRUE, pattern = "*.csv" );
syn_files[[3]] <- list.files( "./output_bsplinesyn_10x12x10_order3_downsample/", full.names = TRUE, pattern = "*.csv" );
syn_files[[4]] <- list.files( "./output_original_SyN/", full.names = TRUE, pattern = "*.csv" );

legendEntries <- c( "SyN[3]", "SyN[5]", "BSyN[10x12x10,3]", "OriginalSyN[3]" );

colors <- c( "red", "green", "blue", "orange" );

times <- list();
metrics0 <- list();
metrics1 <- list();
metrics2 <- list();
metrics3 <- list();

metrics <- list();

labelOverlaps <- list();

for( i in 1:length( syn_files ) )
  {
  print( syn_files[[i]][1] );

  labelOverlaps[[i]] <- read.csv( syn_files[[i]][2], header = FALSE );

  if( i < length( syn_files ) )
    {
    times[[i]] <- read.csv( syn_files[[i]][1], header = FALSE );
    metrics0[[i]] <- read.csv( syn_files[[i]][3], header = FALSE );
    metrics1[[i]] <- read.csv( syn_files[[i]][4], header = FALSE );
    metrics2[[i]] <- read.csv( syn_files[[i]][5], header = FALSE );
    metrics3[[i]] <- read.csv( syn_files[[i]][6], header = FALSE );

    metrics[[i]] <- cbind( metrics0[[i]], metrics1[[i]], metrics2[[i]], metrics3[[i]] );
    }
  }

startIterations <- c( ncol( metrics0[[1]] ), ncol( metrics1[[1]] ),
  ncol( metrics2[[1]] ), ncol( metrics3[[1]] ) );
totalNumberOfIterations <- sum( startIterations );
iterations <- 1:totalNumberOfIterations;

pchNumbers <- c( 19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31 );

# plot the metric evolution over the iterations per level


##################################################

pdf( file = "/Users/ntustison/Desktop/LabelOverlapNIREP.pdf", 8, 6 );

plot( c( 1, nrow( labelOverlaps[[1]] ) ), c( 0.5, 0.85 ), type = "n",
  xlab = expression( paste( bold( "Label" ) ) ),
  ylab = expression( paste( bold( "Target Overlap" ) ) ),
  frame.plot = FALSE );
# title( main = "NIREP Evaluation", col.main = "black", font.main = 2 );

for( i in 1:length( labelOverlaps ) )
  {
  numberOfLabels <- nrow( labelOverlaps[[i]] );
  lines( 1:numberOfLabels, labelOverlaps[[i]][,2],, col = colors[i], lty = 1, lwd = 2 )
  }

legend( "bottomleft", legendEntries,
#   lty = 1:5,
  col = colors[1:length(syn_files)],
#   title = expression( paste( bold( "Canine Subject" ) ) ),
  bty = "n",
  lwd = 1.5 );
axis( 1, tick = TRUE, lwd = 2, font = 2 );
axis( 2, tick = TRUE, lwd = 2, font = 2 );

dev.off();

#################################################################


pdf( file = "/Users/ntustison/Desktop/MetricEvolution.pdf", 4, 5 );
#
plot( c( 1, 300 ), c( -0.975, -0.81 ), type = "n",
  xlab = expression( paste( bold( "Iteration" ) ) ),
  ylab = expression( paste( bold( "CC metric value" ) ) ),
  frame.plot = FALSE );
title( main = "Metric evolution", col.main = "black", font.main = 2 );
#

for( i in 1:length( metrics ) )
  {
  for( j in 1:4 )
    {
    startIteration <- 1;
    if( j > 1 )
      {
      for( k in 1:(j-1) )
        {
        startIteration <- startIteration + startIterations[k];
        }
      }
    endIteration <- startIteration + startIterations[j] - 1;
#
    print( startIteration );
    print( endIteration );
#
#
    mean_metric <- c();
    n_metric <- c();
    for( j in startIteration:endIteration )
      {
      index <- j - startIteration + 1;
      mean_metric[index] <- 0;
      n_metric[index] <- 0;
      for( k in 1:nrow( metrics[[i]] ) )
        {
        if( metrics[[i]][k,j] != 0 )
          {
          mean_metric[index] <- mean_metric[index] + metrics[[i]][k,j];
          n_metric[index] <- n_metric[index] + 1;
          }
        }
      }

    for( j in 1:length( n_metric ) )
      {
      if( n_metric[j] > 5 )
        {
        mean_metric[j] <- mean_metric[j] / n_metric[j];
        }
      }
    mean_metric <- mean_metric[ which( n_metric > 5 ) ];
#
    lines( seq( startIteration, startIteration + length( mean_metric ) - 1, by = 1 ), mean_metric, col = colors[i], lty = 1, lwd = 3 );
#
#
    lines( c( startIteration, startIteration ), c( -0.975, -0.81 ), col = "black", lty = 2, lwd = 1 );
#
    }
  }
# legend( "topright", legendEntries[1:3],
# #   lty = 1:5,
#   col = colors[1:length(syn_files)],
# #   title = expression( paste( bold( "Canine Subject" ) ) ),
#   bty = "n",
#   lwd = 1.5 );
axis( 1, tick = TRUE, lwd = 2, font = 2 );
axis( 2, tick = TRUE, lwd = 2, font = 2 );
#
dev.off();

##################################################

pdf( file = "/Users/ntustison/Desktop/Convergences.pdf", 4, 5 );

plot( c( 1, 300 ), c( 0, 1 ), type = "n",
  xlab = expression( paste( bold( "Iteration" ) ) ),
  ylab = expression( paste( bold( "% converged" ) ) ),
  frame.plot = FALSE );
title( main = "Convergence", col.main = "black", font.main = 2 );

for( i in 1:length( metrics ) )
  {
  for( j in 1:4 )
    {
    startIteration <- 1;
    if( j > 1 )
      {
      for( k in 1:(j-1) )
        {
        startIteration <- startIteration + startIterations[k];
        }
      }
    endIteration <- startIteration + startIterations[j] - 1;
    n_metric <- c();
    for( j in startIteration:endIteration )
      {
      index <- j - startIteration + 1;
      n_metric[index] <- 0;
      for( k in 1:nrow( metrics[[i]] ) )
        {
        if( metrics[[i]][k,j] != 0 )
          {
          n_metric[index] <- n_metric[index] + 1;
          }
        }
      }
    lines( seq( startIteration, startIteration + length( n_metric ) - 1, by = 1 ), 1.0 - n_metric / 240, col = colors[i], lty = 1, lwd = 3 );
    lines( c( startIteration, startIteration ), c( 0, 1 ), col = "black", lty = 2, lwd = 1 );
    }
  }
legend( "topright", legendEntries[1:3],
#   lty = 1:5,
  col = colors[1:length(syn_files)],
#   title = expression( paste( bold( "Canine Subject" ) ) ),
  bty = "n",
  lwd = 1.5 );
axis( 1, tick = TRUE, lwd = 2, font = 2 );
axis( 2, tick = TRUE, lwd = 2, font = 2 );
#
dev.off();

##################################################

pdf( file = "/Users/ntustison/Desktop/TimesNIREP.pdf", 8, 5 );

meanTimes <- c();
sdTimes <- c();
for( i in 1:length( times ) )
  {
  meanTimes[i] <- mean( times[[i]]/3600 );
  sdTimes[i] <- sd( times[[i]]/3600 );
  }

error.bar <- function(x, y, upper, lower=upper, length=0.1,...) {
  if( length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper) )
    {
    stop("vectors must be same length")
    }
  arrows( x, y+upper, x, y-lower, angle=90, code=3, length=length, ...)
  }



barx <- barplot( meanTimes, names.arg = 1:length( times ),
  ylab = expression( paste( bold( "Elapsed times (hours)" ) ) ),
  ylim = c( 0, 20 )
   );
error.bar( barx, meanTimes, 1.96 * sdTimes / sqrt( length( meanTimes ) ) );
# title( main = "SyN NIREP Evaluation", col.main = "black", font.main = 2 );

dev.off();

