library( "earth" )


get.used.pred.names <- function( obj ) # obj is an earth object
{
  any1 <- function(x) any(x != 0)    # like any but no warning if x is double
  names(which(apply(obj$dirs[obj$selected.terms, , drop=FALSE], 2, any1)))
}

rm( list = ls( all = TRUE ) )

meandiffusion <- read.csv( 'meandiffusion.csv', header = TRUE )
volumes <- read.csv( 'volume.csv', header = TRUE )
thickness <- read.csv( 'thickness.csv', header = TRUE )
demographics <- read.csv( 'demographics.csv', header = TRUE )

meandiffusion <- meandiffusion[demographics$SEX == 2,];
volumes <- volumes[demographics$SEX == 2,];
thickness <- thickness[demographics$SEX == 2,];
demographics <- demographics[demographics$SEX == 2,];


agePredictionData <- data.frame( md = meandiffusion, th = thickness );

trainingIndices <- createDataPartition( demographics$AGE, p = 1/2, list = FALSE );

agePredictionData.training <- agePredictionData[trainingIndices,];
agePredictionData.testing <- agePredictionData[-trainingIndices,]

fit <- earth( x = agePredictionData.training, y = demographics$AGE[trainingIndices], degree = 2, pmethod = "exhaustive"  );

age.real <- demographics$AGE[-trainingIndices];
age.predict <- predict( fit, agePredictionData.testing );

pdf( file = "brian_results_model.pdf", 8, 8 );
plot( c( 20, 90 ), c( 20, 90 ), type = "n",
    xlab = expression( paste( bold( "Real Age (years)" ) ) ),
    ylab = expression( paste( bold( "Predicted Age (years)" ) ) ),
    frame.plot = FALSE  );
lines( c( 20, 90 ), c( 20, 90 ), lwd = 3, lty = 3, col = "black" );
title( main = paste( "Cortical thickness prediction of age"), col.main = "black", font.main = 2 );
axis( 1, tick = TRUE, lwd = 2, font = 2 );
axis( 2, tick = TRUE, lwd = 2, font = 2 );
points( x = age.real, y = age.predict, pch = 16, lwd = 1.5, cex = 1.0, col = "darkred" );
dev.off();

print( mean( abs( age.real - age.predict ), na.rm = TRUE ) );
