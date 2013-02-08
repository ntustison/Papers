library( MASS )
library( bestglm )
library( caret )

meandiffusion <- read.csv( 'meandiffusion.csv' )
volumes <- read.csv( 'volume.csv' )
thickness <- read.csv( 'thickness.csv' )
demographics <- read.csv( 'demographics.csv' )

# create the total model

trainingIndices <- createDataPartition( demographics$AGE, p = 1/2, list = FALSE );

training.data <- data.frame( th = thickness[trainingIndices,],
#                              md = meandiffusion[trainingIndices,],
                             sex = demographics$SEX[trainingIndices],
                             age = demographics$AGE[trainingIndices] );
testing.data <- data.frame( th = thickness[-trainingIndices,],
#                             md = meandiffusion[-trainingIndices,],
                             sex = demographics$SEX[-trainingIndices],
                             age = demographics$AGE[-trainingIndices] );

model.total <- bestglm( training.data, IC = "CV" );

# create the three sub-models

subjects <- list();
subjects[[1]] <- ( demographics$AGE[trainingIndices] >= 20 & demographics$AGE[trainingIndices] <= 40 );
subjects[[2]] <- ( demographics$AGE[trainingIndices] >= 40 & demographics$AGE[trainingIndices] <= 60 );
subjects[[3]] <- ( demographics$AGE[trainingIndices] >= 60 );

models.sub <- list();
for( i in seq( 1, 3, by = 1 ) )
  {
  cat( "constructing model ", i, "\n", sep = "" );
  dem <- demographics[subjects[[i]],];
  vol <- volumes[subjects[[i]],];
  th <- thickness[subjects[[i]],];
  md <- meandiffusion[subjects[[i]],];

  training.data <- data.frame( th = th,
#                              md = md,
                               sex = dem$SEX,
                               age = dem$AGE );
  models.sub[[i]] <- bestglm( training.data, IC = "CV", t = 125 );
  }

pdf( file = "brian_results_model.pdf", 8, 8 );
plot( c( 20, 90 ), c( 20, 90 ), type = "n",
    xlab = expression( paste( bold( "Real Age (years)" ) ) ),
    ylab = expression( paste( bold( "Predicted Age (years)" ) ) ),
    frame.plot = FALSE  );
lines( c( 20, 90 ), c( 20, 90 ), lwd = 3, lty = 3, col = "black" );
title( main = paste( "Cortical thickness prediction of age"), col.main = "black", font.main = 2 );
axis( 1, tick = TRUE, lwd = 2, font = 2 );
axis( 2, tick = TRUE, lwd = 2, font = 2 );

age.real <- rep( 0, nrow( testing.data ) );
age.predicted <- rep( 0, nrow( testing.data ) );
age.correct <- rep( 1, nrow( testing.data ) );

for( i in seq( 1, nrow( testing.data ), by = 1 ) )
  {
  age.real[i] <- testing.data[i,]$age;

  age.predicted[i] <- predict( model.total$BestModel, newdata = testing.data[i,] );

  whichsubmodel <- 1;
  if( age.predicted[i] >= 40 && age.predicted[i] <= 60 )
    {
    whichsubmodel <- 2;
    }
  else if( age.predicted[i] >= 60 )
    {
    whichsubmodel <- 3;
    }

  age.predicted[i] <- predict( models.sub[[whichsubmodel]]$BestModel, newdata = testing.data[i,] );

  if( ( whichsubmodel == 1 && age.real[i] > 40 ) ||
      ( whichsubmodel == 2 && ( age.real[i] < 40 || age.real[i] > 60 ) ) ||
      ( whichsubmodel == 3 && age.real[i] < 60 ) )
    {
    age.correct[i] <- 0;
    }
  cat( age.real[i], " <--> ", age.predicted[i], " = ", age.correct[i], "\n", sep = "" );
  }
points( age.real[age.correct == 0], age.predicted[age.correct == 0], pch = 21, lwd = 1.5, cex = 1.0, col = "darkred" );
points( age.real[age.correct == 1], age.predicted[age.correct == 1], pch = 23, lwd = 1.5, cex = 1.0, col = "navyblue" );
cat( "Accuracy: ", mean( abs( age.real - age.predicted ) ), "\n", sep = "" );

dev.off();
