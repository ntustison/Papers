library( MASS )
library( glmulti )
library( caret )

rm( list = ls( all = TRUE ) )

meandiffusion <- read.csv( 'meandiffusion.csv' )
volumes <- read.csv( 'volume.csv' )
thickness <- read.csv( 'thickness.csv' )
demographics <- read.csv( 'demographics.csv' )

mylabels <- c( "whole", "L occipital", "R occipital",
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


wholeThickness <- rowSums( as.matrix( volumes ) * as.matrix( thickness ) ) / rowSums( volumes );
wholeMeanDiffusion <- rowSums( as.matrix( volumes ) * as.matrix( meandiffusion ) ) / rowSums( volumes );

accuracy <- rep( 0, 33 );
for( n in 0:100 )
  {
  cat( "Run ", n, "\n", sep = "" );
  trainingIndices <- createDataPartition( demographics$AGE, p = 1/2, list = FALSE )

  training.data <- data.frame( sex = demographics$SEX[trainingIndices],
                               age = demographics$AGE[trainingIndices],
                               wth = wholeThickness[trainingIndices],
                               th = thickness[trainingIndices,],
                               wmd = wholeMeanDiffusion[trainingIndices],
                               md = meandiffusion[trainingIndices,]
                             );

  testing.data <- data.frame( sex = demographics$SEX[-trainingIndices],
                              age = demographics$AGE[-trainingIndices],
                              wth = wholeThickness[-trainingIndices],
                              th = thickness[-trainingIndices,],
                              wmd = wholeMeanDiffusion[-trainingIndices],
                              md = meandiffusion[-trainingIndices,]
                            );


  models <- list();
  count <- 1;

#   models[[count]] <- rlm( age ~ sex + wth + wmd + I( wth^0.5 ) + I( wmd^0.5 ) + wth:wmd, data = training.data, maxit = 100 );
#   models[[count]] <- lm( age ~ sex + wth + wmd + I( wth^0.5 ) + I( wmd^0.5 ) + wth:wmd, data = training.data );
  models[[count]] <- rlm( age ~ sex + wth + wmd, data = training.data );
  count <- count + 1;

  for( i in 1:32 )
    {
#     formulaString <- paste( "age ~ sex + ", names( testing.data )[3+i],  " + ", names( testing.data )[3+33+i], " + ",
#                             names( testing.data )[3+i],  ":", names( testing.data )[3+33+i], " + ",
#                             "I(", names( testing.data )[3+i], "^0.5) + I(", names( testing.data )[3+33+i], "^0.5)", sep = "" );
    formulaString <- paste( "age ~ sex + ", names( testing.data )[3+i],  " + ", names( testing.data )[3+33+i], sep = "" );
    models[[count]] <- rlm( as.formula( formulaString ), data = training.data, maxit = 100 )
    count <- count + 1;
    }


  for( i in 1:length( models ) )
    {
    age.real <- testing.data$age;
    age.predicted <- predict( models[[i]], newdata = testing.data );

    accuracy[i] <- accuracy[i] + mean( abs( age.real - age.predicted ), na.rm = TRUE );
    }
  }

orderedAccuracy <- sort( accuracy, index.return = TRUE );

for( i in 1:length( accuracy ) )
  {
  cat( "accuracy (", mylabels[orderedAccuracy$ix[i]], "): ",  orderedAccuracy$x[i]/100, "\n", sep = "" );
  }



# pdf( file = "glm_model.pdf", 8, 8 );
# plot( c( 20, 90 ), c( 20, 90 ), type = "n",
#     xlab = expression( paste( bold( "Real Age (years)" ) ) ),
#     ylab = expression( paste( bold( "Predicted Age (years)" ) ) ),
#     frame.plot = FALSE  );
# lines( c( 20, 90 ), c( 20, 90 ), lwd = 3, lty = 3, col = "black" );
# title( main = paste( "Cortical thickness prediction of age"), col.main = "black", font.main = 2 );
# axis( 1, tick = TRUE, lwd = 2, font = 2 );
# axis( 2, tick = TRUE, lwd = 2, font = 2 );
#
# age.real <- rep( 0, nrow( testing.data ) );
# age.predicted <- rep( 0, nrow( testing.data ) );
# age.correct <- rep( 1, nrow( testing.data ) );
#
# points( age.real[age.correct == 0], age.predicted[age.correct == 0], pch = 21, lwd = 1.5, cex = 1.0, col = "darkred" );
# points( age.real[age.correct == 1], age.predicted[age.correct == 1], pch = 23, lwd = 1.5, cex = 1.0, col = "navyblue" );
# cat( "Accuracy: ", mean( abs( age.real - age.predicted ) ), "\n", sep = "" );
#
# dev.off();
