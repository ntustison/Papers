library( kernlab )
library( ggplot2 )

trainingData <- read.csv( "trainingBrainSegmentationPosteriors2Projections.csv" )
testingData <- read.csv( "testingBrainSegmentationPosteriors2Projections.csv" )
# trainingData <- read.csv( "trainingCorticalThicknessProjections.csv" )
# testingData <- read.csv( "testingCorticalThicknessProjections.csv" )

# Remove ID, gender, and site

drops <- c( "ID", "SEX", "SITE" )

trainingData <- trainingData[, !( names( trainingData ) %in% drops )]
testingData <- testingData[, !( names( trainingData ) %in% drops )]

brainAgeVM <- rvm( AGE ~ ., data = trainingData, type = "regression",
                    kernel = "polydot", kpar = list( degree = 1 ),
                    verbosity = 1, tol = .Machine$double.eps,
                    minmaxdiff = 1e-3, cross = 0, fit = TRUE,
                    na.action = na.omit )
predictedAge <- predict( brainAgeVM, testingData )

correlation <- cor( testingData$AGE, predictedAge )
cat( "Correlation (true vs. predicted age): ", correlation, "\n", sep = '' )


plotData <- data.frame( TrueAge = testingData$AGE, PredictedAge = predictedAge )

brainAgeRegression <- lm( testingData$AGE ~ 1 + predictedAge )



brainAgePlot <- ggplot( plotData, aes( x = TrueAge, y = PredictedAge ) ) +
                stat_smooth( colour = "navyblue", formula = y ~ 1 + x, method = "lm",
                  size = 1, n = 1000, level = 0.95, se = TRUE, fullrange = TRUE, fill = 'black', alpha = 0.25 ) +
                geom_point( colour = "darkred", size = 4, alpha = 0.75 ) +
                scale_x_continuous( "True age", breaks = seq( 10, 80, by = 10 ), labels = seq( 10, 80, by = 10 ), limits = c( 10, 80 ) ) +
                scale_y_continuous( "Predicted age", breaks = seq( 10, 80, by = 10 ), labels = seq( 10, 80, by = 10 ), limits = c( 10, 80 ) ) +
                ggtitle( "True vs. Predicted Age" )
ggsave( filename = paste( "brainAge.pdf", sep = "" ), plot = brainAgePlot, width = 6, height = 6, units = 'in' )



# plot the results (Bland Altman)
#
# ageDifference <- ( testingData$AGE - predictedAge )
# meanAgeDifference <- mean( ageDifference )
# stdAgeDifference <- sd( ageDifference )
#
# upperLineIntercept <- meanAgeDifference + 0.95 * stdAgeDifference
# lowerLineIntercept <- meanAgeDifference - 0.95 * stdAgeDifference
#
# plotData <- data.frame( TrueAge = testingData$AGE, AgeDifference = ageDifference )
#
# brainAgePlot <- ggplot( plotData, aes( x = TrueAge, y = AgeDifference ) ) +
#                 geom_point( colour = "darkred", size = 4, alpha = 0.75 ) +
#                 geom_hline( aes( yintercept = upperLineIntercept ), colour = "navyblue", linetype = "dashed" ) +
#                 geom_hline( aes( yintercept = lowerLineIntercept ), colour = "navyblue", linetype = "dashed" ) +
#                 scale_x_continuous( "True age", breaks = seq( 10, 80, by = 10 ), labels = seq( 10, 80, by = 10 ), limits = c( 10, 80 ) ) +
#                 scale_y_continuous( "True - predicted age", breaks = seq( -20, 20, by = 5 ), labels = seq( -20, 20, by = 5 ), limits = c( -20, 20 ) ) +
#                 ggtitle( "True vs. Predicted Age" )
# ggsave( filename = paste( "brainAge.pdf", sep = "" ), plot = brainAgePlot, width = 8, height = 6, units = 'in' )
