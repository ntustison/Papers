# load required libraries
library( boot )
library( ANTsR )
library( ggplot2 )
library( scales )
library( igraph )
library( pheatmap )
library( randomForest )
library( reshape )
library( e1071 )

sigma <- 5
ages <- seq( 10, 80, by = 1 )

#################################################################
##
##   function definitions
##     * makeGraph - returns a vector of scores or weights for
##                   each vertex in the network based on
##     * calculateCorrelationMatrix -
##
#################################################################

makeGraph <- function( myrsfnetworkcorrs )
  {
  correlationThreshold <- 0.000001

  numberOfNeighbors <- nrow( myrsfnetworkcorrs )
  if( numberOfNeighbors == 0 )
    {
    return( 0 )
    }

  edgeWeights <- c( 0 )
  for( x in c( 1:numberOfNeighbors ) )
    {
    for( y in c( x:numberOfNeighbors ) )
      {
      if( myrsfnetworkcorrs[x,y] > correlationThreshold & myrsfnetworkcorrs[x,y] < 1 )
        {
        edgeWeights <- c( edgeWeights, myrsfnetworkcorrs[x,y] )
        }
      }
    }

  edgeWeights <- edgeWeights[2:length( edgeWeights )]
  adjacencyMatrix <- as.matrix(
    ( myrsfnetworkcorrs > correlationThreshold & myrsfnetworkcorrs < 1 ),
    nrow = numberOfNeighbors, ncol = nnumberOfNeighbors )
  g1 <- graph.adjacency( adjacencyMatrix, mode = c( "undirected" ) )

#   gmetric <- betweenness( g1, normalized = T, weights = 1/edgeWeights )
#   gmetric <- closeness( g1, normalized = T, weights = 1/edgeWeights )
#   gmetric <- page.rank( g1, weights = 1/edgeWeights )$vector
#   gmetric <- degree( g1 )

  gmetric <- transitivity( g1, type = "local", isolates = c( "zero" ) )
  return( gmetric )
  }

calculateCorrelationMatrix <- function( mat, weights, nuis )
  {
  if( missing( mat ) | missing( weights ) )
    {
    print( args( calculateCorrelationMatrix ) )
    return( 1 )
    }
  correlationMatrix <- matrix( rep( NA, ncol( mat ) * ncol( mat ) ), ncol = ncol( mat ) )
  for ( x in 1:ncol( mat ) )
    {
    for ( y in 1:ncol( mat ) )
      {
      correlationMatrix[x,y] <- corr( cbind( mat[,x], mat[,y] ), w = weights / max( weights ) )
      }
    }
  return( correlationMatrix )
  }

#################################################################
##
##   main routine
##
#################################################################

resultsIXI <- read.csv( 'labelresultsI.csv' )
resultsKirby <- read.csv( 'labelresultsK.csv' )
resultsNKI <- read.csv( 'labelresultsN.csv' )
resultsOasis <- read.csv( 'labelresultsO.csv' )

resultsCombined <- rbind( resultsIXI, resultsKirby, resultsNKI, resultsOasis );
resultsCombined$SITE <- as.factor( resultsCombined$SITE )

corticalLabels <- c( "L occipital", "R occipital",
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
                     "L postcentral",       "R postcentral" )


numberOfAges <- length( ages )
numberOfLabels <- length( corticalLabels )
thicknessColumns <- 6:ncol( resultsCombined )

weightedAges <- rep( NA, length( ages ) )
corrs <- rep( NA, numberOfAges )

tstatisticMatrix <- matrix( rep( NA, numberOfLabels * numberOfAges ), ncol = numberOfAges )
pvalueMatrix <- matrix( rep( NA, numberOfLabels * numberOfAges ), ncol = numberOfAges )
networkMales <- matrix( rep( NA, numberOfLabels * numberOfAges ), ncol = numberOfAges )
networkFemales <- matrix( rep( NA, numberOfLabels * numberOfAges ), ncol = numberOfAges )
networkAll <- matrix( rep( NA, numberOfLabels * numberOfAges ), ncol = numberOfAges )
rownames(networkAll)<-corticalLabels
colnames(networkAll)<-ages
count <- 1
for( age in ages )
  {
  # calculate the age differences of the entire cohort
  # with the current age.  Used to calculate weights (i.e.
  # weights closer to the current age are weighted more
  # heavily

  ageDifference <- ( resultsCombined$AGE - age )
  resultsSubsetBasedOnAgeDifference <- subset( resultsCombined, abs( ageDifference ) < 5000*sigma )
  thicknessValues <- resultsSubsetBasedOnAgeDifference[,thicknessColumns]

  ageDifference <- ( resultsSubsetBasedOnAgeDifference$AGE - age )
  cweights <- exp( -1.0 * ageDifference^2 / sigma^2 )
  cweights <- cweights / sum( cweights )

  weightedAges[count] <- sum( resultsSubsetBasedOnAgeDifference$AGE * cweights )
  correlationThicknessMatrix <- calculateCorrelationMatrix( thicknessValues , weights = cweights  )

  ############################################
  #
  # permutation testing on sex differences
  #
  ############################################

  maximumNumberOfPermutations <- 1

  initialNetworkDifference <- c();

  permutationCount <- rep( 0, numberOfLabels )
  for( permutation in 0:maximumNumberOfPermutations )
    {
    myth <- residuals( lm( as.matrix( thicknessValues ) ~ resultsSubsetBasedOnAgeDifference$SITE ) )
    samplesSex <- resultsSubsetBasedOnAgeDifference$SEX
    if( permutation > 0 & permutation < maximumNumberOfPermutations )
      {
      samplesSex <- sample( resultsSubsetBasedOnAgeDifference$SEX )
      }

    temp <- calculateCorrelationMatrix( myth, cweights)
    subnet0 <- reduceNetwork( temp, N = 200 )
    networkAll[,count] <- makeGraph( subnet0$network )
    
    males <- samplesSex == 1
    temp <- calculateCorrelationMatrix( myth[males,], cweights[males] )
    subnet0 <- reduceNetwork( temp, N = 200 )
    networkMales[,count] <- makeGraph( subnet0$network )

    females <- samplesSex == 2
    temp <- calculateCorrelationMatrix( myth[females,], cweights[females] )
    subnet0 <- reduceNetwork( temp, N = 200 )
    networkFemales[,count] <- makeGraph( subnet0$network )

    networkDifference <- abs( networkMales[,count] -  networkFemales[,count] )
    
    if( permutation == 0 | permutation == maximumNumberOfPermutations )
      {
      initialNetworkDifference <- networkDifference
      } else {
      permutationCount <- permutationCount + as.numeric( networkDifference >= initialNetworkDifference )
      }
    }
  cat( "Age: ", age, "\n", sep = '' );
  cat( "  p-value (permutation testing) =", permutationCount / maximumNumberOfPermutations, "\n", sep = ' ' )

  ############################################
  #
  #  Perform t-statistic on age vs. thickness label
  #
  ############################################

  for( corticalLabel in 1:length( corticalLabels ) )
    {
    genderTest <- summary( lm( thicknessValues[,corticalLabel] ~ SEX + SITE + VOLUME, weights = cweights, data = resultsSubsetBasedOnAgeDifference ) )

    # get t-statistic and p-value on SEX significance (respectively)
    tstatisticMatrix[corticalLabel,count] <- coef( genderTest )[2,3]
    pvalueMatrix[corticalLabel,count] <- coef( genderTest )[2,4]
    }

  corrs[count] <- mean( cor( thicknessValues ) )

  subnet0 <- reduceNetwork( correlationThicknessMatrix, N = 200 )
  corrs[count] <- mean( subnet0$network[subnet0$network > 0] )

  meanthicknessValues <- mean( apply( thicknessValues, FUN = mean, MARGIN = 2 ) )

  count <- count+1
  }

##################################
#
#  Create plots
#
##################################

# (weighted) age vs. average thickness network plot

corrsPlotData <- data.frame( weightedAges = weightedAges, correlationValues = corrs )
corrsPlot <- ggplot( corrsPlotData, aes( x = weightedAges, y = correlationValues ) ) +
             geom_point( colour = "darkred", size = 4, alpha = 0.75 ) +
             stat_smooth( colour = "navyblue", formula = y ~ 1 + x + I(x^2) + I(x^3) + I(x^4), method = "lm",
                          size = 1, n = 1000, level = 0.95, se = TRUE, fullrange = TRUE, fill = 'black', alpha = 0.25 ) +
             scale_x_continuous( "Age", breaks = seq( 10, 80, by = 10 ), labels = seq( 10, 80, by = 10 ), limits = c( 10, 80 ) ) +
             scale_y_continuous( "Correlation", breaks = seq( 0.7, 1, by = 0.05 ), labels = seq( 0.7, 1, by = 0.05 ), limits = c( 0.685, 1 ) ) +
             ggtitle( "Average thickness network correlation vs. age" )
ggsave( filename = paste( "averageThicknessNetworkWithAge.pdf", sep = "" ), plot = corrsPlot, width = 8, height = 6, units = 'in' )

# (weighted) age vs. average thickness network plot

qvalueMatrix <- p.adjust( pvalueMatrix, method = "bonferroni" )
qvalueMatrix <- matrix( as.numeric( qvalueMatrix < 0.01 ) , nrow = nrow( tstatisticMatrix ) )
qvalueMatrix[ ( qvalueMatrix == 1 ) & ( tstatisticMatrix < 0 ) ] <- ( -1 )

qvalueData <- data.frame( qvalueMatrix )
colnames( qvalueData ) <- ages
qvalueData$CorticalLabels <- factor( corticalLabels, levels = rev( corticalLabels ) )

qvaluePlot <- ggplot( melt( qvalueData ), aes( x = variable, y = CorticalLabels, fill = value ) ) +
              geom_tile( colour = "darkred" ) +
              scale_fill_gradientn( name = "q-value", colours = heat.colors( 7 ) ) +
              scale_x_discrete( 'Age' ) +
              scale_y_discrete( 'Cortical Labels' )
ggsave( filename = "qvalueHeatMap.pdf", plot = qvaluePlot, width = 10, height = 6, units = 'in' )


tstatisticData <- data.frame( tstatisticMatrix )
colnames( tstatisticData ) <- ages
tstatisticData$CorticalLabels <- factor( corticalLabels, levels = rev( corticalLabels ) )

tstatisticPlot <- ggplot( melt( tstatisticData ) ) +
              geom_tile( aes( x = variable, y = CorticalLabels, fill = value ), colour = "gray50", size = 0 ) +
              scale_fill_gradientn( name = "t-statistic", colours = heat.colors( 7 ) ) +
              scale_x_discrete( 'Age' ) +
              scale_y_discrete( 'Cortical Labels' )
ggsave( filename = "tstatisticHeatMap.pdf", plot = tstatisticPlot, width = 10, height = 6, units = 'in' )


networkMaleData <- data.frame( networkMales )
colnames( networkMaleData ) <- ages
networkMaleData$CorticalLabels <- factor( corticalLabels, levels = rev( corticalLabels ) )

networkMalePlot <- ggplot( melt( networkMaleData ) ) +
               geom_tile( aes( x = variable, y = CorticalLabels, fill = value ), colour = "gray50", size = 0 ) +
               scale_fill_gradientn( name = "correlation\nvalues", colours = heat.colors( 7 ) ) +
               scale_x_discrete( 'Age' ) +
               scale_y_discrete( 'Cortical Labels' ) +
               ggtitle( "Male network" )
ggsave( filename = "maleNetwork.pdf", plot = networkMalePlot, width = 10, height = 6, units = 'in' )


networkFemaleData <- data.frame( networkFemales )
colnames( networkFemaleData ) <- ages
networkFemaleData$CorticalLabels <- factor( corticalLabels, levels = rev( corticalLabels ) )

networkFemalePlot <- ggplot( melt( networkFemaleData ) ) +
               geom_tile( aes( x = variable, y = CorticalLabels, fill = value ), colour = "gray50", size = 0 ) +
               scale_fill_gradientn( name = "correlation\nvalues", colours = heat.colors( 7 ) ) +
               scale_x_discrete( 'Age' ) +
               scale_y_discrete( 'Cortical Labels' ) +
               ggtitle( "Female network" )
ggsave( filename = "femaleNetwork.pdf", plot = networkFemalePlot, width = 10, height = 6, units = 'in' )


pvals<-rep(NA,nrow(networkAll))
for ( n in 1:32 )
  {
  dd<-summary( lm( networkFemales[n,] ~ I(ages) + I(ages)^2 ) )
  dd<-summary( lm( networkMales[n,] ~  I(ages) + I(ages)^2 ) )
  dd<-summary( lm( networkAll[n,] ~  I(ages) + I(ages^2) ) )
  pvals[n]<-coefficients(dd)[3,4]
  }
print( "transitivity with age" ) 
print( p.adjust(pvals,method='BH') )
