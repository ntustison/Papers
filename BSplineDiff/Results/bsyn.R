library( ggplot2 )

stats <- read.csv( "overlapMeasures3.csv", header = TRUE )

# difference <- stats$bsyn - stats$syn
# average <- 0.5 * ( stats$syn + stats$bsyn )
#
# dmean <- data.frame( x = -c( -Inf, Inf ), y = mean( difference ) )
# pStd <- data.frame( x = -c( -Inf, Inf ), y = mean( difference ) + 1.96 * sd( difference ) )
# mStd <- data.frame( x = -c( -Inf, Inf ), y = mean( difference ) - 1.96 * sd( difference ) )
#
# plotData <- data.frame( cbind( SyN = stats$syn, BSyN = stats$bsyn, Difference = difference, Average = average  ) )
#
# myPlot <- ggplot( plotData, aes( x = Average, y = Difference ) ) +
#           geom_line( aes( x, y ), linetype = 1, dmean ) +
#           geom_line( aes( x, y ), linetype = 2, pStd ) +
#           geom_line( aes( x, y ), linetype = 2, mStd ) +
#           geom_point( aes( colour = Difference ), size = 3 ) +
#           geom_point( colour = "grey90", size = 1 ) +
#           scale_x_continuous( "Average Dice", breaks = seq( 0.6, 0.8, by = 0.05 ), labels = seq( 0.6, 0.8, by = 0.05 ), limits = c( 0.6, 0.8 ) ) +
#           scale_y_continuous( "BSyN Dice - SyN Dice", breaks = seq( -0.025, 0.05, by = 0.025 ), labels = seq( -0.025, 0.05, by = 0.025 ), limits = c( -0.025, 0.05 ) ) +
# #           theme( legend.justification = c( 0, 0 ), legend.position = c( 0, 0 ) ) +
#           scale_shape( solid = FALSE ) +
#           ggtitle( paste(  "SyN vs. BSyN", sep = "" ) )
# ggsave( filename = paste( "syn.pdf", sep = "" ), plot = myPlot, width = 8, height = 5, units = 'in' )


plotLabels <- c( rep.int( "SyN", length( stats$syn ) ), rep.int( "BSplineSyN", length( stats$bsyn ) ) )

plotData <- data.frame( Type = plotLabels, Dice = c( stats$syn, stats$bsyn ) )

myPlot <- ggplot( plotData,  aes( x = factor( Type ), y = Dice ) ) +
          geom_boxplot( aes( fill = Type ) ) +
          ggtitle( "MICCAI 2012 Multi-Label Atlas Challenge Data" ) +
          scale_x_discrete( "" )

ggsave( filename = paste( "syn.pdf", sep = "" ), plot = myPlot, width = 8, height = 5, units = 'in' )

