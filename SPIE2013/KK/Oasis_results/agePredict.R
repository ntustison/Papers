library( ggplot2 )

results <- read.csv( file = paste( "labelresults", ".csv", sep = '' ), header = TRUE );

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

male_color = "darkred";
female_color = "navyblue";

age <- results$AGE;
gender <- cut( results$SEX, breaks = c( 0.5, 1.5, 2.5 ), label = c( "male", "female" ) );

for( i in 1:32 )
  {
  th <- results[,i+4];

  gender[which( gender == 1 )] <- 'male';
  gender[which( gender == 2 )] <- 'female';

  plotData <- data.frame( cbind( Age = age, Thickness = th, Gender = gender ) )
  plotData <- transform( plotData, Gender = factor( Gender ) );

  thickPlot <- ggplot( plotData, aes( x = Age, y = Thickness, group = Gender ) ) +
#                stat_smooth( aes( group = Gender, colour = Gender ), formula = y ~ 1 + x, method = "lm", size = 1, n = 1000, level = 0.95, se = TRUE, fullrange = FALSE, fill = 'black', alpha = 0.5 ) +
               stat_smooth( aes( group = Gender, colour = Gender ), formula = y ~ 1 + x + I(x^2), method = "lm", size = 1, n = 1000, level = 0.95, se = TRUE, fullrange = TRUE, fill = 'black', alpha = 0.5 ) +
#                geom_smooth( aes( group = Gender, colour = Gender ), formula = y ~ 1 + x + I(x^2), method = "lm", size = 1, n = 1000, level = 0.95, se = TRUE, fill = 'black', alpha = 0.5 ) +
               geom_point( data = plotData, aes( colour = Gender, shape = Gender ), size = 3 ) +
               scale_x_continuous( "Age (years)", breaks = seq( 20, 90, by = 10 ), labels = seq( 20, 90, by = 10 ), limits = c( 20, 90 ) ) +
               scale_y_continuous( "Thickness (mm)", breaks = seq( 0, 6, by = 1 ), labels = seq( 0, 6, by = 1 ), limits = c( 0, 6 ) ) +
               scale_colour_manual( values = c( "navyblue", "darkred" ), breaks = c( 1, 2 ), labels = c( "Male", "Female" ) ) +
               scale_shape_manual( values = c( 18, 16 ), breaks = c( 1, 2 ), labels = c( "Male", "Female" ) ) +
#                scale_shape_discrete( breaks = c( 1, 2 ), labels = c( "Male", "Female" ) ) +
               theme( legend.justification = c( 0, 0 ), legend.position = c( 0, 0 ) ) +
               ggtitle( paste( "Cortical thickness (", cortical_labels[i], ")", sep = "" ) )
  ggsave( filename = paste( "label", i, "_results.pdf", sep = "" ), plot = thickPlot, width = 8, height = 6, units = 'in' )
  }
