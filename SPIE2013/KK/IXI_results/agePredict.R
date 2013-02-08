library( ggplot2 )

demographics <- read.csv( file = "demographics.csv", header = TRUE );
# volumes <- read.csv( file = "volume.csv", header = TRUE );
thickness <- read.csv( file = "thickness.csv", header = TRUE );

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

nSamples <- 100;
male_color = "darkred";
female_color = "navyblue";

age <- demographics$AGE;
gender <- cut( demographics$SEX, breaks = c( 0.5, 1.5, 2.5 ), label = c( "male", "female" ) );

sampled.th <- seq( min( thickness ), max( thickness ) , length = nSamples );
sampled.age <- seq( 20, 90, length = nSamples );

for( i in 1:32 )
  {
  th <- thickness[,i];

#   n <- lm( formula = age ~ th + gender );
#   n <- lm( formula = age ~ md + th + I(th^2) + vo + th:vo + gender );
#   n <- lm( formula = age ~ md + th + vo + gender );

    n <- lm( formula = th ~ age + I(age^2) + gender + age:gender );

#   coeffs <- coef( n );
#   pvalues <- anova( n )[,"Pr(>F)"];

  male.data <- list();
  female.data <- list();

  male.data$x <- data.frame( age = sampled.age, th = sampled.th, gender = rep( "male", nSamples ) );
  female.data$x <- data.frame( age = sampled.age, th = sampled.th, gender = rep( "female", nSamples ) );

  male.data$y <- predict( n, male.data$x, interval = "p", level = 0.95 );
  female.data$y <- predict( n, female.data$x, interval = "p", level = 0.95 );

  gender[which( gender == 1 )] <- 'male';
  gender[which( gender == 2 )] <- 'female';

  plotData <- data.frame( cbind( Age = age, Thickness = th, Gender = gender ) )
  plotData <- transform( plotData, Gender = factor( Gender ) );

  thickPlot <- ggplot( plotData, aes( x = Age, y = Thickness, group = Gender ) ) +
               stat_smooth( aes( group = Gender, colour = Gender ), formula = y ~ x + I(x^2), method = "lm", size = 1, n = 1000, level = 0.95, se = TRUE, fullrange = TRUE, fill = 'black', alpha = 0.5 ) +
#                geom_smooth( aes( group = Gender, colour = Gender ), formula = y ~ x + I(x^2), method = "lm", size = 1, n = 1000, level = 0.95, se = TRUE, fill = 'black', alpha = 0.5 ) +
               geom_point( data = plotData, aes( colour = Gender, shape = Gender ), size = 3 ) +
               scale_x_continuous( "Age (years)", breaks = seq( 20, 90, by = 10 ), labels = seq( 20, 90, by = 10 ), limits = c( 20, 90 ) ) +
               scale_y_continuous( "Thickness (mm)", breaks = seq( 0, 5, by = 1 ), labels = seq( 0, 5, by = 1 ), limits = c( 0, 5 ) ) +
               scale_colour_manual( values = c( "navyblue", "darkred" ), breaks = c( 1, 2 ), labels = c( "Male", "Female" ) ) +
               scale_shape_manual( values = c( 18, 16 ), breaks = c( 1, 2 ), labels = c( "Male", "Female" ) ) +
#                scale_shape_discrete( breaks = c( 1, 2 ), labels = c( "Male", "Female" ) ) +
               theme( legend.justification = c( 0, 0 ), legend.position = c( 0, 0 ) ) +
               ggtitle( paste( "Cortical thickness (", cortical_labels[i], ")", sep = "" ) )
  ggsave( filename = paste( "yylabel", i, "_results.pdf", sep = "" ), plot = thickPlot, width = 8, height = 6, units = 'in' )


#   pdf( file = paste( "yylabel", i, "_results.pdf", sep = "" ), 8, 6 );
#
#   plot( c( 20, 90 ), c( 0, 6 ), type = "n",
#     xlab = expression( paste( bold( "Age (years)" ) ) ),
#     ylab = expression( paste( bold( "Thickness (mm)" ) ) ),
#     frame.plot = FALSE  );
#   title( main = paste( "Cortical thickness results (", cortical_labels[i], ")", sep = "" ), col.main = "black", font.main = 2 );
#
#   # plot male
#   points( age[gender == "male"], th[gender == "male"], pch = 22, lwd = 0.5, cex = 0.75, col = male_color );
#   lines( male.data$x[,1], male.data$y[,1], lty=1, col = male_color, lwd = 3 );
#   lines( male.data$x[,1], male.data$y[,2], col = male_color , lty = 2, lwd = 2 )
#   lines( male.data$x[,1], male.data$y[,3], col = male_color , lty = 2, lwd = 2 )
#
#   # plot female
#   points( age[gender == "female"], th[gender == "female"], pch = 23, lwd = 0.5, cex = 0.75, col = female_color );
#   lines( female.data$x[,1], female.data$y[,1], lty=1, col = female_color, lwd = 3 );
#   lines( female.data$x[,1], female.data$y[,2], col = female_color , lty = 2, lwd = 2 )
#   lines( female.data$x[,1], female.data$y[,3], col = female_color , lty = 2, lwd = 2 )
#
#   legend( "topright", c( "Female", "Male" ),
#      col = c( female_color, male_color ),
#    pch = c( 23, 22 ),
#    title = expression( paste( bold( "Gender" ) ) ),
#   # box.lwd = 2,
#     bty = "n",
#     lwd = 2 );
#   axis( 1, tick = TRUE, lwd = 2, font = 2 );
#   axis( 2, tick = TRUE, lwd = 2, font = 2 );
#   dev.off();

  }
