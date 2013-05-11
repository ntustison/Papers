resultsX <- list()
resultsX[[1]] <- read.csv( file = "../IXI_results/labelresults.csv", header = TRUE )
resultsX[[2]] <- read.csv( file = "../Kirby_results/labelresults.csv", header = TRUE )
resultsX[[3]] <- read.csv( file = "../NKI_results/labelresults.csv", header = TRUE )
resultsX[[4]] <- read.csv( file = "../Oasis_results/labelresults.csv", header = TRUE )


results <- rbind( resultsX[[1]], resultsX[[2]], resultsX[[3]], resultsX[[4]] )

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
                      "L postcentral",       "R postcentral" )



for( i in seq( 1, 32, by = 2 ) )
  {
  left <- results[,4+i];
  right <- results[,4+i+1];

  asym <- t.test( left, right, alternative = "two.sided", paired = TRUE )

  meandiff <- mean( left - right )

  prefix <- ''
  if( asym$p.value <= 0.01 )
    {
    prefix <- '**'
    } else if( asym$p.value <= 0.1 )
    {
    prefix <- '*'
    }


  cat( prefix, sub( "L ", '', cortical_labels[i] ), ": p-value = ", asym$p.value,
    ", mean difference = ", meandiff, "\n", sep = "" )
  }
