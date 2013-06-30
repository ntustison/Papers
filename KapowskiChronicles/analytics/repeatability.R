kirby <- read.csv( "labelResultsK_pairwise.csv" )
oasis <- read.csv( "labelResultsO_pairwise.csv" )

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

#
# The repeat scan is listed at every other row, i.e.
#
#    subject 1, scan 1, ...
#    subject 1, scan 2, ...
#    subject 2, scan 1, ...
#    subject 2, scan 2, ...
#
# We reorganize to list all the first scans in the first
# half of the data frame followed by all the repeat scans
# and remove the 'ID' column
#

subjectID <- paste0( "S", c( 1:( 0.5 * ( nrow( kirby ) + nrow( oasis ) ) ) ) )
subjectID <- rbind( subjectID, subjectID )
subjectID <- as.vector( tmp )

scanOrder <- rep( c( "First", "Repeat" ), 0.5 * ( nrow( kirby ) + nrow( oasis ) ) )

repeatabilityDataFrame <- data.frame( rbind( kirby, oasis ) )
repeatabilityDataFrame$ID <- as.factor( subjectID )
repeatabilityDataFrame$ScanOrder <- as.factor( scanOrder )

# thickness.FirstScan <- c()
# thickness.RepeatScan <- c()
# for( i in 6:37 )
#   {
#   thickness.FirstScan <- c( thickness.FirstScan, repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "First" ), i] )
#   thickness.RepeatScan <- c( thickness.RepeatScan, repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "Repeat" ), i] )
#   }
# myTest <- t.test( thickness.FirstScan,
#                   thickness.RepeatScan,
#                   paired = TRUE )

qvalues <- c()
meanDifferences <- c()
for( i in 6:37 )
  {
  myTest <- t.test( repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "First"),i],
                    repeatabilityDataFrame[which( repeatabilityDataFrame$ScanOrder == "Repeat"),i],
                    paired = TRUE, conf.int = TRUE )
  qvalues[i-5] <- myTest$p.value
  meanDifferences[i-5] <- myTest$estimate
  }
qvalues <- p.adjust( qvalues, method = "fdr" )

for( i in 1:32 )
  {
  cat( corticalLabels[i], ": mean of the differences = ", meanDifferences[i],
       " (q-value = ", qvalues[i], ")\n", sep = '' )
  }

