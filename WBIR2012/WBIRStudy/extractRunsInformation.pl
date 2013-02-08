#! /usr/bin/perl/ -w

my $baseDir = '/Users/ntustison/Data/Public/LONI/LPBA40/WBIRStudy/';
my $dataDir = $baseDir . $ARGV[0];

my @elapsedTimes = ();
my @metricValues = ();

my @metricFileNames = ( "${dataDir}/metricValuesLevel0.csv", "${dataDir}/metricValuesLevel1.csv",
  "${dataDir}/metricValuesLevel2.csv", "${dataDir}/metricValuesLevel3.csv" );

for( my $i = 0; $i < @metricFileNames; $i++ )
  {
  open( FILE, ">$metricFileNames[$i]" );
  close( FILE );
  }

for( my $i = 1; $i <= 10; $i++ )
  {
  for( my $j = 1; $j <= 10; $j++ )
    {
    if( $i == $j )
      {
      next;
      }

    my $i_idx = sprintf( "%02d", $i );
    my $j_idx = sprintf( "%02d", $j );

    my $filename = "S${j_idx}toS${i_idx}_log.txt";

    if( ! -e "${dataDir}/${filename}" )
      {
      next;
      }

    print "Doing ${filename}\n";

    open( FILE, "<${dataDir}/${filename}" );

    my @filecontents = <FILE>;

    close( FILE );
    # first get start of stage line numbers

    my @stageLineNumbers = ();

    for( my $k = 200; $k < @filecontents; $k++ )
      {
      if( $filecontents[$k] =~ m/^Stage/ )
        {
        push( @stageLineNumbers, $k );
        }
      }

    print "@stageLineNumbers\n";
    my @iterations = ();

    for( my $k = $stageLineNumbers[1]; $k < @filecontents; $k++ )
      {
      if( $filecontents[$k] =~ m/Elapsed time \(stage/ )
        {
        my @tmp = split( ': ', $filecontents[$k] );
        chomp( $tmp[1] );
        push( @elapsedTimes, scalar( $tmp[1] ) );
        }
      if( $filecontents[$k] =~ m/iterations/ && $filecontents[$k] =~ m/x/ )
        {
        my @tmp = split( ' = ', $filecontents[$k] );
        print "$filecontents[$k]\n";
        chomp( $tmp[1] );
        @iterations = split( 'x', $tmp[1] );
        }
      }

    my @startLines = ();
    for( my $k = 0; $k < @iterations; $k++ )
      {
      for( my $m = $stageLineNumbers[1]; $m < @filecontents; $m++ )
        {
        if( $filecontents[$m] =~ m/Current level = ${k}/ )
          {
          push( @startLines, $m );
          }
        }
      }
    my $end = @filecontents - 1;
    push( @startLines, $end );

    for( my $k = 0; $k < @startLines - 1; $k++ )
      {
      my @metricValuesPerLevel = ();

      for( my $m = $startLines[$k]; $m < $startLines[$k+1]; $m++ )
        {
        if( $filecontents[$m] =~ m/Iteration/ )
          {
          my @tmp = split( ' ', $filecontents[$m] );
          chomp( $tmp[5] );
          my $metric = $tmp[5];
          $metric =~ s/,//;

          push( @metricValuesPerLevel, $metric );
          }
        }

      for( my $m = @metricValuesPerLevel; $m < $iterations[$k]; $m++ )
        {
        push( @metricValuesPerLevel, 0 );
        }
      open( FILE, ">>$metricFileNames[$k]" );
      my $metricString = join( ',', @metricValuesPerLevel );
      print FILE "${metricString}\n";
      close( FILE );

      push( @metricValues, @metricValuesPerLevel );
      }
    }
  }


open( FILETIMES, ">${dataDir}/elapsedTimes.csv" );
for( $i = 0; $i < @elapsedTimes; $i++ )
  {
  print FILETIMES "$elapsedTimes[$i]\n";
  }
close( FILETIMES );


# `cp ${baseDir}/../native_labels_renamed/S01.nii.gz ${dataDir}/averageTargetOverlapOnS01.nii.gz`;
#
# open( FILE, ">${dataDir}/labelOverlapMeasures.csv" );
#
# for( $n = 0; $n < @averageTargetOverlap; $n++ )
#   {
#   $averageTargetOverlap[$n] /= $count;
#   print FILE "$labels[$n],$averageTargetOverlap[$n]\n";
#   `UnaryOperateImage 3 ${dataDir}/averageTargetOverlapOnS01.nii.gz r 0 ${dataDir}/averageTargetOverlapOnS01.nii.gz $labels[$n] $averageTargetOverlap[$n]`;
#   }
# close( FILE );
#
