#! /usr/bin/perl/ -w

my $baseDir = '/Users/ntustison/Data/Public/LONI/LPBA40/WBIRStudy/';
my $dataDir = $baseDir . $ARGV[0];


my @labels = ();
my @averageTargetOverlap = ();
my $count = 0;

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
    print "Doing S${j_idx}toS${i_idx}\n";

#    my @args = ( 'antsApplyTransforms',
#                 '-d', 3,
#                 '-i', "${baseDir}/native_labels_renamed/S${j_idx}.nii.gz",
#                 '-r', "${baseDir}/native_labels_renamed/S${i_idx}.nii.gz",
#                 '-o', "${baseDir}/WBIRStudy/S${j_idx}toS${i_idx}_WarpedLabels.nii.gz",
#                 '-n', 'NearestNeighbor',
#                 '-t', "${baseDir}/WBIRStudy/S${j_idx}toS${i_idx}_1Warp.nii.gz",
#                 '-t', "${baseDir}/WBIRStudy/S${j_idx}toS${i_idx}_0Affine.mat",
#                 '-v', 0
#               );
#    system( @args ) == 0 || die "Fail.\n";

    my @out = `LabelOverlapMeasures 3 ${dataDir}S${j_idx}toS${i_idx}_WarpedLabels.nii.gz ${baseDir}/../native_labels_renamed/S${i_idx}.nii.gz`;

    for( my $n = 5; $n < @out; $n++ )
      {
      my @labelStats = split( ' ', $out[$n] );;
      if( $i == 1 && $j == 2 )
        {
        $averageTargetOverlap[$n-5] = $labelStats[1];
        $labels[$n-5] = $labelStats[0];
        }
      else
        {
        $averageTargetOverlap[$n-5] += $labelStats[1];
        if( $labels[$n-5] != $labelStats[0] )
          {
          print "Mismatched labels.\n";
          exit;
          }
        }
      }
    $count++;
    }
  }


`cp ${baseDir}/../native_labels_renamed/S01.nii.gz ${dataDir}/averageTargetOverlapOnS01.nii.gz`;

open( FILE, ">${dataDir}/labelOverlapMeasures.csv" );

for( $n = 0; $n < @averageTargetOverlap; $n++ )
  {
  $averageTargetOverlap[$n] /= $count;
  print FILE "$labels[$n],$averageTargetOverlap[$n]\n";
  `UnaryOperateImage 3 ${dataDir}/averageTargetOverlapOnS01.nii.gz r 0 ${dataDir}/averageTargetOverlapOnS01.nii.gz $labels[$n] $averageTargetOverlap[$n]`;
  }
close( FILE );

