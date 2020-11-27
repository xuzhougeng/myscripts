#!/usr/bin/env perl


use List::Util qw[min max];
use POSIX qw(floor);

use warnings;
use strict;

my $v = 3;
my $gtf  = $ARGV[0];

print STDERR "\# " . (localtime) . ": Computing flanking region size for "
    . "AUGUSTUS training genes\n" if ($v > 2);

my $size = 0;
my %gene;

open( GTF, "<", $gtf ) or die("ERROR in file " . __FILE__ ." at line "
    . __LINE__ ."\nCould not open file $gtf!\n");

while (<GTF>) {
    if (m/\tCDS\t/) {
        chomp;
        my @gtfLine = split(/\t/);
        $gtfLine[8] =~ m/;Parent=([^;]+)/;
        if( not( defined( $gene{$1}{'start'} ) ) ) {
            $gene{$1}{'start'} = min( $gtfLine[3], $gtfLine[4] );
        }elsif( $gene{$1}{'start'} > min( $gtfLine[3], $gtfLine[4] ) ) {
            $gene{$1}{'start'} = min( $gtfLine[3], $gtfLine[4] );
        }
        if( not( defined( $gene{$1}{'stop'} ) ) ) {
            $gene{$1}{'stop'} = max($gtfLine[3], $gtfLine[4]);
        }elsif( $gene{$1}{'stop'} < max( $gtfLine[3], $gtfLine[4] ) ) {
            $gene{$1}{'stop'} = max( $gtfLine[3], $gtfLine[4] );
        }
    }
}
close(GTF) or die("ERROR in file " . __FILE__ ." at line ". __LINE__
    ."\nCould not close file $gtf!\n");
my $nGenes   = 0;
my $totalLen = 0;
my $avLen    = 0;
foreach my $key ( keys %gene ) {
    $nGenes++;
    $totalLen += $gene{$key}{'stop'} - $gene{$key}{'start'} +1
}
$avLen = $totalLen / $nGenes;
$size = min( ( floor( $avLen / 2 ), 10000 ) );
if ( $size < 0 ) {
    print STDERR "#*********\n"
            . "# WARNING: \$flanking_DNA has the value $size , which is "
            . "smaller than 0. Something must have gone wrong, there. "
            . "Replacing by value 10000.\n"
            . "#*********\n" if ($v > 0);
    $size = 10000;
}

open( FLANK , ">flank_size.txt" ) or die("Could not open file\n");
print FLANK "$size";
close(FLANK) or die("Could not close file\n");
