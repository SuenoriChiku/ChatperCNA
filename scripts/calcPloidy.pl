#!/usr/bin/env perl
#
# calculation for ploidy according to segment length
#
use strict;
use warnings;
use FindBin;
if ( -f "$FindBin::Bin/utilities.pl" ) {
    require("$FindBin::Bin/utilities.pl");
} else {
    require("utilities.pl");
}

&usage() if ( @ARGV < 2 );
{
    my $aberrFile = shift @ARGV;
    my $segmeFile = shift @ARGV;
    calcPloidy( $aberrFile, $segmeFile );
}
exit(0);

#
# Subroutines
#
# usage
sub usage {
    print STDERR "Usage:\n";
    print STDERR "./calcPloidy.pl ASCAT.aberr ASCAT.segment > ASCAT.abpl\n";
    exit(1);
}
         
sub calcPloidy {
    my ( $aberrFile, $segmeFile ) = @_;

    # $aberr[$i][0] = ID
    # $aberr[$i][1] = aberrantcellfraction
    # $aberr[$i][2] = ploidy
    # $aberr[$i][3] = psi
    my @aberr;
    print STDERR "Reading $aberrFile...\n";
    readTSVWithSkip( $aberrFile, 1, \@aberr );

    # $segHash{ID}[$j][0] = ID
    # $segHash{ID}[$j][1] = chr
    # $segHash{ID}[$j][2] = start (1-offset)
    # $segHash{ID}[$j][3] = end (1-offset)
    # $segHash{ID}[$j][4] = startSNP
    # $segHash{ID}[$j][5] = endSNP
    # $segHash{ID}[$j][6] = number of SNPs
    # $segHash{ID}[$j][7] = nA
    # $segHash{ID}[$j][8] = nB
    # $segHash{ID}[$j][9] = nA+nB
    # $segHash{ID}[$j][10]= over 2-copy LOH flag (1:LOH, 0:not LOH)
    #
    # $segHash{ID}[$j][11]= p or q
    # $segHash{ID}[$j][12]= distance
    my %segHash;
    print STDERR "Reading $segmeFile...\n";
    readTSVHashArray( $segmeFile, \%segHash, 0, 1 );

    # $aberr[$i][0] = ID
    # $aberr[$i][1] = aberrantcellfraction
    # $aberr[$i][2] = ploidy by mean of weighted length
    # $aberr[$i][3] = ploidy by number of SNPs
    # $aberr[$i][4] = psi from ASPCF
    insertLengthPloidy( \%segHash, \@aberr );

    outputPloidyHeader( \*STDOUT, $aberrFile );
    @aberr = sort{$a->[0] cmp $b->[0]} @aberr;
    outputTSV( \*STDOUT, \@aberr );
}

sub insertLengthPloidy {
    my ( $segHash, $aberr ) = @_;

    foreach my $record ( @$aberr ) {
#	my $id = $$record[0];
	if ( $$record[0] =~ /[0-9]{4}$/ ) {
	    $$record[0] .= "D";
	}
	my $id = $$record[0];
	for ( my $j = scalar(@$record)+1; $j > 2; $j-- ) {
	    $$record[$j] = $$record[$j-2];
	}
	if ( exists($$segHash{$id}) ) {
	    getLengthPloidy( $$segHash{$id}, \$$record[2], \$$record[3] );
#	} elsif ( exists($$segHash{$idD}) ) {
#	    getLengthPloidy( $$segHash{$idD}, \$$record[2], \$$record[3] );
	} else {
	    $id =~ s/\./\-/g;
	    if ( exists($$segHash{$id}) ) {
		$$record[0] = $id;
		getLengthPloidy( $$segHash{$id}, \$$record[2], \$$record[3] );
	    } else {
		$$record[2] = "NA";
		$$record[3] = "NA";
	    }
	}
    }
}

sub getLengthPloidy {
    my ( $segHashID, $meanL, $meanN ) = @_;

    $$meanL = "NA";
    $$meanN = "NA";
    my $nLen = 0;
    my $len = 0;
    my $n = 0;
    my $nCount = 0;
    foreach my $record ( @$segHashID ) {
	if ( $$record[1] !~ /[0-9]+/ ) {
	    next;
	}
	my $l = $$record[3] - $$record[2] + 1;
	$nLen += $l * $$record[9];
	$len += $l;
	#
	$n += $$record[9];
	$nCount++;
    }
    if ( $len > 0 ) {
	$$meanL = sprintf( "%1.3f", $nLen/$len );
    }
    if ( $nCount> 0 ) {
	$$meanN = sprintf( "%1.3f", $n/$nCount );
    }
}

sub outputPloidyHeader {
    my ( $out, $aberrFile ) = @_;

    print $out "ID\taberr\tploidyL\tploidyN";
    my $line = getFirstLine( $aberrFile );
    chomp( $line );
    my @terms = split( /\t/, $line );
    if ( scalar(@terms) > 2 ) {
	print $out "\tASCATploidy\tASCATpsi";
    }
#    for ( my $j = 1; $j < @terms; $j++ ) {
#	print $out "\tASCAT$terms[$j]";
#    }
    print $out "\n";
}
