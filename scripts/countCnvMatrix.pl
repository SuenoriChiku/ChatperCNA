#!/usr/bin/env perl
#
# Count copy number for cnvMatrix
#
use strict;
use warnings;
use FindBin;
if ( -f "$FindBin::Bin/utilities.pl" ) {
    require("$FindBin::Bin/utilities.pl");
} else {
    require("utilities.pl");
}

&usage() if ( @ARGV < 1 );
{
    my $cnvFile = shift @ARGV;
    countCnvMatrix( $cnvFile );
}
exit(0);

#
# Subroutines
#
# usage
sub usage {
    print "Usage:\n";
    print "./countCnvMatrix.pl ASCAT.matrix >! ASCAT.count\n";
    exit(1);
}
         
sub countCnvMatrix {
    my ( $cnvFile ) = @_;

    my $in;
    openFile( $cnvFile, \$in );
    my $line = readline($in); # header
    $line =~ s/[\r\n]+\z//;
    print STDERR "Number of sample: ",scalar(split(/\t/,$line))-3,"\n";
    print "chr\tstart\tend\tm2\tm1\tp1\tp2\tp3\to2LOH\tloss\tgain\n";
    while ( my $line = readline($in) ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line eq "" || $line =~ /^#/ ) {
	     next;
	}
	# $terms[0] = chr
	# $terms[1] = start
	# $terms[2] = end
	# $terms[3] = ID-1
	# $terms[4] = ID-2
	# $terms[5] = ...
	my @terms = split( /\t/, $line );
	my @bin = ( 0, 0, 0, 0, 0, 0 );
	my $nEff = 0;
	for( my $i = 3; $i < @terms; $i++ ) {
	    $nEff++;
	    if ( $terms[$i] eq "NA" ) {
		$nEff--;
		next;
	    } elsif ( $terms[$i] == 2 ) {
		next;
	    } elsif ( $terms[$i] == 0 ) {
		$bin[0]++;
	    } elsif ( $terms[$i] == 1 ) {
		$bin[1]++;
	    } elsif ( $terms[$i] == 3 ) {
		$bin[2]++;
	    } elsif ( $terms[$i] == 4 ) {
		$bin[3]++;
	    } elsif ( $terms[$i] >= 5 ) {
		$bin[4]++;
	    } elsif ( $terms[$i] == 1.5 ) {
		$bin[5]++;
	    } else {
		print STDERR "Error. chr$terms[0]:$terms[1]-$terms[2] $i-th count is $terms[$i]\n";
	    }
	}
	if ( $nEff == 0 ) {
	    print STDERR "Warning. chr$terms[0]:$terms[1]-$terms[2] is all NA.\n";
	    next;
	}
	print "$terms[0]\t$terms[1]\t$terms[2]";
	foreach my $count ( @bin ) {
	    printf( "\t%f", $count/$nEff );
	}
	printf( "\t%f\t%f\n", ($bin[0]+$bin[1]+$bin[5])/$nEff, ($bin[2]+$bin[3]+$bin[4])/$nEff );
    }
    close($in);
}
