#!/usr/bin/env perl
#
# Change cnvMatrix from segment-file
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
    my $segmentFile = shift @ARGV;
    makeCNVmatrix( $segmentFile );
}
exit(0);

#
# Subroutines
#
# usage
sub usage {
    print "Usage:\n";
    print "./makeCNVmatrix.pl ASCAT.segment > ASCAT.matrix\n";
    exit(1);
}

sub makeCNVmatrix {
    my ( $segmentFile ) = @_;

    # $cnvData{ID}{chr}[$j][0] = ID-D
    # $cnvData{ID}{chr}[$j][1] = chr
    # $cnvData{ID}{chr}[$j][2] = start
    # $cnvData{ID}{chr}[$j][3] = end
    # $cnvData{ID}{chr}[$j][4] = copy number
    my %cnvData;
    print STDERR "Reading $segmentFile...\n";
    readSegmentFile( $segmentFile, \%cnvData );
    print STDERR "Making cnvMatrix...\n";
    outputCnvMatrix( \*STDOUT, \%cnvData );
}

sub readSegmentFile {
    my ( $segmentFile, $cnvData ) = @_;

    my $in;
    openFile( $segmentFile, \$in );
    my $line = readline($in); # header
    my $a = 4;
    my $b = 5;
    if( scalar(split(/\t/, $line)) > 6 ) {
	$a = 7;
	$b = 8;
    }
    while ( my $line = readline($in) ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line eq "" || $line =~ /^#/ ) {
	     next;
	}
	# $terms[0] = ID               # $terms[0] = sampleID
	# $terms[1] = chr	       # $terms[1] = chr     
	# $terms[2] = start (1-offset) # $terms[2] = startpos
	# $terms[3] = end (1-offset)   # $terms[3] = endpos  
	# $terms[4] = startSNP	       # $terms[4] = nMajor  
	# $terms[5] = endSNP	       # $terms[5] = nMinor  
	# $terms[6] = number of SNPs
	# $terms[7] = nA
	# $terms[8] = nB
	# $terms[9] = nA+nB
	# $terms[10]= over 2-copy LOH flag (1:LOH, 0:not LOH)
	my @terms = split( /\t/, $line );
	if ( $terms[1] =~ /M/ ) {
	    next;
	} elsif ( $terms[1] eq "X" ) {
	    $terms[1] = 23;
	} elsif ( $terms[1] eq "Y" ) {
	    $terms[1] = 24;
	}
	if ( !exists($$cnvData{$terms[0]}) ) {
	    %{$$cnvData{$terms[0]}} = ();
	}
	if ( !exists($$cnvData{$terms[0]}{$terms[1]}) ) {
	    @{$$cnvData{$terms[0]}{$terms[1]}} = ();
	}
	$terms[6] = $terms[$a] + $terms[$b];
	if ( $terms[$b] == 0 && $terms[6] > 1
	    && $terms[1] ne "23" && $terms[1] ne "24" ) {
	    $terms[6] = 1.5;
	}
	push( @{$$cnvData{$terms[0]}{$terms[1]}}, \@terms );
    }
    close($in);
}

# $cnvMatrix[$j][0] = chr
# $cnvMatrix[$j][1] = start
# $cnvMatrix[$j][2] = end
# $cnvMatrix[$j][3] = ID-1
# $cnvMatrix[$j][$i] = ID-i
sub outputCnvMatrix {
    my ( $out, $cnvData ) = @_;

    # $breaks{c}{pos}
    my %breaks;
    foreach my $id ( keys %$cnvData ) {
	foreach my $c ( keys %{$$cnvData{$id}} ) {
	    #print STDERR "$id\t$c\n";
	    my $ref = $$cnvData{$id}{$c};
	    for ( my $j = 0; $j < @$ref; $j++ ) {
		if ( !exists($breaks{$c}{$$ref[$j][2]}) ) {
		    $breaks{$c}{$$ref[$j][2]} = 1;
		}
		if ( !exists($breaks{$c}{$$ref[$j][3]}) ) {
		    $breaks{$c}{$$ref[$j][3]} = 1;
		}
	    }
	}
    }
    # header
    my @ids = sort(keys(%$cnvData));
    print $out "chr\tstart\tend\t";
    outputTSVLine( $out, \@ids );
    # output CNV
    my $nColumn = scalar(@ids)+3;
    foreach my $c ( sort{ $a <=> $b } keys %breaks ) {
	# $pList[$l] = pos
	my @pList = sort{ $a <=> $b } keys %{$breaks{$c}};
	# $pHash{pos} = $l
	my %pHash;
	for( my $l = 0; $l < @pList; $l++ ) {
	    $pHash{$pList[$l]} = $l;
	}
	# $regionIDs{$id}{s}{e} = copy number
	my %regionIDs;
	foreach my $id ( @ids ) {
	    if ( !exists($$cnvData{$id}{$c}) ) {
		next;
	    }
	    my $ref = $$cnvData{$id}{$c};
	    for( my $j = 0; $j < @$ref; $j++ ) {
		for(my $l=$pHash{$$ref[$j][2]};$l < $pHash{$$ref[$j][3]};$l++){
		    $regionIDs{$id}{$pList[$l]}{$pList[$l+1]} = $$ref[$j][6];
		}
	    }
	}
	# output
	my $chr = $c;
	if ( $c eq "23" ) {
	    $chr = "X";
	} elsif ( $c eq "24" ) {
	    $chr = "Y";
	}
	my @preLine = (0) x $nColumn;
	my @oneLine = ();
	my $preFlag = 0;
	for( my $l = 1; $l < @pList; $l++ ) {
	    my $s = $pList[$l-1];
	    my $e = $pList[$l];
	    @oneLine = ($chr, $s, $e );
	    my $flag = 0;
	    foreach my $id ( @ids ) {
		if ( exists($regionIDs{$id}{$s}{$e}) ) {
		    push( @oneLine, $regionIDs{$id}{$s}{$e} );
		    $flag = 1;
		} else {
		    push( @oneLine, "NA" );
		}
	    }
	    if ( compPreLine( \@preLine, \@oneLine ) ) {
		if ( $preFlag ) {
		    outputTSVLine( $out, \@preLine );
		}
		@preLine = @oneLine;
	    }
	    $preFlag = $flag;
	}
	# last line
	if ( $preLine[2] == $oneLine[2] ) {
	    outputTSVLine( $out, \@preLine );
	} else {
	    outputTSVLine( $out, \@oneLine );
	}
    }
}

sub compPreLine {
    my ( $preLine, $oneLine ) = @_;

    for( my $j = 3; $j < @$preLine; $j++ ) {
	if ( $$preLine[$j] ne $$oneLine[$j] ) {
	    return 1;
	}
    }
    $$preLine[2] = $$oneLine[2];
    return 0;
}
