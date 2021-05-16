#!/usr/bin/env perl
#
# Add nProbes, rs# and centromere breaks
#
use strict;
use warnings;
use FindBin;
use Getopt::Long;
if ( -f "$FindBin::Bin/utilities.pl" ) {
    require("$FindBin::Bin/utilities.pl");
} else {
    require("utilities.pl");
}

&usage() if ( @ARGV < 2 );
{
    my $cenFile = "";
    GetOptions('centromere=s'=>\$cenFile);
    my $posFile = shift @ARGV;
    addProbeInfoSegment( $posFile, $cenFile, \@ARGV );
}
exit(0);

#
# Subroutines
#
# usage
sub usage {
    print "Usage:\n";
    print "./addProbeInfoSegment.pl [-c centromere.txts] omniPos.txt ASCAT.segments... > ASCAT.segments.txt\n";
    exit(1);
}
         
sub addProbeInfoSegment {
    my ( $posFile, $cenFile, $segmentFiles ) = @_;

    # $posHash{chr}{pos}[0] = rs
    # $posHash{chr}{pos}[1] = p or q
    # $posHash{chr}{pos}[2] = n
    my %posHash;
    print STDERR "Reading $posFile...\n";
    readPosFile( $posFile, \%posHash );

    if ( -f $cenFile ) {
	print STDERR "Adding $cenFile...\n";
	addPorQ( $cenFile, \%posHash );
    }

    # $segment[$j][0] = ID
    # $segment[$j][1] = chr
    # $segment[$j][2] = start (1-offset)
    # $segment[$j][3] = end (1-offset)
    # $segment[$j][4] = startSNP
    # $segment[$j][5] = endSNP
    # $segment[$j][6] = number of SNPs
    # $segment[$j][7] = nA
    # $segment[$j][8] = nB
    # $segment[$j][9] = nA+nB
    # $segment[$j][10]= over 2-copy LOH flag (1:LOH, 0:not LOH)
    print "ID\tchr\tstart\tend\tstartSNP\tendSNP\tnumSNPs\tnMajor\tnMinor\tcopyNumber\to2LOH\n";
    foreach my $file ( @$segmentFiles ) {
	print STDERR "Reading $file...\n";
	outputSegment( \*STDOUT, $file, \%posHash );
    }
}

sub readPosFile {
    my ( $posFile, $posHash ) = @_;

    my $in;
    openFile( $posFile, \$in );
    my $n = 0;
    while ( my $line = readline($in) ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line eq "" || $line =~ /^#/ || $line =~ /^Chr\tPosition/ ) {
	    next;
	}
	# $terms[0] = rs
	# $terms[1] = chr
	# $terms[2] = pos
	my @terms = split( /\t/, $line );
	my @entry;
	$terms[0] =~ s/^chr//;
	$entry[0] = $terms[0];
	$entry[1] = "NA";
	$entry[2] = $n++;
	$$posHash{$terms[1]}{$terms[2]} = \@entry;
    }
    close($in);
}

sub addPorQ {
    my ( $cenFile, $posHash ) = @_;

    # $centromere{chr} = pos
    my %centromere;
    my $in;
    openFile( $cenFile, \$in );
    while ( my $line = readline($in) ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line eq "" || $line =~ /^#/ ) {
	    next;
	}
	# $terms[0] = chr
	# $terms[1] = pos
	my @terms = split( /\t/, $line );
	$terms[0] =~ s/^chr//;
	if ( exists($centromere{$terms[0]}) ) {
	    print STDERR "Warning. chr$terms[0] are duplicated. Skipping.";
	    next;
	}
	$centromere{$terms[0]} = $terms[1];
    }
    close($in);
    # Add p or q
    foreach my $c ( keys %$posHash ) {
	if ( !exists($centromere{$c}) ) {
	    print STDERR "Error. chr$c is not found in $cenFile.\n";
	    exit( 0 );
	}
	my $flag = 1;
	my $prePos = 0;
	my $n = 0;
	foreach my $pos ( sort{$a <=> $b} keys %{$$posHash{$c}} ) {
	    $n++;
	    if ( $pos < $centromere{$c} ) {
		$$posHash{$c}{$pos}[1] = "p";
		$$posHash{$c}{$pos}[2] = $n;
		$prePos = $pos;
	    } else {
		$$posHash{$c}{$pos}[1] = "q";
		$$posHash{$c}{$pos}[2] = $n;
		if ( $flag ) {
		    $$posHash{$c}{"p"} = $prePos;
		    $$posHash{$c}{"q"} = $pos;
		    $flag = 0;
		}
	    }
	}
    }
}

sub outputSegment {
    my ( $out, $cnvFile, $posHash ) = @_;

    my $in;
    openFile( $cnvFile, \$in );
    my $line = readline($in); # header
    my $prevC = 0;
    while ( my $line = readline($in) ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line eq "" || $line =~ /^#/ ) {
	    print STDERR "Skipping $line\n";
	    next;
	}
	# $terms[0] = sampleID
	# $terms[1] = chr
	# $terms[2] = startpos
	# $terms[3] = endpos
	# $terms[4] = nMajor
	# $terms[5] = nMinor
	my @terms = split( /\t/, $line );
	if ( !exists($$posHash{$terms[1]}{$terms[2]})
	     || !exists($$posHash{$terms[1]}{$terms[3]}) ) {
	    print STDERR "Error. chr$terms[1]:$terms[2]-$terms[3] is not found in SNP-pos-file.\n";
	    exit( 0 );
	}
	# $entry[0] = ID
	# $entry[1] = chr
	# $entry[2] = start (1-offset)
	# $entry[3] = end (1-offset)
	# $entry[4] = startSNP
	# $entry[5] = endSNP
	# $entry[6] = number of SNPs
	# $entry[7] = nA
	# $entry[8] = nB
	# $entry[9] = nA+nB
	# $entry[10]= over 2-copy LOH flag (1:LOH, 0:not LOH)
	my @entry = ($terms[0],$terms[1],$terms[2],$terms[3],0,0,0,$terms[4],$terms[5],$terms[4]+$terms[5],0);
	my $posChr = $$posHash{$terms[1]};
	$entry[4] = $$posChr{$terms[2]}[0];
	$entry[5] = $$posChr{$terms[3]}[0];
	$entry[6] = $$posChr{$terms[3]}[2] - $$posChr{$terms[2]}[2] + 1;
	if ( $entry[8] == 0 && $entry[9] > 1 ) {
	    $entry[10] = 1;
	}
	# Centromere check
	if ( $$posChr{$terms[2]}[1] ne $$posChr{$terms[3]}[1] ) {
	    my @entry2 = @entry;
	    my $p = $$posChr{"p"};
	    my $q = $$posChr{"q"};
	    $entry[3] = $p;
	    $entry[5] = $$posChr{$p}[0];
	    $entry[6] = $$posChr{$p}[2] - $$posChr{$entry[2]}[2] + 1;
	    outputTSVLine( $out, \@entry );
	    $entry2[2] = $q;
	    $entry2[4] = $$posChr{$q}[0];
	    $entry2[6] = $$posChr{$entry2[3]}[2] - $$posChr{$q}[2] + 1;
	    outputTSVLine( $out, \@entry2 );
	} else {
	    outputTSVLine( $out, \@entry );
	}
    }
}
