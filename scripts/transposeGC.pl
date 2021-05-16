#!/usr/bin/env perl
#
# Transpose "bedtools nuc" results
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
    transposeGC( \@ARGV );
}
exit(0);

#
# Subroutines
#
# usage
sub usage {
    print STDERR "Usage:\n";
    print STDERR "./transposeGC.pl *.gcContent > GC_OmniChip.txt\n";
    exit(1);
}
         
sub transposeGC {
    my ( $gcFiles ) = @_;

    my @gcEntry = ();
    outGCheader( $$gcFiles[0], \*STDOUT, \@gcEntry );
    foreach my $gc ( @$gcFiles ) {
	print STDERR "Reading $gc...\n";
	readOutTransposeGC( $gc, \*STDOUT, \@gcEntry );
    }
    outputTSVLine( \*STDOUT, \@gcEntry );
}

sub outGCheader {
    my ( $gcFile, $out, $gcEntry ) = @_;

    my $in;
    openFile( $gcFile, \$in );
    my $fline = readline($in); # Header
    $fline = readline($in); # For pre-process
    # $terms[0] = chr
    # $terms[1] = start
    # $terms[2] = end
    # $terms[3] = rs
    # $terms[4] = target postion
    # $terms[5] = length
    # $terms[6] = pct_at
    # $terms[7] = pct_gc
    # $terms[8] = num_A
    # $terms[9] = num_C
    # $terms[10]= num_G
    # $terms[11]= num_T
    # $terms[12]= num_N
    # $terms[13]= num_oth
    # $terms[14]= seq_len
    my @terms = split(/\t/, $fline);
    $terms[0] =~ s/^chr//;
    @$gcEntry = ($terms[3], $terms[0], $terms[4]);
    my @tags = ("Chr","Position", $terms[5]);
    my $limtM = 1000000;
    while ( my $line = readline($in) ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line =~ /^#/ || $line eq "" ) {
	    next;
	}
	my @terms = split(/\t/, $line);
	if ( $$gcEntry[0] ne $terms[3] ) {
	    outputTSVLine( $out, \@tags );
	    $terms[0] =~ s/^chr//;
	    last;
	}
	if ( $terms[5] >= $limtM ) {
	    $terms[5] /= $limtM;
	    $terms[5] .= "M";
	}
	push( @tags, $terms[5] );
    }
    close($in);
}

sub readOutTransposeGC {
    my ( $gc, $out, $gcEntry ) = @_;
    
    my $in;
    openFile( $gc, \$in );
    while ( my $line = readline($in) ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line =~ /^#/ || $line eq "" ) {
	    next;
	}
	my @terms = split(/\t/, $line);
	if ( $$gcEntry[0] ne $terms[3] ) {
	    outputTSVLine( $out, \@$gcEntry );
	    $terms[0] =~ s/^chr//;
	    @$gcEntry = ($terms[3], $terms[0], $terms[4], $terms[7]);
	} else {
	    push( @$gcEntry, $terms[7] );
	}
    }
    close( $in );    
}
