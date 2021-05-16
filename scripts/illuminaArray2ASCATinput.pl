#!/usr/bin/env perl
#
# GenomeStudio output to ASCAT input
#
use strict;
use warnings;
use File::Basename;
use FindBin;
use Getopt::Long;
if ( -f "$FindBin::Bin/utilities.pl" ) {
    require("$FindBin::Bin/utilities.pl");
} else {
    require("utilities.pl");
}

&usage() if ( @ARGV < 1 );
{
    my $prefix = "ASCAT";
    my $sampleMapFile = "";
    GetOptions('prefix=s' => \$prefix, 'map=s'=>\$sampleMapFile);
    illuminaArray2ASCATinput( $prefix, $sampleMapFile, \@ARGV );
}
exit(0);

#
# Subroutines
#
# usage
sub usage {
    print STDERR "Usage:\n";
    print STDERR "./illuminaArray2ASCATinput.pl [--prefix LM] [--map GSE44297_SampleMap.txt] *array.txt.gz\n";
    exit(1);
}

sub illuminaArray2ASCATinput {
    my ( $prefix, $sampleMapFile, $arrayFiles ) = @_;

    # $idMap{orgID} = newID
    my %idMap;
    if ( -f $sampleMapFile ) {
	print STDERR "Reading $sampleMapFile...\n";
	readTSVHashKeyRecord( $sampleMapFile, \%idMap, 0, 1, 0 );
    }
    # $lbData[$j][$i][0] = Log R Ratio
    # $lbData[$j][$i][1] = B Allele Freq
    my @lbData;
    # $posData{SNP Name}[0] = Chr
    # $posData{SNP Name}[1] = Position
    # $posData{SNP Name}[2] = $j
    my %posData;
    # $sampleIDs[$i] = sample ID
    my @sampleIDs;
    foreach my $aFile ( @$arrayFiles ) {
	print STDERR "Reading $aFile...\n";
	readGenomeStudioOutput( $aFile, \%idMap, \@sampleIDs, \%posData, \@lbData );
    }
    # $chrPos[$l][0] = SNP Name
    # $chrPos[$l][1] = Chr
    # $chrPos[$l][2] = Posision
    # $chrPos[$l][3] = $j
    my @chrPos;
    foreach my $rs ( keys(%posData) ) {
	my $ref = $posData{$rs};
	push(@chrPos,[($rs,$$ref[0],$$ref[1],$$ref[2])]);
    }
    print STDERR "Sorting...\n";
    @chrPos = sort{$a->[1] <=> $b->[1] || $a->[2] <=> $b->[2]} @chrPos;

    my $lrrFile = $prefix . "_LogR.txt";
    print STDERR "Outputing $lrrFile...\n";
    outputASCATinput( $lrrFile, \@sampleIDs, \@chrPos, \@lbData, 0 );

    my $bafFile = $prefix . "_BAF.txt";
    print STDERR "Outputing $bafFile...\n";
    outputASCATinput( $bafFile, \@sampleIDs, \@chrPos, \@lbData, 1 );
}

sub readGenomeStudioOutput {
    my ( $aFile, $idMap, $sampleIDs, $posData, $lbData ) = @_;

    my $i = scalar(@$sampleIDs);
    my $orgID = basename($aFile, ".gz");
    $orgID =~ s/_.*//;
    my $newID = $orgID;
    if ( exists($$idMap{$orgID}) ) {
	$newID = $$idMap{$orgID};
	print STDERR "$orgID is replaced by $newID\n";
    }
    push( @$sampleIDs, $newID );
    
    my $in;
    openFile( $aFile, \$in );
    while ( my $line = readline($in) ) { # Skip header
	if ( $line =~ /\[Data\]/ ) {
	    last;
	}
    }
    # Seach columns
    my $tagLine = readline($in);
    $tagLine =~ s/[\r\n]+\z//;
    my @tags = split(/\t/, $tagLine );
    my ($s, $c, $p, $l, $b ) = (-1) x 5;
    for( my $j = 0; $j < @tags; $j++ ) {
	if ( $tags[$j] eq "SNP Name" ) {
	    $s = $j;
	} elsif ( $tags[$j] eq "Chr" ) {
	    $c = $j;
	} elsif ( $tags[$j] eq "Position" ) {
	    $p = $j;
	} elsif ( $tags[$j] eq "Log R Ratio" ) {
	    $l = $j;
	} elsif ( $tags[$j] eq "B Allele Freq" ) {
	    $b = $j;
	}
    }
    if ( $s < 0 || $c < 0 || $p < 0 || $l < 0 || $b < 0 ) {
	print STDERR "Error. Input $aFile format is illegal.\n";
	print STDERR "\tSNP Name:$s, Chr:$c, Position:$p, Log R Ratio: $l, B Allele Freq:$b\n";
	return;
    }
    # Body
    while ( my $line = readline($in) ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line =~ /^#/ || $line eq "" ) {
	    next;
	}
	# $terms[$s] = SNP Name
	# $terms[$c] = Chr
	# $terms[$p] = Position
	# $terms[$l] = Log R Ratio
	# $terms[$b] = B Allele Freq
	my @terms = split(/\t/, $line);
	if ( $terms[$c] eq "X" ) {
	    $terms[$c] = 23;
	} elsif ( $terms[$c] eq "Y" ) {
	    $terms[$c] = 24;
	} elsif ( $terms[$c] eq "XY" || $terms[$c] eq "MT"
		  || $terms[$c] eq "0" || $terms[$p] eq "0" ) {
	    #print STDERR "Skipping chr$terms[$c]: $terms[$s]\t$terms[$c]:$terms[$p]\t$terms[$l]\t$terms[$b]\n";
	    next;
	}
	if ( !exists($$posData{$terms[$s]}) ) {
	    my $j = scalar(@$lbData);
	    @{$$posData{$terms[$s]}} = ($terms[$c],$terms[$p],$j);
	    for( my $k = 0; $k < $i; $k++ ) {
		$$lbData[$j][$k][0] = "NaN";
		$$lbData[$j][$k][1] = "NaN";
	    }
	}
	my $j = $$posData{$terms[$s]}[2];
	$$lbData[$j][$i][0] = $terms[$l];
	$$lbData[$j][$i][1] = $terms[$b];
    }
    close( $in );
}

sub outputASCATinput {
    my ( $outFile, $sampleIDs, $chrPos, $lbData, $p ) = @_;

    open( OUT, "> $outFile") || die "Cannot open $outFile\n";
    print OUT "SNPid\tchr\tpos\t";
    outputTSVLine( \*OUT, $sampleIDs );
    foreach my $record ( @$chrPos ) {
	my $chr = $$record[1];
	if ( $$record[1] == 23 ) {
	    $chr = "X";
	} elsif ( $$record[1] == 24 ) {
	    $chr = "Y";
	}
	print OUT "$$record[0]\t$chr\t$$record[2]";
	my $j = $$record[3];
	for( my $i = 0; $i < @{$$lbData[$j]}; $i++ ) {
	    print OUT "\t$$lbData[$j][$i][$p]";
	}
	print OUT "\n";
    }
    close(OUT);
}
