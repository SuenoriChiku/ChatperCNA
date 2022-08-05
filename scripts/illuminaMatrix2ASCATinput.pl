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

&usage() if ( @ARGV < 3 );
{
    my $prefix = "ASCAT";
    my $sampleMapFile = "";
    my $posFile = "";
    GetOptions('prefix=s'=>\$prefix, 'map=s'=>\$sampleMapFile, 'pos=s'=>\$posFile);
    my $csvFile = shift @ARGV;
    illuminaMatrix2ASCATinput( $prefix, $sampleMapFile, $posFile, $csvFile );
}
exit(0);

#
# Subroutines
#
# usage
sub usage {
    print STDERR "Usage:\n";
    print STDERR "./illuminaMatrix2ASCATinput.pl [--prefix LM] [--map GSE44297_SampleMap.txt] --pos GC_Illumina660k.txt GSE44297_Genotype_Matrix_processed.csv.gz\n";
    exit(1);
}

sub illuminaMatrix2ASCATinput {
    my ( $prefix, $sampleMapFile, $posFile, $csvFile ) = @_;

    if ( $posFile eq "" ) {
	usage();
    }
    # $idMap{orgID} = newID
    my %idMap;
    if ( -f $sampleMapFile ) {
	print STDERR "Reading $sampleMapFile...\n";
	readTSVHashKeyRecord( $sampleMapFile, \%idMap, 0, 1, 0 );
    }
    # $posData{rs}[0] = chr
    # $posData{rs}[1] = pos
    # $posData{rs}[2] = $j
    my %posData;
    print STDERR "Reading $posFile...\n";
    readSNPposFile( $posFile, \%posData );
    
    # $lbData[$j][$i][0] = Log R Ratio
    # $lbData[$j][$i][1] = B Allele Freq
    # $lbData[$j][$i][2] = Genotype
    my @lbData;
    # $sampleIDs[$i] = sample ID
    my @sampleIDs;
    print STDERR "Reading $csvFile...\n";
    readGenomeStudioMatrix( $csvFile, \%idMap, \@sampleIDs, \%posData, \@lbData );

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

    my $lrrFile = $prefix . "_Tumor_LogR.txt";
    print STDERR "Outputting $lrrFile...\n";
    outputASCATinput( $lrrFile, \@sampleIDs, \@chrPos, \@lbData, 0 );
    my $bafFile = $prefix . "_Tumor_BAF.txt";
    print STDERR "Outputting $bafFile...\n";
    outputASCATinput( $bafFile, \@sampleIDs, \@chrPos, \@lbData, 1 );
    # Germline
    $lrrFile = $prefix . "_Germline_LogR.txt";
    $bafFile = $prefix . "_Germline_BAF.txt";
    print STDERR "Outputting $lrrFile and $bafFile...\n";
    outputGermlineInput( $lrrFile, $bafFile, \@sampleIDs, \@chrPos, \@lbData );
}

sub readSNPposFile {
    my ( $posFile, $posData ) = @_;

    my $in;
    openFile( $posFile, \$in );
    my $line = readline($in); # Header
    # Body
    while ( my $line = readline($in) ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line =~ /^#/ || $line eq "" ) {
	    next;
	}
	# $terms[0] = SNP Name
	# $terms[1] = Chr
	# $terms[2] = Position
	# $terms[3] = Probe
	# $terms[4] = 200bp
	my @terms = split(/\t/, $line);
	if ( $terms[1] eq "X" ) {
	    $terms[1] = 23;
	} elsif ( $terms[1] eq "Y" ) {
	    $terms[1] = 24;
	} elsif ( $terms[1] eq "XY" || $terms[1] eq "MT"
		  || $terms[1] eq "0" || $terms[2] eq "0" ) {
	    #print STDERR "Skipping chr$terms[1]: $terms[$s]\t$terms[1]:$terms[$p]\t$terms[$l]\t$terms[$b]\n";
	    next;
	}
	if ( exists($$posData{$terms[0]}) ) {
	    print STDERR "Error. $terms[0] are duplicated. Skipping.\n";
	    next;
	}
	@{$$posData{$terms[0]}} = ($terms[1],$terms[2],"NA");
    }
    close( $in );
}

sub readGenomeStudioMatrix {
    my ( $csvFile, $idMap, $sampleIDs, $posData, $lbData ) = @_;

    my $in;
    openFile( $csvFile, \$in );
    # Seach columns
    my $tagLine = readline($in);
    $tagLine =~ s/[\r\n]+\z//;
    my @tags = split(/,/, $tagLine );
    # $jList[$i][0] = Log R Ratio $j
    # $jList[$i][1] = B Allele Freq $j
    # $jList[$i][2] = Genotype $j
    my @jList;
    for( my $j = 1; $j < @tags; ) {
	my $orgID = $tags[$j];
	$orgID =~ s/\..*//;
	my $newID = $orgID;
	if ( exists($$idMap{$orgID}) ) {
	    $newID = $$idMap{$orgID};
	    print STDERR "$orgID is replaced by $newID\n";
	}
	push( @$sampleIDs, $newID );
	my @entry = (-1)x3;
	for( ; $j < @tags; $j++ ) {
	    if ( $tags[$j] !~ /$orgID/ ) {
		last;
	    } elsif ( $tags[$j] =~ /GType/ ) {
		$entry[2] = $j;
	    } elsif ( $tags[$j] =~ /Log R Ratio/ ) {
		$entry[0] = $j;
	    } elsif ( $tags[$j] =~ /B Allele Freq/ ) {
		$entry[1] = $j;
	    }
	}
	if ( $entry[0] < 0 || $entry[1] < 0 || $entry[2] < 0 ) {
	    print STDERR "Error. Input $csvFile format is illegal: $orgID $newID\n";
	    print STDERR "GType: $entry[2], Log R Ratio: $entry[0], B Allele Freq: $entry[1]\n";
	    exit(1);
	}
	push( @jList, \@entry );
    }
    # Body
    while ( my $line = readline($in) ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line =~ /^#/ || $line eq "" ) {
	    next;
	}
	# $terms[0] = Name
	# $terms[1] = Asingh-660W_A01.GType
	# $terms[2] = Asingh-660W_A01.Score
	# $terms[3] = Asingh-660W_A01.Theta
	# $terms[4] = Asingh-660W_A01.R
	# $terms[5] = Asingh-660W_A01.B Allele Freq
	# $terms[6] = Asingh-660W_A01.Log R Ratio
	# $terms[7] = Asingh-660W_B01.GType
	# $terms[8] = ....
	my @terms = split(/,/, $line);
	if ( !exists($$posData{$terms[0]}) ) {
	    next;
	}
	my $j = scalar(@$lbData);
	$$posData{$terms[0]}[2] = $j;
	for( my $i = 0; $i < @jList; $i++ ) {
	    $$lbData[$j][$i][0] = $terms[$jList[$i][0]];
	    $$lbData[$j][$i][1] = $terms[$jList[$i][1]];
	    $$lbData[$j][$i][2] = $terms[$jList[$i][2]];
	}
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
    close( OUT);    
}

sub outputGermlineInput {
    my ( $lrrFile, $bafFile, $sampleIDs, $chrPos, $lbData ) = @_;

    open( LRR, "> $lrrFile") || die "Cannot open $lrrFile\n";
    print LRR "SNPid\tchr\tpos\t";
    outputTSVLine( \*LRR, $sampleIDs );
    
    open( BAF, "> $bafFile") || die "Cannot open $bafFile\n";
    print BAF "SNPid\tchr\tpos\t";
    outputTSVLine( \*BAF, $sampleIDs );
    foreach my $record ( @$chrPos ) {
	my $chr = $$record[1];
	if ( $$record[1] == 23 ) {
	    $chr = "X";
	} elsif ( $$record[1] == 24 ) {
	    $chr = "Y";
	}
	print LRR "$$record[0]\t$chr\t$$record[2]";
	print BAF "$$record[0]\t$chr\t$$record[2]";
	my $j = $$record[3];
	for( my $i = 0; $i < @{$$lbData[$j]}; $i++ ) {
	    if ( $$lbData[$j][$i][2] eq "AA" ) {
		print LRR "\t0.0";
		print BAF "\t0.0";
	    } elsif ( $$lbData[$j][$i][2] eq "AB" ) {
		print LRR "\t0.0";
		print BAF "\t0.5";
	    } elsif ( $$lbData[$j][$i][2] eq "BB" ) {
		print LRR "\t0.0";
		print BAF "\t1.0";
	    } else {
		#print STDERR "\t$$lbData[$j][$i][2]";
		print LRR "\tNaN";
		print BAF "\tNaN";
	    }
	}
	print LRR "\n";
	print BAF "\n";
    }
    close(LRR);
    close(BAF);
}    
