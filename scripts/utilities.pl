use strict;
use warnings;
#
# Utilities for perl
#
sub openFile {
    my ( $inputFile, $in ) = @_;

    if ( ! -f $inputFile ) {
	print STDERR "Error. $inputFile is not found.\n";
	exit(1);
    }
    if ( $inputFile =~ /\.gz$/ ) {
	open( $$in, "gzip -dc $inputFile |" ) || die "Cannot open $inputFile\n";
    } elsif ( $inputFile =~ /\.bz2$/ ) {
	open( $$in, "bzip2 -dc $inputFile |") || die "Cannot open $inputFile\n";
    } elsif ( $inputFile =~ /\.zip$/ ) {
	open( $$in, "unzip -p $inputFile |" ) || die "Cannot open $inputFile\n";
    } else {
	open( $$in, "< $inputFile" ) || die "Cannot open $inputFile\n";
    }
}

sub readTSV {
    my ( $inputFile, $entries ) = @_;

    my $in;
    openFile( $inputFile, \$in );
    while ( my $line = readline($in) ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line ne "" ) {
	    push( @$entries, [split( /\t/, $line )] );
	}
    }
    close( $in );
}

sub readTSVWithSkip {
    my ( $inputFile, $numSkip, $entries ) = @_;

    open( IN, "< $inputFile" ) || die "Cannot open $inputFile file\n";
    while ( $numSkip > 0 ) {
	#print "$numSkip\n";
	my $line = <IN>;
	$numSkip--;
    }
    while ( my $line = <IN> ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line ne "" ) {
	    push( @$entries, [split( /\t/, $line )] );
	}
    }
    close( IN );
}

sub readTSVWithComment {
    my ( $inputFile, $entries, $commentSymbol ) = @_;

    open( IN, "< $inputFile" ) || die "Cannot open $inputFile file\n";
    while ( my $line = <IN> ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line ne "" && $line !~ /^$commentSymbol/ ) {
	    push( @$entries, [split( /\t/, $line )] );
	}
    }
    close( IN );
}

sub readTSVPosOnly {
    my ( $inputFile, $entries, $p, $numSkip ) = @_;

    open( IN, "< $inputFile" ) || die "Cannot open $inputFile file\n";
    while ( $numSkip > 0 ) {
	my $line = <IN>;
	$numSkip--;
    }
    while ( my $line = <IN> ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line eq "" || $line =~ /^#/ ) {
	     next;
	}
	my @terms = split(/\t/, $line);
	$terms[$p] =~ s/^\s+//;
	$terms[$p] =~ s/\s+$//;
	push( @$entries, $terms[$p] );
    }
    close( IN );
}

sub readList {
    my ( $inputFile, $entries ) = @_;

    open( IN, "< $inputFile" ) || die "Cannot open $inputFile file\n";
    while ( my $line = <IN> ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line eq "" || $line =~ /^#/ ) {
	     next;
	}
	push( @$entries, $line );
    }
    close( IN );
}

sub readHeader {
    my ( $inputFile, $entries ) = @_;

    open( IN, "< $inputFile" ) || die "Cannot open $inputFile file\n";
    my $line = <IN>;
    $line =~ s/[\r\n]+\z//;
    $line =~ s/^#//;
    @$entries = split( /\t/, $line );
    close( IN );
}

sub readTSVHash {
    my ( $inputFile, $entries, $hashPosition ) = @_;

    open( IN, "< $inputFile" ) || die "Cannot open $inputFile file\n";
    while ( my $line = <IN> ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line ne "" && $line !~ /^#/ ) {
	    my @terms = split( /\t/, $line );
	    if ( !exists($$entries{ $terms[$hashPosition] }) ) {
		$$entries{ $terms[$hashPosition] } = \@terms;
	    } else {
		print STDERR "Warning. $terms[$hashPosition] already exists.\n";
	    }
	}
    }
    close( IN );
}

sub readTSVHashWithoutCheck {
    my ( $inputFile, $entries, $hashPosition ) = @_;

    open( IN, "< $inputFile" ) || die "Cannot open $inputFile file\n";
    while ( my $line = <IN> ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line ne "" && $line !~ /^#/ ) {
	    my @terms = split( /\t/, $line );
	    $$entries{ $terms[$hashPosition] } = \@terms;
	}
    }
    close( IN );
}

sub readTSVHashWithSkip {
    my ( $inputFile, $entries, $hashPosition, $numSkip ) = @_;

    open( IN, "< $inputFile" ) || die "Cannot open $inputFile file\n";
    while ( $numSkip > 0 ) {
	#print "$numSkip\n";
	my $line = <IN>;
	$numSkip--;
    }
    while ( my $line = <IN> ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line ne "" && $line !~ /^#/ ) {
	    my @terms = split( /\t/, $line );
	    if ( !exists($$entries{ $terms[$hashPosition] }) ) {
		$$entries{ $terms[$hashPosition] } = \@terms;
	    } else {
		print STDERR "Warning. $terms[$hashPosition] already exists.\n";
	    }
	}
    }
    close( IN );
}

sub readTSVHashPosOnly {
    my ( $inputFile, $entries, $pos ) = @_;

    open( IN, "< $inputFile" ) || die "Cannot open $inputFile file\n";
    while ( my $line = <IN> ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line ne "" && $line !~ /^#/ ) {
	    my @terms = split( /\t/, $line );
	    if ( !exists($$entries{ $terms[$pos] }) ) {
		$$entries{ $terms[$pos] } = 0;
	    } else {
		print STDERR "Error. $terms[$pos] already exists.\n";
	    }
	}
    }
    close( IN );
}

sub readTSVHashKeyRecord {
    my ( $inputFile, $entries, $p, $q, $numSkip ) = @_;

    open( IN, "< $inputFile" ) || die "Cannot open $inputFile file\n";
    while ( $numSkip > 0 ) {
	my $line = <IN>;
	$numSkip--;
    }
    while ( my $line = <IN> ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line eq "" || $line =~ /^#/ ) {
	     next;
	}
	my @terms = split(/\t/, $line);
	$terms[$p] =~ s/^\s+//;
	$terms[$p] =~ s/\s+$//;
	$terms[$q] =~ s/^\s+//;
	$terms[$q] =~ s/\s+$//;
	if ( exists($$entries{$terms[$p]}) ) {
	    print STDERR "Warning. $terms[$p] already exists. $$entries{$terms[$p]} vs $terms[$q]\n";
	}
	$$entries{$terms[$p]} = $terms[$q];
    }
    close( IN );
}

sub readTSVHashArray {
    my ( $inputFile, $entries, $hashPosition, $numSkip ) = @_;

    open( IN, "< $inputFile" ) || die "Cannot open $inputFile file\n";
    while ( $numSkip > 0 ) {
	#print "$numSkip\n";
	my $line = <IN>;
	$numSkip--;
    }
    while ( my $line = <IN> ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line eq "" || $line =~ /^#/ ) {
	    next;
	}
	my @terms = split( /\t/, $line );
	push( @{$$entries{$terms[$hashPosition]}}, \@terms );
    }
    close( IN );
}

sub readSpaceSeparated {
    my ( $inputFile, $entries ) = @_;

    open( IN, "< $inputFile" ) || die "Cannot open $inputFile file\n";
    while ( my $line = <IN> ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line eq "" || $line =~ /^#/ ) {
	     next;
	}
	my @terms = split( /\s+/, $line );
	while ( $terms[0] eq "" ) {
	    shift( @terms );
	}
	while ( $terms[-1] eq "" ) {
	    pop( @terms );
	}
	push( @$entries, [@terms] );
    }
    close( IN );
}

sub readSpaceSeparatedHash {
    my ( $inputFile, $entries, $hashPosition ) = @_;

    open( IN, "< $inputFile" ) || die "Cannot open $inputFile file\n";
    while ( my $line = <IN> ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line eq "" || $line =~ /^#/ ) {
	    next;
	}
	$line =~ s/^\s+//;
	$line =~ s/\s+$//;
	my @terms = split( /\s+/, $line );
	$$entries{ $terms[$hashPosition] } = \@terms;
    }
    close( IN );
}

sub readSpaceSeparatedWithSkip {
    my ( $inputFile, $numSkip, $entries ) = @_;

    open( IN, "< $inputFile" ) || die "Cannot open $inputFile file\n";
    while ( $numSkip > 0 ) {
	#print "$numSkip\n";
	my $line = <IN>;
	$numSkip--;
    }
    while ( my $line = <IN> ) {
	if ( $line ne "" ) {
	    $line =~ s/[\r\n]+\z//;
	    my @terms = split( /\s+/, $line );
	    #print "length ", length($line), ", size ", scalar(@terms), "\n";
	    while ( $terms[0] eq "" ) {
		shift( @terms );
	    }
	    while ( $terms[-1] eq "" ) {
		pop( @terms );
	    }
	    push( @$entries, [@terms] );
	}
    }
    close( IN );
}

sub readSpaceSeparatedWithComment {
    my ( $inputFile, $entries, $commentSymbol ) = @_;

    open( IN, "< $inputFile" ) || die "Cannot open $inputFile file\n";
    while ( my $line = <IN> ) {
	if ( $line ne "" && $line !~ /^$commentSymbol/ ) {
	    $line =~ s/[\r\n]+\z//;
	    my @terms = split( /\s+/, $line );
	    #print "length ", length($line), ", size ", scalar(@terms), "\n";
	    while ( $terms[0] eq "" ) {
		shift( @terms );
	    }
	    while ( $terms[-1] eq "" ) {
		pop( @terms );
	    }
	    push( @$entries, [@terms] );
	}
    }
    close( IN );
}

sub readCSVWithSkip {
    my ( $inputFile, $numSkip, $entries ) = @_;

    open( IN, "< $inputFile" ) || die "Cannot open $inputFile file\n";
    while ( $numSkip > 0 ) {
	#print "$numSkip\n";
	my $line = <IN>;
	$numSkip--;
    }
    while ( my $line = <IN> ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line ne "" ) {
	    push( @$entries, [split( /,/, $line )] );
	}
    }
    close( IN );
}

sub readCSVHashWithSkip {
    my ( $inputFile, $entries, $hashPosition, $numSkip ) = @_;

    open( IN, "< $inputFile" ) || die "Cannot open $inputFile file\n";
    while ( $numSkip > 0 ) {
	#print "$numSkip\n";
	my $line = <IN>;
	$numSkip--;
    }
    while ( my $line = <IN> ) {
	$line =~ s/[\r\n]+\z//;
	if ( $line ne "" && $line !~ /^#/ ) {
	    my @terms = split( /,/, $line );
	    if ( !exists($$entries{ $terms[$hashPosition] }) ) {
		$$entries{ $terms[$hashPosition] } = \@terms;
	    } else {
		print STDERR "Warning. $terms[$hashPosition] already exists.\n";
	    }
	}
    }
    close( IN );
}

sub splitcsv {
    my ( $csvstr ) = @_;

    $csvstr .= ',';
    $csvstr =~ s/("([^"]|"")*"|[^,]*),/$1$;/g;
    $csvstr =~ s/"([^$;]*)"$;/$1$;/g;
    $csvstr =~ s/""/"/g;

    return split(/$;/, $csvstr);
}

sub pushLineData {
    my ( $entries, $line ) = @_;

    $line =~ s/[\r\n]+\z//;
    my @terms = split( /\s+/, $line );
    while ( $terms[0] eq "" ) {
	shift( @terms );
    }
    while ( $terms[-1] eq "" ) {
	pop( @terms );
    }
    push( @$entries, [@terms] );
}

sub pushLine {
    my ( $entries, $line ) = @_;

    $line =~ s/[\r\n]+\z//;
    my @terms = split( /\s+/, $line );
    while ( $terms[0] eq "" ) {
	shift( @terms );
    }
    while ( $terms[-1] eq "" ) {
	pop( @terms );
    }
    foreach my $entry ( @terms ) {
	push( @$entries, $entry );
    }
}

sub outputValues {
    my ( $out, $outData ) = @_;

    foreach my $v ( @$outData ) {
	print $out "$v\n";
    }
}

# ex.
# outputTSV( \*STDOUT, \@extractEntries );
sub outputTSV {
    my ( $out, $outputEntries ) = @_;

    foreach my $entry ( @$outputEntries ) {
	outputTSVLine( $out, $entry );
    }
}

sub outputTopTSV {
    my ( $out, $num, $outputEntries ) = @_;

    for( my $j = 0; $j < $num; $j++ ) {
	outputTSVLine( $out, $$outputEntries[$j] );
    }
}

sub outputTSVLine {
    my ( $out, $entry ) = @_;

    print $out "$$entry[0]";
    for ( my $i = 1; $i < @$entry; $i++ ) {
#	if ( !defined($$entry[$i]) ) {
#	    print STDERR "Undefined: $$entry[0]\t$i\n";
#	    next;
#	}
	print $out "\t$$entry[$i]";
    }
    print $out "\n";
}

sub outputTSVLineWithoutCR {
    my ( $out, $entry ) = @_;

    for ( my $i = 0; $i < @$entry; $i++ ) {
	print $out "\t$$entry[$i]";
    }
}

sub outputTSVLineOnly {
    my ( $out, $entry ) = @_;

    print $out "$$entry[0]";
    for ( my $i = 1; $i < @$entry; $i++ ) {
	print $out "\t$$entry[$i]";
    }
}

sub outputSSV {
    my ( $out, $outputEntries ) = @_;

    foreach my $entry ( @$outputEntries ) {
	outputSSVLine( $out, $entry );
    }
}

sub outputSSVLine {
    my ( $out, $entry ) = @_;

    print $out "$$entry[0]";
    for ( my $i = 1; $i < @$entry; $i++ ) {
	print $out " $$entry[$i]";
    }
    print $out "\n";
}

sub outputSSVLineWithoutCR {
    my ( $out, $entry ) = @_;

    print $out "$$entry[0]";
    for ( my $i = 1; $i < @$entry; $i++ ) {
	print $out " $$entry[$i]";
    }
}

sub outputCSV {
    my ( $out, $outputEntries ) = @_;

    foreach my $entry ( @$outputEntries ) {
	outputCSVLine( $out, $entry );
    }
}

sub outputCSVLine {
    my ( $out, $entry ) = @_;

    print $out "$$entry[0]";
    for ( my $i = 1; $i < @$entry; $i++ ) {
	print $out ",$$entry[$i]";
    }
    print $out "\n";
}

sub outputHashArray {
    my ( $out, $hashEntries ) = @_;

    foreach my $key ( sort keys %$hashEntries ) {
	outputTSVLine( $out, $$hashEntries{$key} );
    }
}

sub outputHashKeys {
    my ( $out, $hashEntries ) = @_;

    foreach my $key ( sort keys %$hashEntries ) {
	print $out "$key\n";
    }
}

sub outputLatexTable {
    my ( $out, $outputEntries ) = @_;

    foreach my $entry ( @$outputEntries ) {
	outputLatexTableLine( $out, $entry );
    }
}

sub outputLatexTableLine {
    my ( $out, $entry ) = @_;
    
    print $out "$$entry[0]";
    for ( my $i = 1; $i < @$entry; $i++ ) {
	print $out " & $$entry[$i]";
    }
    print $out "\n";
}

sub outputTSVNA {
    my ( $out, $num ) = @_;

    for ( my $i = 0; $i < $num; $i++ ) {
	print $out "\tNA";
    }
}

sub outputCSVNA {
    my ( $out, $num ) = @_;

    for ( my $i = 0; $i < $num; $i++ ) {
	print $out ",NA";
    }
}

sub outputHeader {
    my ( $out, $inputFile ) = @_;

    my $in;
    openFile( $inputFile, \$in );
    while ( my $line = readline($in) ) {
	if ( $line =~ /^#/ ) {
	    print $out "$line";
	} else {
	    last;
	}
    }
    close( $in );
}

sub getHeaders {
    my ( $inputFile, $entries ) = @_;

    my $in;
    openFile( $inputFile, \$in );
    while ( my $line = readline($in) ) {
	if ( $line =~ /^#/ ) {
	     $line =~ s/[\r\n]+\z//;
	     push( @$entries, $line );
	} else {
	    last;
	}
    }
    close( $in );
}

sub outputFirstLine {
    my ( $out, $inputFile ) = @_;

    my $in;
    openFile( $inputFile, \$in );
    my $line = readline($in);
    print $out "$line";
    close( $in );
}

sub getFirstLine {
    my ( $inputFile ) = @_;

    my $in;
    openFile( $inputFile, \$in );
    my $line = readline($in);
    #chomp( $line );
    $line =~ s/[\r\n]+\z//;
    close( $in );
    return $line;
}

#
# ExtractNames
#
# Usage:
# &extractNames( $inputName, \$cutExtension, \$basename, \$header );
#
# ex.
# $inputName = hoge/huge/hege.sjis.txt
#
# $cutExtension = hoge/huge/hege.sjis
# $baseName    = hege.sjis.txt
# $header      = hege.sjis
#
sub extractNames {
    my ( $inputName, $cutExtension, $baseName, $header ) = @_;

    $$cutExtension = &cutExtension( $inputName );
    $$baseName = &basename2( $inputName );
    $$header = &cutExtension( $$baseName );
}

# Delete the last extension
sub cutExtension {
    my ( $inputName ) = @_;
    my $cutExtension = $inputName;
    $cutExtension =~ s/^(.+)\.(.+?)$/$1/; # Delete the last extension
    return $cutExtension;
}

# Delete head directory names
sub basename2 {
    my ( $inputName ) = @_;
    my $baseName = $inputName;
    $baseName =~ s/^(.+)\/(.+?)$/$2/; # Delete the head directories
    return $baseName;
}

# Extract only header
sub extractHeader {
    my ( $inputName ) = @_;
    my $baseName = &basename2( $inputName );
    my $header = &cutExtension( $baseName );
    return $header;
}

sub getExtension {
    my ( $inputName ) = @_;
    $inputName =~ s/^(.+)\.(.+?)$/$2/;
    return $inputName;
}

sub dirName {
    my ( $inputName ) = @_;
    $inputName =~ s/^(.+)\/(.+?)$/$1/;
    return $inputName;
}

sub getDigitValue {
    my ( $x, $digit ) = @_;
    my $l = $digit+3;
    #printf( "%${l}.${digit}f\n", $x );
    return sprintf( "%${l}.${digit}f", $x );
}

sub checkMvFile {
    my ( $fileName ) = @_;

    if ( -f $fileName ) {
	my $mvFile = extractHeader($fileName) . "2." . getExtension($fileName);
	#print STDERR "Warning. $fileName already exists. The original file is moved to $mvFile.\n";
	rename( $fileName, $mvFile );
    }
}

sub getEffDigitValue {
    my ( $x, $digit ) = @_;

    if ( $x !~ /^[\.0-9e-]+$/ ) {
	return( $x );
    }
    return sprintf( "%.${digit}f", $x );
#    my $l = $digit+3;
#    my $str = sprintf( "%${l}.${digit}f", $x );
#    my $pL = 0;
#    if ( $str =~ /\./ ) {
#	$str =~ s/\.[0-9]+//;
#	$pL = $digit - length($str) + 1;
#    }
#    return sprintf( "%${l}.${pL}f", $x );
}

sub getPaddingValue {
    my ( $x, $digit ) = @_;
    return sprintf( "%0${digit}d", $x );
}

sub get3CommaNumber {
    my ( $x ) = @_;

    my $num = $x;
    if ( $num =~ /^[-+]?\d\d\d\d+/g ) {
	my $i = 0;
	my $j = 0;
	for ( $i = pos($num) - 3, $j = $num =~ /^[-+]/; $i > $j; $i -= 3){
	    substr($num, $i, 0) = ',';
	}
    }
    return $num;
}

sub getLatexScienceFormat {
    my ( $x, $digit ) = @_;

    my $num = sprintf( "%.${digit}e", $x );
    $num =~ s/e/ \\times 10\^\{/;
    $num =~ s/([-+])0+/$1/;
    $num = "\$$num\}\$";

    return $num;
}

sub getLatexScienceFormatF {
    my ( $x, $digit ) = @_;

    my $num = "NA";
    my $l = $digit+3;
    if ( $x < 0.01 ) {
	$num = sprintf( "%.${digit}e", $x );
	$num =~ s/e/ \\times 10\^\{/;
	$num =~ s/([-+])0+/$1/;
	$num = "\$$num\}\$";
    } elsif ( $x < 0.1 ) {
	$l++;
	$digit++;
	$num = sprintf( "%${l}.${digit}f", $x );
    } else {
	$num = sprintf( "%${l}.${digit}f", $x );
    }
    return $num;
}

sub ceil {
    my ( $x ) = @_;
    my $a = 0;
    if ( $x > 0 && $x != int ($x) ) {
	$a = 1;
    }
    return int($x + $a);
}

sub max {
    my ( $a, $b ) = @_;
    return ($a > $b) ? $a : $b;
}

sub min {
    my ( $a, $b ) = @_;
    return ($a < $b) ? $a : $b;
}

sub rmHeadTailSpaces {
    my ( $record ) = @_;
    $$record =~ s/^\s+//;
    $$record =~ s/\s+$//;
}

1;
