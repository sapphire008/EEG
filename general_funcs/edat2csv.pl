#!/usr/bin/perl -w
# This script reads in a edat (in text format) file, and exports a reduced set of data in a CSV format.
# The source edat file name is passed in using the command line.
# An usage example would be "./ConArrow_to_csv.pl foo". Where "foo" is the name of the edat file to be processed.
# The file produced will have the name "foo.csv".
# 
# The script works by matching edat column header names against a
# list of header names, to be extracted, contained within this script. (@outheaders)
#
# If you wish to modify the columns of data to be extracted just motify the list
# contained witin the braces "( )" on the line marked @outheaders.
#
# Note!!! If you do modify this file please rename it before you save it!!!
# You do not want to muck things up for the person who 
# is expecting this code to work in the "normal" manner.
# Dennis Thompson 8/21/06
@outheaders = qw(PaceProbeBoxed.RT ProbeBoxed.RT PaceProbeBoxed.ACC ProbeBoxed.ACC Condition ProbeTime);

# test for input file on command line 
$len = scalar @ARGV;
if( $len < 1 ) {die "Need to name input file\n"};

#open input file
$infile =  $ARGV[0];
open INPUT, $infile or die "Can't open output file $outfile\n";

#open output file
@outfile = split "\\.", $infile;
$len = scalar @outfile;
#print "$len";

if($len == 2){$outfile = ">$outfile[0].csv";}
else { $outfile = ">./$outfile[$len-2].csv";}

open OUTPUT, $outfile or die "Can't open output file\n";


#drop first 3 lines in file, keep 4th line (is headers) 
for (my $i = 0; $i < 4; $i++){
	$line = <INPUT>;
}

# split line into separate strings using tab delimiter as marker
local $/ = "\r\n";  #in case dos format
chomp $line;
local $/ = "\n"; #set back to default
chomp $line;

@inheaders = split "\t", $line;

# match the header strings to find the column numbers to extract
$n = 0;
foreach my $out (@outheaders){
	$i = 0;
	foreach my $in (@inheaders){
	 	if ( $in eq $out ){
                      #  print "$in  $out $i \n";
			$index[$n] = $i;
			$header{$i} = $out;
			$n++;
		}
		$i++;
	}
}


#do a numeric sort of the index to get headers in the correct order
@newindex = sort { $a <=> $b } @index;

#foreach (@newindex){print"$header{$_} \n";}

#print the headers to the output file
$len = scalar @newindex; # len is the number of headers
$i=0;
foreach (@newindex) {
        $tmp = $header{$_};
        $tmp =~ s/\.//g;
        $tmp =~ s/\]//g;
        $tmp =~ s/\[//g;
        print OUTPUT "$tmp";
	if( $i < $len-1 ){print OUTPUT ",";}
	$i++;
}
print OUTPUT "\n";

# print data to file
while(<INPUT>){
        local $/ = "\r\n";  #in case dos format
        chomp $_;
        local $/ = "\n"; #set back to default
        chomp $_;

	chomp $_;
	@data = split "\t", $_;
	 $k = 0;
	 $i = 0;
	foreach (@data){
		# if index number matches column to be extracted
		if ( $i == $newindex[$k]){ 
			# first check for null data -- replace with NaN
			if ($data[$i] eq 'NULL'){ $data[$i] = 'NaN';}
			# if index == last column print just the data
			if( $k == ($len-1)){print OUTPUT "$data[$i]";}	
			# if index < last column print data and comma 
			if( $k < ($len-1)){ $k++; print OUTPUT "$data[$i],";}	
		}
		$i++;
	}
	print OUTPUT "\n";
}

close OUTPUT;
close INPUT;
