#!/usr/bin/perl
use strict;
use warnings;

if ($#ARGV != 2) {
    print "Args should specify 1/ regexpr for 14th column, 2/ regexpr for start of 10th column, and 3/ output file prefix\n";
    exit;
}

my $regex_14th_col = $ARGV[0];
if ($regex_14th_col !~ /^\^BC:Z:/) {
    print "The provided regular expression for the 14th column does not start with '^BC:Z:'. Adding it automatically.\n";
    $regex_14th_col = '^BC:Z:' . $regex_14th_col;
}

my $output1 = "$ARGV[2]\_1.fq";

if (-e $output1) {
    unlink($output1);
}

open(my $out1, ">>", $output1) or die "Cannot open output file: $!";

while (<STDIN>) {
    my @cur = split(/\t/, $_);
    # 4 sam flag = single-end read (this should be checked depending on your SAM file specifics)
    if ($cur[1] == 4 && $cur[13] =~ /$regex_14th_col/ && $cur[9] =~ /^$ARGV[1]/) { 
        # Trim the second regular expression from the start of the read
        $cur[9] =~ s/^$ARGV[1]//;

        # Extract the first 10 nucleotides and their quality scores from the trimmed read
        my $first_10_nucleotides = substr($cur[9], 0, 10);

        # Trim the first 10 nucleotides and their quality scores from the read
        $cur[9] = substr($cur[9], 10);
        
        # Trim the quality scores similar to the read (4nt eBC and 10nt UMI)
        $cur[10] = substr($cur[10], 14);

        # Append the first 10 nucleotides to the read ID
        my $new_read_id = $cur[0] . ":" . $first_10_nucleotides;

        say $out1 "@" . $new_read_id;
        say $out1 $cur[9];
        say $out1 "+";
        say $out1 $cur[10];
    }
}
close($out1);