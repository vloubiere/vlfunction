#!/usr/bin/perl
use strict;
use warnings;

# Check if the correct number of arguments are provided
if ($#ARGV != 1) {
    die "Usage: $0 <regular expression> <output file prefix>\n";
}

# Get the output file prefix
my $output1 = "$ARGV[1]_1.fq";

# Remove existing output file if it exists
if (-e $output1) {
    unlink($output1) or die "Cannot delete existing output file: $!";
}

# Open output file for writing
open(my $out1, ">>", $output1) or die "Cannot open output file $output1: $!";

# Process the SAM file from standard input
while (<STDIN>) {
    chomp;
    my @fields = split(/\t/, $_);

    # Skip header lines
    next if ($fields[0] =~ /^\@/);

    my $flag = $fields[1];

    # Only process single-end reads (not part of a pair)
    if (!($flag & 0x1)) {  # Flag 0x1 indicates the read is paired
        # Check if the read matches the regular expression in the 14th column
        if ($fields[13] =~ /$ARGV[0]/) {
            # Write the read to the output file in FASTQ format
            say $out1 "@" . $fields[0];
            say $out1 $fields[9];
            say $out1 "+";
            say $out1 $fields[10];
        }
    }
}

# Close the output file
close($out1);

print "Processing complete.\n";
