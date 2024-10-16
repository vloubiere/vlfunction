#!/usr/bin/perl
use strict;
use warnings;

# Check if the correct number of arguments is provided
if ($#ARGV != 2) {
    die "Args should specify 1/ regexpr for 14th column, 2/ regexpr for start of 10th column, and 3/ output file prefix\n";
}

# Read and adjust the first argument
my $regex_14th_col = $ARGV[0];
if ($regex_14th_col !~ /^\^BC:Z:/) {
    print "The provided regular expression for the 14th column does not start with '^BC:Z:'. Adding it automatically.\n";
    $regex_14th_col = '^BC:Z:' . $regex_14th_col;
}

# Prepare the output file name
my $output1 = "$ARGV[2]\_1.fq";

# Delete the output file if it already exists
if (-e $output1) {
    unlink($output1) or die "Cannot delete existing output file: $!";
}

# Open the output file for writing
open(my $out1, ">>", $output1) or die "Cannot open output file: $!";

# Variable to track if any paired-end reads are found
my $paired_reads_found = 0;

# Process the input from stdin
while (<STDIN>) {
    my @cur = split(/\t/, $_);

    # Check if the read is the first in a pair and unmapped (FLAG = 77)
    if ($cur[1] == 77 && $cur[13] =~ /$regex_14th_col/ && $cur[9] =~ /^$ARGV[1]/) {
        # Mark that we have found at least one paired read
        $paired_reads_found = 1;

        # Trim the second regular expression from the start of the read
        $cur[9] =~ s/^$ARGV[1]//;

        # Extract the first 10 nucleotides from the trimmed read
        my $first_10_nucleotides = substr($cur[9], 0, 10);

        # Trim the first 10 nucleotides from the read
        $cur[9] = substr($cur[9], 10);

        # Trim the corresponding quality scores
        $cur[10] = substr($cur[10], 14);

        # Append the first 10 nucleotides to the read ID
        my $new_read_id = $cur[0] . ":" . $first_10_nucleotides;

        # Write the FASTQ formatted output
        print $out1 "@" . $new_read_id . "\n";
        print $out1 $cur[9] . "\n";
        print $out1 "+\n";
        print $out1 $cur[10] . "\n";
    }
}

# Close the output file
close($out1);

# Check if no paired reads were found and print a warning
if (!$paired_reads_found) {
    die "No paired-end reads found in the BAM file. Please check if the BAM file is single-ended.";
}
