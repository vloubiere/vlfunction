#!/usr/bin/perl
use strict;
use warnings;

# Check if two index sequences are given as arguments
if (@ARGV != 2) {
    die "Usage: $0 <index1_prefix> <index2_prefix>\n";
}

# Get the index sequences from command line arguments
my ($index1_prefix, $index2_prefix) = @ARGV;

# Process input from the pipe or standard input
while (<STDIN>) {
    # Read four lines at a time for each FASTQ entry
    my $header = $_;
    my $sequence = <STDIN>;
    my $plus_line = <STDIN>;
    my $quality = <STDIN>;

    # Parse the header to extract the information
    # Example header format:
    # @LH00402:64:22MYG2LT3:8:1101:1074:1064 1:N:0:GATCAGATCTGGGCC+ACGTCCTGGTGTAGA
    if ($header =~ /^@(.+)\s+1:N:0:([A-Za-z]+)\+([A-Za-z]+)/) {
        my ($read_info, $index1, $index2) = ($1, $2, $3);

        # Check if the beginning of the indexes matches the provided prefixes
        if ($index1 =~ /^$index1_prefix/ && $index2 =~ /^$index2_prefix/) {
            # Print the FASTQ entry if the read passes the filter
            print $header;
            print $sequence;
            print $plus_line;
            print $quality;
        }
    }
}

