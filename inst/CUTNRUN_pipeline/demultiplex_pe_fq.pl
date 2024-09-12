#!/usr/bin/perl
use strict;
use warnings;

# Check if two input files and two index sequences are given
if (@ARGV != 4) {
    die "Usage: $0 <R1_fastq> <R2_fastq> <index1_prefix> <index2_prefix>\n";
}

# Get the input files and index sequences from command line arguments
my ($r1_file, $r2_file, $index1_prefix, $index2_prefix) = @ARGV;

# Open the R1 and R2 files
open my $r1_fh, '<', $r1_file or die "Could not open R1 file '$r1_file': $!";
open my $r2_fh, '<', $r2_file or die "Could not open R2 file '$r2_file': $!";

# Process the two files simultaneously
while (my $r1_header = <$r1_fh>) {
    # Read four lines for R1
    my $r1_sequence = <$r1_fh>;
    my $r1_plus_line = <$r1_fh>;
    my $r1_quality = <$r1_fh>;

    # Read four lines for R2
    my $r2_header = <$r2_fh>;
    my $r2_sequence = <$r2_fh>;
    my $r2_plus_line = <$r2_fh>;
    my $r2_quality = <$r2_fh>;

    # Parse the R1 header to extract the index sequences
    if ($r1_header =~ /^@(.+)\s+1:N:0:([A-Za-z]+)\+([A-Za-z]+)/) {
        my ($read_info, $index1, $index2) = ($1, $2, $3);

        # Check if the beginning of the indexes match the provided prefixes
        if ($index1 =~ /^$index1_prefix/ && $index2 =~ /^$index2_prefix/) {
            # Print the R1 and R2 FASTQ entries if the read passes the filter
            print $r1_header;
            print $r1_sequence;
            print $r1_plus_line;
            print $r1_quality;

            print $r2_header;
            print $r2_sequence;
            print $r2_plus_line;
            print $r2_quality;
        }
    }
}

# Close the filehandles
close $r1_fh;
close $r2_fh;
