#!/usr/bin/perl
use strict;
use warnings;

# Check if correct number of arguments are provided
if (@ARGV != 3) {
    die "Usage: $0 <i5 index> <i7 index> <output file prefix>\n";
}

# Get the command-line arguments
my $i5_index = $ARGV[0];
my $i7_index = $ARGV[1];
my $output_prefix = $ARGV[2];

# Open the output file
open(my $out, '>', "${output_prefix}_1.fq") or die "Cannot open output file: $!";

# Process the SAM file from standard input
while (my $line = <STDIN>) {
    chomp $line;

    # Skip header lines
    next if ($line =~ /^\@/);

    my @fields = split("\t", $line);

    # Check the i5 and i7 indexes in the fields array
    if ($fields[11] =~ /$i5_index/ && $fields[13] =~ /$i7_index/) {
        # Write the read to the output FASTQ file
        print $out "@" . $fields[0] . "\n";
        print $out $fields[9] . "\n";
        print $out "+\n";
        print $out $fields[10] . "\n";
    }
}

# Close the output file
close($out);

print "Processing complete.\n";