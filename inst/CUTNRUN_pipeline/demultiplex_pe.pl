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

# Open the output files
open(my $out1, '>', "${output_prefix}_1.fq") or die "Cannot open output file 1: $!";
open(my $out2, '>', "${output_prefix}_2.fq") or die "Cannot open output file 2: $!";

# Process the SAM file from standard input
while (my $line = <STDIN>) {
    chomp $line;

    # Skip header lines
    next if ($line =~ /^\@/);

    my @fields = split("\t", $line);

    # Check if this is the first mate
    if ($fields[1] & 0x40) {
        # Check i5 and i7 indexes
        if ($fields[11] =~ /$i5_index/ && $fields[13] =~ /$i7_index/) {
            # Write the first mate to _1.fq in FASTQ format
            print $out1 "@" . $fields[0] . "\n";
            print $out1 $fields[9] . "\n";
            print $out1 "+\n";
            print $out1 $fields[10] . "\n";

            # Read the next line for the second mate
            my $second_mate = <STDIN>;
            chomp $second_mate;
            my @second_fields = split("\t", $second_mate);

            # Write the second mate to _2.fq in FASTQ format
            print $out2 "@" . $second_fields[0] . "\n";
            print $out2 $second_fields[9] . "\n";
            print $out2 "+\n";
            print $out2 $second_fields[10] . "\n";
        }
    }
}

# Close the output files
close($out1);
close($out2);

print "Processing complete.\n";