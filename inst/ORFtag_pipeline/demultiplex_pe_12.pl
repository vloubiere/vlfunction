#!/usr/bin/perl
use strict;
use warnings;

# Check if the correct number of arguments is provided
if ($#ARGV != 1) {
    die "Usage: $0 <regular expression> <output file prefix>\n";
}

my $output1 = "$ARGV[1]_1.fq";
my $output2 = "$ARGV[1]_2.fq";

# Remove existing output files if they exist
if (-e $output1) {
    unlink($output1) or die "Cannot delete existing output file $output1: $!";
}
if (-e $output2) {
    unlink($output2) or die "Cannot delete existing output file $output2: $!";
}

# Open output files for writing
open(my $out1, ">>", $output1) or die "Cannot open output file $output1: $!";
open(my $out2, ">>", $output2) or die "Cannot open output file $output2: $!";

# Variable to track if any paired-end reads are found
my $paired_reads_found = 0;

# Process the input SAM file from STDIN
while (my $line = <STDIN>) {
    chomp $line;

    # Skip header lines
    next if ($line =~ /^\@/);

    my @fields = split(/\t/, $line);

    # Check for paired-end reads using the SAM flag
    my $flag = $fields[1];

    # Check if this is the first mate in a pair (flag 0x40)
    if ($flag & 0x40) {
        if ($fields[11] =~ /$ARGV[0]/) {
            $paired_reads_found = 1;  # Mark that paired-end reads were found

            # Write the first mate to _1.fq in FASTQ format
            say $out1 "@" . $fields[0];
            say $out1 $fields[9];
            say $out1 "+";
            say $out1 $fields[10];
        }
    }
    # Check if this is the second mate in a pair (flag 0x80)
    elsif ($flag & 0x80) {
        if ($fields[11] =~ /$ARGV[0]/) {
            $paired_reads_found = 1;  # Mark that paired-end reads were found

            # Write the second mate to _2.fq in FASTQ format
            say $out2 "@" . $fields[0];
            say $out2 $fields[9];
            say $out2 "+";
            say $out2 $fields[10];
        }
    }
}

# Close the output files
close($out1);
close($out2);

# Check if no paired-end reads were found
if (!$paired_reads_found) {
    die "No paired-end reads detected. Please check if the input SAM file contains paired-end reads.";
}

print "Processing complete.\n";