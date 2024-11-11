#!/usr/bin/perl
use strict;
use warnings;

# Usage: perl script_name.pl tar_file_path barcode i5 fq_prefix

# Get arguments
my ($tar_file, $barcode, $i5, $fq_prefix) = @ARGV;

# Open output files with .gz compression
open(my $out1, "| gzip > ${fq_prefix}_1.fq.gz") or die "Cannot open output file 1: $!";
open(my $out2, "| gzip > ${fq_prefix}_2.fq.gz") or die "Cannot open output file 2: $!";

# Open tar files for read 1 and read 2
open(my $r1, "tar -xOzf $tar_file --wildcards '*_R1_*.fastq.gz' | zcat |") or die "Cannot open tar file for R1: $!";
open(my $r2, "tar -xOzf $tar_file --wildcards '*_R2_*.fastq.gz' | zcat |") or die "Cannot open tar file for R2: $!";

# Create regex patterns for barcode and i5 with flexible matching
my $barcode_regex = qr/:N:0:$barcode/;         # Matches barcode part after N:0:
my $i5_regex = qr/\+$i5/;                      # Matches i5 part after '+'

# Debug counters
my $total_reads = 0;
my $matched_reads = 0;

# Read paired-end files line by line
while (my $r1_line = <$r1>) {
    my $r2_line = <$r2>;

    # Count total reads processed
    $total_reads++;

    # Process the header line to check if it matches barcode and i5
    if ($r1_line =~ /:N:0:$barcode.*?\+$i5/) {  # Adjusted regex to match any characters after barcode and i5
        # Count matched reads
        $matched_reads++;

        # Capture only the part of the read ID before the first space (to exclude barcode information)
        my ($read_id1) = $r1_line =~ /^@(\S+)/;
        my ($read_id2) = $r2_line =~ /^@(\S+)/;

        # Write matching read pairs to output files with trimmed headers
        print $out1 "\@$read_id1\n";
        print $out1 scalar(<$r1>);  # sequence
        print $out1 scalar(<$r1>);  # +
        print $out1 scalar(<$r1>);  # quality

        print $out2 "\@$read_id2\n";
        print $out2 scalar(<$r2>);  # sequence
        print $out2 scalar(<$r2>);  # +
        print $out2 scalar(<$r2>);  # quality
    } else {
        # Skip lines for non-matching reads
        <$r1> for (1..3);
        <$r2> for (1..3);
    }
}

# Close files
close $r1;
close $r2;
close $out1;
close $out2;

# Print debugging info
print "Total reads processed: $total_reads\n";
print "Total reads matched: $matched_reads\n";
