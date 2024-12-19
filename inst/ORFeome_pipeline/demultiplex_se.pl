#!/usr/bin/perl
use strict;
use warnings;
use IO::Zlib;

# Check arguments
if (@ARGV != 2) {
    die "Usage: $0 <i7 index> <output file prefix> < <SAM file>\n";
}

# Parse arguments
my $i7_index = $ARGV[0];
my $output_prefix = $ARGV[1];
my $output_file = "${output_prefix}_1.fq.gz";

# Quick check to ensure the SAM file contains single-end reads
my $contains_single_end = 0;
my $line_count = 0;
while (<STDIN>) {
    chomp;
    next if /^@/; # Skip header lines
    my @fields = split("\t");
    my $flag = $fields[1];

    # Check if the read is single-end (flag 0x1 is NOT set)
    if (!($flag & 0x1)) {
        $contains_single_end = 1;
        last;
    }

    last if ++$line_count >= 100; # Stop checking after 100 lines
}

# Error out if no single-end reads are found
unless ($contains_single_end) {
    die "Error: Input SAM file does not appear to contain single-end reads.\n";
}

# Reopen STDIN to process the file again
seek(STDIN, 0, 0);

# Open output file for writing
my $output_fh = IO::Zlib->new($output_file, "wb") or die "Could not open $output_file: $!\n";

# Process the SAM input line by line
while (<STDIN>) {
    chomp;
    next if /^@/; # Skip header lines

    my @fields = split("\t");
    my $read_id = $fields[0];
    my $flag = $fields[1];
    my $seq = $fields[9];
    my $qual = $fields[10];

    # Automatically detect the column containing the "BC:Z:" tag
    my ($bc_tag) = grep { /^BC:Z:/ } @fields[11..$#fields];
    next unless $bc_tag; # Skip if no BC:Z tag is found

    # Check if the tag matches the given i7 index
    if ($bc_tag =~ /^BC:Z:$i7_index/) {
        # Write the read to the output file
        print $output_fh "\@$read_id\n$seq\n+\n$qual\n";
    }
}

# Close output file
$output_fh->close;

print "Done! Output written to $output_file\n";
