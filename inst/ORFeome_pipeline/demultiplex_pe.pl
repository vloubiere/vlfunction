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
my $mate1_file = "${output_prefix}_1.fq.gz";
my $mate2_file = "${output_prefix}_2.fq.gz";

# Quick check to ensure the SAM file contains paired-end reads
my $contains_paired_end = 0;
my $line_count = 0;
while (<STDIN>) {
    chomp;
    next if /^@/; # Skip header lines
    my @fields = split("\t");
    my $flag = $fields[1];

    # Check if the read is part of a pair (flag 0x1 is set)
    if ($flag & 0x1) {
        $contains_paired_end = 1;
        last;
    }

    last if ++$line_count >= 100; # Stop checking after 100 lines
}

# Error out if no paired-end reads are found
unless ($contains_paired_end) {
    die "Error: Input SAM file does not appear to contain paired-end reads.\n";
}

# Open output files for writing
my $mate1_out = IO::Zlib->new($mate1_file, "wb") or die "Could not open $mate1_file: $!\n";
my $mate2_out = IO::Zlib->new($mate2_file, "wb") or die "Could not open $mate2_file: $!\n";

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
        # Determine mate (based on flag)
        if ($flag & 0x40) {
            # Mate 1
            print $mate1_out "\@$read_id\n$seq\n+\n$qual\n";
        } elsif ($flag & 0x80) {
            # Mate 2
            print $mate2_out "\@$read_id\n$seq\n+\n$qual\n";
        }
    }
}

# Close output files
$mate1_out->close;
$mate2_out->close;

print "Done! Output written to $mate1_file and $mate2_file\n";
