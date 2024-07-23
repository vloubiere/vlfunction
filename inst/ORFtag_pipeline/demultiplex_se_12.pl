#!/usr/bin/perl
use strict;
use warnings;

if ($#ARGV != 1) {
	print "Args should specify 1/ regexpr and 2/ output file prefix\n";
	exit;
}

my $output1= "$ARGV[1]\_1.fq";

if (-e $output1) 
{
    unlink($output1);
}

open(my $out1, ">>", $output1) or die "...";

while(<STDIN>)
{
    my @cur = split(/\t/, $_);
    if($cur[11] =~ $ARGV[0])
    {
      say $out1 "@" . $cur[0];
      say $out1 $cur[9];
      say $out1 "+";
      say $out1 $cur[10];
    }
}
close($out1);