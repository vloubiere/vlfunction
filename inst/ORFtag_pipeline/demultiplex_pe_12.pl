#!/usr/bin/perl
use strict;
use warnings;

if ($#ARGV != 1) {
	print "Args should specify 1/ regexpr and 2/ output file prefix\n";
	exit;
}

my $output1= "$ARGV[1]\_1.fq";
my $output2= "$ARGV[1]\_2.fq";

if (-e $output1) 
{
    unlink($output1);
}
if (-e $output2) 
{
    unlink($output2);
}

open(my $out1, ">>", $output1) or die "...";
open(my $out2, ">>", $output2) or die "...";

my $mate= 0;
while(<STDIN>)
{
    my @cur = split(/\t/, $_);
    if($mate == 0)
    {
        if($cur[11] =~ $ARGV[0])
        {
            say $out1 "@" . $cur[0];
            say $out1 $cur[9];
            say $out1 "+";
            say $out1 $cur[10];
            $mate= 1;
        }
    }else{
            say $out2 "@" . $cur[0];
            say $out2 $cur[9];
            say $out2 "+";
            say $out2 $cur[10];
        $mate= 0;
    }
}
close($out1);
close($out2);