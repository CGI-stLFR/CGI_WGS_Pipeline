#!/usr/bin/perl

if(@ARGV != 2)
{
        print "Usage:   perl mark_dups_txt.pl <input_metrics> <output_metrics> \n";
        print "Example: perl mark_dups_txt.pl hg001_dedup_metrics.txt hg001_dedup_metrics2.txt \n";
        exit(1);
}

$n = 0;
$j = -1;
open IN,"$ARGV[0]" or die "Can't open file";  # input_metrics
open OUT,"> $ARGV[1]" or die "Can't write file";
while(<IN>){
    $n++;
    if($n==2 || $n==3){
        @A;
        $j++;
        chomp;
        @t = split;
        for($i = 0; $i < @t; $i++){
                $A[$j][$i] = $t[$i];
        }
    }
}
for($i = 0; $i < @t; $i++){
    print OUT "$A[0][$i]\t$A[1][$i+1]\n";
}

close(IN);
close(OUT);
