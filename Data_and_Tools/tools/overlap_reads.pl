#!/usr/bin/perl

if(@ARGV != 2)
{
        print "Usage:   perl overlap_reads.pl <reads_on_ref> <overlap_file> \n";
        print "Example: perl overlap_reads.pl read_on_ref.pos overlap_reads.plot \n";
        exit(0);
}

$b1="";
$b3="";
open IN,"$ARGV[0]" or die "Can't open file";  # reads_on_ref
open OUT,"> $ARGV[1]" or die "Can't write file";
while(<IN>){
    chomp;
    @t=split;
    if(($t[1] eq $b1) && ($t[3] eq $b3)){
        $d=$t[2]-$pos;
         print OUT "$d\n";
    }
    $b1=$t[1];
    $b3=$t[3];
    $pos=$t[2];
}
close(IN);
close(OUT);
