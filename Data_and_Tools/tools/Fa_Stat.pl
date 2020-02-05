#!/usr/bin/perl
use strict;

if(@ARGV != 2)
{      
        print "Example: perl Fa_Stat.pl main.fa Chr_ID\n";     
        exit(0);
}

my @line;
my $len = 0;
my $LEN = 1;
my $ID;
my $Chr_num = 0;


#######################################
open IN, $ARGV[0] or die "Can't open file";   # genome.fa
open OUT,">$ARGV[1]" or die "Can't open file";     
while(<IN>)
{
    chomp;
    @line = split;
    
    if(/^>/)
    {
    	if($len != 0)
    	{              
           $Chr_num ++;
    	     print OUT "$ID\t$Chr_num\t$len\t$LEN\n";     	
    	}
    	
        $LEN = $LEN + $len;
    	  $ID  = substr($line[0],1,length($line[0])-1);
    	  $len = 0;
    }
    else
    {
       $len = $len + length($line[0]);
    }

}

# last one
       $Chr_num ++;
       print OUT "$ID\t$Chr_num\t$len\t$LEN\n";     		

close(IN);
close(OUT);
print "done!!!\n";
