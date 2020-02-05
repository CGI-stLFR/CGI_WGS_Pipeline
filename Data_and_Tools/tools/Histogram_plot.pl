#!/usr/bin/perl
use strict;

# 按区间统计一列数据的直方图

if(@ARGV != 5)
{
	print "Usage:   perl Histogram_plot.pl <data column> <min> <step> <max> <output>\n";
	print "Example: perl Histogram_plot.pl InsertSize 0 20 1000 InsertSize.plot\n";			
	exit(0);
}

open IN1,$ARGV[0] or die "Can't open file";
open OUT,">".$ARGV[4] or die "Can't write file";
my $min  = $ARGV[1];
my $step = $ARGV[2];
my $max  = $ARGV[3];
my $total_num = $max/$step;
my $i = 0;
my @index;
my @line;
my @frequency;


   while($i <= $total_num)
   {
      $frequency[$i] = 0;
      $index[$i] = $i*$step;
   	  $i++;   	
   }


   while(<IN1>)
   { 
   	  chomp;
  	  @line=split;  	    
      $frequency[int($line[0]/$step)] ++ ; #关键语句，分类
 	    
   }
    

   $i = 0;
   while($i <= $total_num)
   {
   	  if($index[$i] >= $min)
   	  {  print OUT "$index[$i]\t$frequency[$i]\n";  	
   	  }
   	  $i ++;
   }   

close IN1;
close OUT;

