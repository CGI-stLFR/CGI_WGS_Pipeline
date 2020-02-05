#!/usr/bin/perl
use strict;

if(@ARGV != 4)
{
        print "Usage:   perl mate_pairs_distance.pl <Chr_ID> <aln.sam>  <read_len> <output_prefix> \n";	
        print "Example: perl mate_pairs_distance.pl Chr_ID aln.sam 100 mate_pairs_distance \n";
        exit(0);
}


my $readlen  = $ARGV[2];

my @line;
my (%hash_chr,%hash_chr_reverse);
my @ID;
my $barcode = 0;
my $chr_num = 0;
my ($p,$i,$j);
my (@A1,@A2,@B1,@B2,@N);
my ($chr1,$chr2);
my $dist;
my $cigar;
my @C;
my ($pos1,$pos2);
my $n = 0;
my $check = 0;

open OUT,"> $ARGV[3].plot" or die "Can't write file";

############### 
open IN,"$ARGV[0]" or die "Can't open file";  # Chr_ID
while(<IN>)
{
	chomp;
  @line = split;	
  
  $hash_chr{$line[0]} = $line[1];
  $hash_chr_reverse{$line[1]} = $line[0];
}

$chr_num = keys %hash_chr;
close(IN);

############### 
open IN2,"$ARGV[1]" or die "Can't open file";  # aln.sam
while(<IN2>)
{		
	if(!(/^@/))
	{		
		$n++;
		chomp;
    @line = split;
		
        @ID   = split(/\#/,$line[0]);
        
        if($barcode != $ID[1])
        {
        	  if($barcode != 0)
        	  {
               for($p=1;$p<=$chr_num;$p++)
               {              	 
               	 @B1 = sort {$a <=> $b} @{$A1[$p]}[1..$#{$A1[$p]}]; 
               	 @B2 = sort {$a <=> $b} @{$A2[$p]}[1..$#{$A2[$p]}];     
               	 
               	 if(scalar @B1 > 0)
               	 {
               	   for($j=0;$j<@B2;$j++)
               	   {
               	      for($i=0;$i<@B1;$i++)
               	      {
               	      	 if($B1[$i] > $B2[$j])
               	      	 {
               	      	    $dist = $B1[$i] - $B2[$j];
               	      	    print OUT "$dist\n";
               	      	    last;
               	      	 }
               	      }
                 	 } 
                	               	   	 
                 }
               }
 
            }
            
            undef @A1;    undef @A2;
            undef @B1;    undef @B2;            
            for($p=1;$p<=$chr_num;$p++)
            {  $N[$p]    = 0;
            	 $A1[$p][0] = 0;  
            	 $A2[$p][0] = 0;             
            }            
            $barcode = $ID[1];
        }
      
        if($n % 2 == 1)
        {
           if($line[2] ne "*")
           {   	
               $chr1 = $hash_chr{$line[2]};
               
               $cigar =  $line[5];
               @C = split(/S/,$cigar);
               if($C[0]=~ /^\d+$/)  # is number?
               {  $pos1 =  $line[3] - $C[0];
               }
               else
               {  $pos1 =  $line[3];
               }                
           }
           else 
           {   $check ++;
           }
        }
        else   # ($n % 2 == 0)
        {
           if($line[2] ne "*")
           {   	
               $chr2 = $hash_chr{$line[2]};
               
               $cigar =  $line[5];
               @C = split(/S/,$cigar);
               if($C[0]=~ /^\d+$/)  # is number?
               {  $pos2 =  $line[3] - $C[0];
               }
               else
               {  $pos2 =  $line[3];
               } 
           }
           else
           {   $check ++;
           }
               
                
           if($check == 0 && ($chr1 == $chr2))  # chr ne "*"
           {
           	   $N[$chr1] ++;
           	   
           	   if($pos1 < $pos2)
           	   {
           	   	   $A1[$chr1][$N[$chr1]] = $pos1;  
           	   	   $A2[$chr2][$N[$chr2]] = $pos2 + $readlen -1;           	   	
           	   }
           	   else
           	   {
            	   	 $A1[$chr2][$N[$chr2]] = $pos2;  
           	   	   $A2[$chr1][$N[$chr1]] = $pos1 + $readlen -1;             	   	           	   	
           	   }
           	
           }
           
           $check = 0;
        }
  }
}

# last one
               	 if(scalar @B1 > 0)
               	 {
               	   for($j=0;$j<@B2;$j++)
               	   {
               	      for($i=0;$i<@B1;$i++)
               	      {
               	      	 if($B1[$i] > $B2[$j])
               	      	 {
               	      	    $dist = $B1[$i] - $B2[$j];
               	      	    print OUT "$dist\n";
               	      	    last;
               	      	 }
               	      }
                 	 } 
                	               	   	 
                 }
               
close(IN2);

###############

close(OUT);    
   
