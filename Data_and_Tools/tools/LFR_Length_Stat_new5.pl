#!/usr/bin/perl
use strict;

# 处理 CIGAR 串，修正 softclip 造成的 reads 比对起始位置偏移

if(@ARGV != 6)
{
        print "Usage:   perl LFR_Length_Stat.pl <Chr_ID> <aln.sam>  <MAPQ_limit> <read_len> <genome_len> <output_prefix> \n";	
        print "Example: perl LFR_Length_Stat.pl Chr_ID aln.sam 60 100 3137161264  Reads_on_Ref \n";
        exit(0);
}

my $MAPQ     = $ARGV[2];
my $readlen  = $ARGV[3];
my $genome   = $ARGV[4];
my @line;
my (%hash_chr,%hash_chr_reverse);
my @ID;
my $barcode = 0;
my $chr_num = 0;
my ($i,$j);
my (@A,@B,@N);
my $chr;
my $start;
my $cigar;
my @C;
my $pos;


open OUT1,"> $ARGV[5]" or die "Can't write file";

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
		chomp;
    @line = split;
		
		if($line[4] >= $MAPQ)
		{
        @ID   = split(/\#/,$line[0]);
        
        if($barcode != $ID[1])
        {
        	  if($barcode != 0)
        	  {
               for($i=1;$i<=$chr_num;$i++)
               {              	 
               	 @B = sort {$a <=> $b} @{$A[$i]}[1..$#{$A[$i]}];    
               	 
               	 if(scalar @B > 0)
               	 {
               	   for($j=0;$j<@B;$j++)
               	   {
               	      $start = $B[$j];               	   	  
               	   	  print OUT1 "$readlen\t$hash_chr_reverse{$i}\t$start\t$barcode\n";
                 	 } 
                	               	   	 
                 }
               }
 
            }
            
            undef @A;
            for($i=1;$i<=$chr_num;$i++)
            {  $N[$i]    = 0;
            	 $A[$i][0] = 0;               
            }            
            $barcode = $ID[1];
        }
      
        if($line[2] ne "*")
        {   	
               $chr = $hash_chr{$line[2]};
               $N[$chr] ++;
               
               $cigar =  $line[5];
               @C = split(/S/,$cigar);
               if($C[0]=~ /^\d+$/)  # is number?
               {  $pos =  $line[3] - $C[0];
               }
               else
               {  $pos =  $line[3];
               } 
               	
               $A[$chr][$N[$chr]] = $pos; 
        }
    } 
  }
}

# last one
               for($i=1;$i<=$chr_num;$i++)
               {              	 
               	 @B = sort {$a <=> $b} @{$A[$i]}[1..$#{$A[$i]}];    
               	 
               	 if(scalar @B > 0)
               	 {
               	   for($j=0;$j<@B;$j++)
               	   {
               	      $start = $B[$j];               	   	  
               	   	  print OUT1 "$readlen\t$hash_chr_reverse{$i}\t$start\t$barcode\n";
                 	 } 
                	               	   	 
                 }
               }
               
close(IN2);

###############

close(OUT1);    
   
