#!usr/bin/perl
use strict;
my $usage =<<USAGE;
Function: compute the rate which dulicate reads from the same barcode.
         (1)Duplicate barcode rate in duplicate reads: Duplicate barcode number/ Duplicate reads number.

example: perl $0 <L0.sort.rmdup.bam>
USAGE
die $usage if(@ARGV != 2);

my $file = @ARGV[0];
my $path = @ARGV[1];

open (IN,"$file") or die ("can't open the $file!\n");
open (OUT,">$path/PE_dup_reads") or die ("can't creat the file!\n");
open (OUT2, ">$path/dup_info") or die ("can't create the file!\n");
open (OUT3, ">$path/duplicate_rate") or die("can't create the file!\n");

my %locate_hash;
my ($key,$value);
my $dup_num = 0;
my $dup_PE = 0;

print "step1: start to create the hash";
while(<IN>){
  chomp;
  next if(/\@HD|\@SQ|\@PG/);
  my @line = split;
  if($line[1]>=1024){
    my @ID = split /#/,$line[0];
    my $chrID = "$line[2]:$line[3]";
    if(!exists($locate_hash{$chrID}{$ID[1]}{$ID[0]})){
      $locate_hash{$chrID}{$ID[1]}{$ID[0]} = $line[9];
    }else{
      $dup_PE ++;
      print OUT "$chrID\t$ID[1]\t$ID[0]\t1:$locate_hash{$chrID}{$ID[1]}{$ID[0]}\t2:$line[9]\n";
    }
    $dup_num ++;
  }
}

print OUT2 "#duplicate location\tbar_in_reads_rate\tdupbar_in_reads_rate\tun_snp_in_bar_rate\n";
print "step2:caculate the duplicate rate\n";

my $bar_total = 0;
my $reads_total = 0;
my $dup_reads = 0;
my $un_snp_total = 0;
my $bar_reads_total =0;

#print OUT2 "dupplicate location\tbar_num\treads_num\tdup_bar\tun_snp\n";
while(($key,$value)= each %locate_hash){
  my $bar_num = keys%$value;
  my $total_reads = 0;
  my $false_dup = 0;
  my $un_snp_num = 0;
  while(my($k1,$v1) = each %$value){
    my $reads_num = keys%$v1;
    $total_reads = $total_reads + $reads_num;
    if($reads_num > 1){
      $dup_reads = $dup_reads + $reads_num;
      my %seq_hash;
      while(my($k2,$v2) = each %$v1){
        if(!exists($seq_hash{$v2})){
           $seq_hash{$v2} = 1;
        }else{
           $seq_hash{$v2} += 1;
        }
      }
      while(my($k3,$v3) = each %seq_hash){
        if($v3 > 1){
         $un_snp_num += $v3;
        }
      }  
    }else{
      $false_dup ++;
    }
  }
 my $bar_in_reads_rate = $bar_num/$total_reads;
 my $dupbar_in_reads_rate = 1 - ($false_dup / $total_reads);
 my $bar_reads = $total_reads - $false_dup;
 my $un_snp_in_bar_rate = $un_snp_num;
 if($bar_reads != 0){
    $un_snp_in_bar_rate = $un_snp_num / $bar_reads;
 }
 $un_snp_total = $un_snp_total + $un_snp_num;
 $bar_total = $bar_total + $bar_num;
 $reads_total = $reads_total + $total_reads;
 $bar_reads_total = $bar_reads_total + $bar_reads;
 print OUT2 "$key\t$bar_in_reads_rate\t$dupbar_in_reads_rate\t$un_snp_in_bar_rate\n";
}
my $dup_bar_rate = $bar_total/$reads_total;
my $real_dup_bar_rate = $dup_reads / $reads_total;
my $un_snp_dup_rate = $un_snp_total / $bar_reads_total;

print OUT3 "duplicate reads number:$reads_total\n";
print OUT3 "duplicate barcode number:$bar_total\n";
#print OUT3 "total duplicate in different reads:$dup_reads\n";
print OUT3 "bar_rate in dup reads:$dup_bar_rate\n";
#print OUT3 "real_dup_bar_rate in dup reads:$real_dup_bar_rate\n";
#print OUT3 "un_snp_dup_rate in dup_bar:$un_snp_dup_rate\n";

close IN;
close OUT;
close OUT2;
close OUT3;
