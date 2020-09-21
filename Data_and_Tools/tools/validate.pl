#! /usr/bin/perl
use strict;

my %marker;
open VAR, $ARGV[0] or die $!;
while(<VAR>){
	chomp;
	next if /^#/;
	my ($ref, $pos, $ref_base, $alt_base) = (split /\s+/, $_)[0,1,3,4];
	
	next if /LowQual/; #discard low quality variants
	
#	next unless length $ref_base == 1; #part of indels, especially deletions
	my $info = (split /\s+/, $_)[-1];
	my $depth = (split /:/, $info)[2];
	if ($info=~/1\/2/ or $info=~ /2\/1/ or $info=~/1\|2/ or $info=~/2\|1/){ #het-het
		$alt_base =~ /(\w+),(\w+)/;
		my ($alle_1, $alle_2) = ($1, $2);
		my $new_alle_1 = alle ($alle_1, $ref_base);
		my $new_alle_2 = alle ($alle_2, $ref_base);
		$marker{$ref}{$pos} = "$new_alle_1:$new_alle_2";
	}elsif ($info=~/0\/1/ or $info=~/0\|1/ or $info=~/1\|0/){ #het-ref
		my $new_alle = alle ($alt_base, $ref_base);
		$ref_base =~ /^([A|T|G|C])/;
		$ref_base = $1;
		$marker{$ref}{$pos} = "$ref_base:$new_alle";
	}elsif ($info=~/1\/1/ or $info=~ /1\|1/){ #homo
		my $new_alle = alle ($alt_base, $ref_base);
		$ref_base =~ /^([A|T|G|C])/;
        $ref_base = $1;
        #$marker{$ref}{$pos} = "$ref_base:$new_alle";
	}else{ #missed other type?
		warn "missed other type:\n$_\n";
	}
}

sub alle {
        my ($alt, $ref) = @_;
        my $alt_len = length $alt;
        my $ref_len = length $ref;
        if ($alt_len < $ref_len){
                my $indel_len = $ref_len - $alt_len;
                $alt = "D$indel_len";
        }elsif ($alt_len > $ref_len){
                my $indel_len = $alt_len - $ref_len;
                $alt = "I$indel_len";
        }else{
                $alt =~ /^([A|T|G|C])/; #GA G,AA
                $alt = $1;
        }
        return $alt;
}

my @files =`find $ARGV[1] -name "*.log"`;
chomp @files;
my %Len;
my %Num;
my %hash1;
my %hash2;
my %Var;
my %VarNum;
my %info;
my ($Total_Len, %Total_Num, $Total_Var, %Total_VarNum, %Total_info);
foreach my $file (@files){
	my $basename = (split /\//, $file)[-1];
	$basename =~ /(.+)(\.log)/;
	my $chr = $1;
	$/ = "break";
	<IN>;
	my $total_len;
	my %hash;
	open IN, $file or die $!;
	while(<IN>){
		chomp;
		my @lines = split /\n/, $_;
		my $line1 = $lines[1];
		my $line2 = $lines[2];
		my @pos1 = split /;/, $line1;
		my @pos2 = split /;/, $line2;
		next unless @pos1 >1;
		my $min;
		my $max;
		my $var_num;
		for (my $i=0; $i<@pos1; $i++){
			$var_num++;
			my $pos_1 = $pos1[$i];
			my $pos_2 = $pos2[$i];
			my ($pos,$alle);
			$pos_1 =~ /(.+)_(\w+)_(\w+)/;
			($pos,$alle) = ($2,$3);
			$hash1{$chr}{$pos} =$alle;
			$pos_2 =~ /(.+)_(\w+)_(\w+)/;
			($pos,$alle) = ($2,$3);
			$hash2{$chr}{$pos} = $alle;
			$min = $pos if $i==0;
			$min = $pos if $pos <$min;
			$max = $pos if $pos >$max;
		}
		my $len = $max- $min;
		my $raw_num;
		foreach my $pos (sort {$a<=>$b} keys %{$marker{$chr}}){
			if ($pos >= $min and $pos <= $max){
				$raw_num++;
			}elsif ($pos >$max){
				last;
			}
		}
        if ($raw_num == 0){
            warn "$chr has log variants but none in VCF\n";
            next;
        }
		my $adjust_len = $len * ($var_num/$raw_num);
		$Len{$chr} += $len;
		$Num{$chr}{$len}++;
	
		$Var{$chr} += $var_num;

		$info{$chr}{$adjust_len} += $var_num;

		$Total_Len += $len;
		$Total_Num{$len}++;
		$Total_VarNum{$adjust_len} += $var_num;
		
		$Total_Var+= $var_num;
	}
	close IN;
	$/ = "\n";
}

my $total_ac_len;
foreach my $len (sort {$b<=>$a} keys %Total_Num){
	$total_ac_len += $len * $Total_Num{$len};
	if ($total_ac_len > 0.5*$Total_Len){
		print "N50: \t",$len, "\n";
		last;
	}
}

my $total_ac_var;
foreach my $adjust_len (sort {$b<=>$a} keys %Total_VarNum){
	$total_ac_var += $Total_VarNum{$adjust_len};
	if ($total_ac_var >0.5*$Total_Var){
		print "AN50:\t",$adjust_len, "\n";
		last;
	}
}
foreach my $var_num (sort {$b<=>$a} keys %Total_VarNum){
#	print "max block snp frac:",$var_num/$Total_Var,"\n";
#	print "$var_num\t$Total_Var\n";
	last;
}

open CHR, ">$ARGV[1]/stat.by.chromosome" or die;
print CHR "N50 by chromosome:\n";
foreach my $chrn (sort {$a<=>$b} keys %marker){
	my $ac_len;
	foreach my $len (sort {$b<=>$a} keys %{$Num{$chrn}}){
		my $num = $Num{$chrn}{$len};
		$ac_len += $len * $num;
		if ($ac_len >0.5*$Len{$chrn}){
			print CHR "$chrn:",$len, "\n";
			last;
		}
	}
}

print CHR "AN50 by chromosome:\n";
foreach my $chrn (sort {$a <=> $b} keys %marker){
	my $ac_var;
	foreach my $adjust_len (sort {$b<=>$a} keys %{$info{$chrn}}){
		$ac_var += $info{$chrn}{$adjust_len};
        if ($ac_var >0.5*$Var{$chrn}){
	        print CHR "$chrn:",$adjust_len, "\n";
            last;
        }
	}
}
close CHR;
$/ = "\n";

