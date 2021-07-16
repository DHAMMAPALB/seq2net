#!/usr/bin/perl
#use strict;
#use warnings;

$file=$ARGV[0];
$ref_genome=$ARGV[1];
$dir=$ARGV[2];

print "Gene neighborhood ...\n";
$line=`grep -n $ref_genome $file`;
my @line_n=split(/:/,$line);
$ref_n=$line_n[0]-1;
open(FILE,$file);
my @A=<FILE>;
close FILE;
$A=join("",@A);
$A=~s/\r//ig;
$A=~s/.faa\n/.csv\n/ig;
my @B=split(/\n/,$A);
@B=sort(@B);
$n=@B;
$o_table="$dir/Orthologs/Orthologs_table.txt";
open(FILE2,$o_table);
my @A2=<FILE2>;
close FILE2;
$A2=join("",@A2);
$A2=~s/\r//ig;
my @B2=split(/\n/,$A2);
$n2=@B2;
mkdir "$dir/Gene_Neighborhood";
my @ppis;
open(OUT,">$dir/Gene_Neighborhood/GN_PPIs.txt");
$fig=int($n2/100);
$fig2=$fig;
for($i=1;$i<$n2;$i++)
{
	if($fig2 eq $i){print ".";$fig2=$fig2+$fig;}
	my @C2=split(/\t/,$B2[$i]);
	$nC2=@C2;
	for($j=$i+1;$j<$n2;$j++)
	{
		my @D2=split(/\t/,$B2[$j]);
		$count=0;
		my @grp1=`grep $C2[0] $B[$ref_n]`;
		my @grp2=`grep $D2[0] $B[$ref_n]`;
		$grp1[0]=~s/\"//ig;
		$grp2[0]=~s/\"//ig;
		my @cor1=split(/,/,$grp1[0]);
		my @cor2=split(/,/,$grp2[0]);
		$diff=abs($cor2[2]-$cor1[3]);
		if($cor1[6] ne "") { $p1=$cor1[6]; } else { $p1=$cor1[7]; }
		if($cor2[6] ne "") { $p2=$cor2[6]; } else { $p2=$cor2[7]; }
		if($diff<=300)
		{
			$count++;
		}
		for($k=0;$k<$nC2;$k++)
		{
			if($k!=$ref_n)
			{
				if($C2[$k+1] ne 0 && $D2[$k+1] ne 0)
				{
					my @E2=split(/,/,$C2[$k+1]);
					$nE2=@E2;
					my @F2=split(/,/,$D2[$k+1]);
					$nF2=@F2;
					for($l=0;$l<$nE2;$l++)
					{
						my @grp1=`grep $E2[$l] $B[$k]`;
						$ngrp1=@grp1;
						if($ngrp1==0){last};
						my @cor1=split(/,/,$grp1[0]);
						$count2=0;
						for($m=0;$m<$nF2;$m++)
						{
							my @grp2=`grep $F2[$m] $B[$k]`;
							$ngrp2=@grp2;
							if($ngrp2==0){last};
							my @cor2=split(/,/,$grp2[0]);
							$diff=abs($cor2[2]-$cor1[3]);
							if($diff<=300)
							{
								$count++;
								$count2++;
								last;
							}
						}
						if($count2>0)
						{
							last;
						}
					}
				}
			}
			if($count>1)
			{
				last;
			}
		}
		if($count>1)
		{
			print OUT "$C2[0]\t$D2[0]\n";
			push(@ppis,"$p1\t$p2");
		}
		$count=0;
	}
}
close OUT;
$nppis=@ppis;
open(OUT,">$dir/Gene_Neighborhood/GN_PPIs_gs.txt");
print OUT "Protein_1\tProtein_2\n";
for($i=0;$i<$nppis;$i++)
{
	print OUT "$ppis[$i]\n";
}
close OUT;
print "\nDone.\n";
