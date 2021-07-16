#!/usr/bin/perl
#use strict;
#use warnings;

$ref_genome=$ARGV[0];
$dir=$ARGV[1];

print "Phylogenetic profile ...\n";
$o_table="$dir/Orthologs/Orthologs_table.txt";
open(FILE,$o_table);
my @A=<FILE>;
close FILE;
$A=join("",@A);
$A=~s/\r//ig;
my @B=split(/\n/,$A);
$n=@B;
$ref_genome_t=$ref_genome;
$ref_genome_t=~s/.faa$/.csv/ig;
open(FILE2,$ref_genome_t);
my @A2=<FILE2>;
close FILE2;
$A2=join("",@A2);
$A2=~s/\r//ig;
$A2=~s/\"//ig;
my @B2=split(/\n/,$A2);
$n2=@B2;
for($i=0;$i<$n2;$i++)
{
	my @C2=split(/,/,$B2[$i]);
	if($C2[6] ne "")
	{
		push(@gene_symbols,"$C2[8]\t$C2[6]");
	}
	else
	{
		push(@gene_symbols,"$C2[8]\t$C2[7]");
	}
}

mkdir "$dir/Phylogenetic_Profile";
open(OUT,">$dir/Phylogenetic_Profile/PP_PPIs.txt");
$fig=int($n/100);
$fig2=$fig;
for($i=1;$i<$n;$i++)
{
	if($fig2 eq $i){print ".";$fig2=$fig2+$fig;}
	my @C=split(/\t/,$B[$i]);
	$nC=@C;
	my @grp1=grep{/$C[0]\t/ig}@gene_symbols;
	$ngrp1=@grp1;
	if($ngrp1 ne 0){ my @ppi1=split(/\t/,$grp1[0]); $p1=$ppi1[1];} else { $p1=0; }
	for($j=$i+1;$j<$n;$j++)
	{
		my @D=split(/\t/,$B[$j]);
		$count1=0;$count0=0;
		for($k=1;$k<$nC;$k++)
		{
			if($C[$k] ne 0 && $D[$k] ne 0)
			{
				$count1++;
			}
			if($C[$k] eq 0 && $D[$k] eq 0)
			{
				$count0++;
			}
		}
		$count=$count1+$count0;
		if($count==$nC-1)
		{			
			print OUT "$C[0]\t$D[0]\n";
			my @grp2=grep{/$D[0]\t/ig}@gene_symbols;
			$ngrp2=@grp2;
			if($ngrp2 ne 0){ my @ppi2=split(/\t/,$grp2[0]); $p2=$ppi2[1];} else { $p2=0; }
			if($p1 ne 0 && $p2 ne 0)
			{
				push(@ppis,"$p1\t$p2");
			}
		}
		$count=0;
	}
}
close OUT;
$nppis=@ppis;
open(OUT,">$dir/Phylogenetic_Profile/PP_PPIs_gs.txt");
print OUT "Protein_1\tProtein_2\n";
for($i=0;$i<$nppis;$i++)
{
	print OUT "$ppis[$i]\n";
}
close OUT;
print "\nDone.\n";
