#!usr/bin/perl;
#use strict;
#use warnings;
use List::Util qw(uniq);

$ref_genome=$ARGV[0];
$dir=$ARGV[1];

print  "Gene fusion ...\n";
$o_table="$dir/OrthoParalogs/Orthologs_table.txt";
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

my @genomes=split(/\t/,$B[0]);
$fig=int($n/100);
$fig2=$fig;
for($i=1;$i<$n;$i++)
{
	if($fig2 eq $i){print ".";$fig2=$fig2+$fig;}
	my @C=split(/\t/,$B[$i]);
	$nC=@C;
	for($j=1;$j<$nC;$j++)
	{
		if($C[$j] ne 0)
		{
			my @orthos;
			for($k=1;$k<$n;$k++)
			{
				my @D=split(/\t/,$B[$k]);
				push(@orthos,$D[$j]);
			}
			$northos=@orthos;
#			print @orthos;getc();
			for($k=0;$k<$northos;$k++)
			{
				if($k ne $i-1)
				{
					if($C[$j] eq $orthos[$k])
					{
						my @p2=split(/\t/,$B[$k+1]);
						push(@ppis,"$C[0]-$p2[0]\t$orthos[$k]\t$genomes[$j]");
					}
					else
					{
						my @x1=split(/\, /,$C[$j]);
						my @x2=split(/\, /,$orthos[$k]);
						$nx1=@x1;$nx2=@x2;
						for($l=0;$l<$nx1;$l++)
						{
							for($m=0;$m<$nx2;$m++)
							{
								if($x1[$l] eq $x2[$m])
								{
									my @p3=split(/\t/,$B[$k+1]);
									push(@ppis,"$C[0]-$p3[0]\t$x2[$m]\t$genomes[$j]");				
								}
							}
						}
					}
				}
			}
		}
	}
}
@ppis=sort(@ppis);
@ppis=uniq(@ppis);
$nppis=@ppis;
$ngrp=0;
for($i=0;$i<$nppis;$i++)
{
	my @ss=split(/\t/,$ppis[$i]);
	my @tt=split(/-/,$ss[0]);
	$C8="$tt[1]-$tt[0]\t$ss[1]";
	my @grp=grep{/$C8/i}@ppis[$i+1..$nppis-1];
	$ngrp=@grp;
	if($ngrp==0)
	{
		push(@ppis2,$ppis[$i]);
	}
	$ngrp=0;
}
@ppis2=sort(@ppis2);
$nppis2=@ppis2;

mkdir "$dir/Gene_Fusion";
open(OUT,">$dir/Gene_Fusion/Composite_proteins.txt");
print OUT "Component_proteins\tComposite_protein(s)\tOrganism\n";
for($i=0;$i<$nppis2;$i++)
{
	print OUT "$ppis2[$i]\n";
}
close OUT;

open(FILE3,"$dir/OrthoParalogs/RBBH_p.txt");
my @A3=<FILE3>;
close FILE3;
$A3=join("",@A3);
$A3=~s/\r//ig;
my @B3=split(/\n/,$A3);
$nB3=@B3;
$ngrp1=0;$ngrp2=0;
open(OUT,">$dir/Gene_Fusion/GF_PPIs.txt");
for($i=0;$i<$nppis2;$i++)
{
	my @ss=split(/\t/,$ppis2[$i]);
	my @tt=split(/-/,$ss[0]);
	my @grp1=grep{/$tt[0]\t$tt[1]/i}@B3;
	$ngrp1=@grp1;
	my @grp2=grep{/$tt[1]\t$tt[0]/i}@B3;
	$ngrp2=@grp2;
	if($ngrp1==0 && $ngrp2==0)
	{
		print OUT "$tt[1]\t$tt[0]\n";
		push(@ppis3,"$tt[1]\t$tt[0]");
	}
	$ngrp1=0;$ngrp2=0;	
}
close OUT;
$nppis3=@ppis3;
open(OUT,">$dir/Gene_Fusion/GF_PPIs_gs.txt");
$ngrp1=0;$ngrp2=0;
print OUT "Protein_1\tProtein_2\n";
for($i=0;$i<$nppis3;$i++)
{
	my @C3=split(/\t/,$ppis3[$i]);
	my @grp1=grep{/$C3[0]/i}@gene_symbols;
	$ngrp1=@grp1;
	my @grp2=grep{/$C3[1]/i}@gene_symbols;
	$ngrp2=@grp2;
	if($ngrp1!=0 && $ngrp2!=0)
	{
		my @p1=split(/\t/,$grp1[0]);
		my @p2=split(/\t/,$grp2[0]);
		print OUT "$p1[1]\t$p2[1]\n";
	}
}
close OUT;
print "\nDone.\n";
