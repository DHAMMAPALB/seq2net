#!/usr/bin/perl
#use strict;
#use warnings;
use List::Util qw(uniq);

$ref_genome=$ARGV[0];
$ppi_genome=$ARGV[1];
$ppi_list=$ARGV[2];
$e_value=$ARGV[3];
$coverage=$ARGV[4];
$dir=$ARGV[5];

print "Interlog ...\n";
mkdir "$dir/Interlog";
`makeblastdb -in $ref_genome -dbtype prot`;
`makeblastdb -in $ppi_genome -dbtype prot`;
my @out_q=split(/\//,$ref_genome);$nout_q=@out_q;
my @out_s=split(/\//,$ppi_genome);$nout_s=@out_s;
$out_f="$dir/Interlog/$out_q[$nout_q-1]-vs-$out_s[$nout_s-1]";
`blastp -query $ref_genome -subject $ppi_genome -out $out_f -evalue $e_value -outfmt \"6 qseqid sseqid pident evalue bitscore\"`;
$out_r="$dir/Interlog/$out_s[$nout_s-1]-vs-$out_q[$nout_q-1]";
`blastp -query $ppi_genome -subject $ref_genome -out $out_r -evalue $e_value -outfmt \"6 qseqid sseqid pident evalue bitscore\"`;

open(FILE3,"$out_f");my @A3=<FILE3>;close FILE3;$nA3=@A3;
open(FILE4,"$out_r");my @A4=<FILE4>;close FILE4;$nA4=@A4;

$out_rbbh="$dir/Interlog/RBBH_i.txt";
$ngrp=0;
open(OUT,">$out_rbbh");
for($j=0;$j<$nA3;$j++)
{
	my @B3=split(/\t/,$A3[$j]);my @grp=grep{/$B3[1]\t$B3[0]\t/ig}@A4;$ngrp=@grp;
	if($ngrp!=0)
	{
		my @C3=split(/\t/,$grp[0]);
		if($B3[2]>=$coverage && $C3[2]>=$coverage){ print OUT "$B3[0]\t$B3[1]\n"; push(@pp,"$B3[0]\t$B3[1]"); } $ngrp=0;
	}
}
close OUT;
@pp=uniq(@pp);
$npp=@pp;
for($i=0;$i<$npp;$i++)
{
	my @ss=split(/\t/,$pp[$i]);
	push(@RV_i,$ss[0]);
	push(@ortho_i,$ss[1]);
}
$nnn=@RV_i;
$out_rbbh_i_t="$dir/Interlog/RBBH_i_table.txt";
open(OUT,">$out_rbbh_i_t");
print OUT "$RV_i[0]\t";
for($l=0,$k=1;$l<$nnn;$l++,$k++)
{
	if($RV_i[$l] eq $RV_i[$k]) {print OUT "$ortho_i[$l],";}
	else { print OUT "$ortho_i[$l]\n$RV_i[$k]\t";}
}
close OUT;

open(FILE5,$ppi_list);
my @A5=<FILE5>;
close FILE5;
$A5=join("",@A5);
$A5=~s/\r//ig;
my @B5=split(/\n/,$A5);
$fig=int($nnn/100);
$fig2=$fig;
for($i=0;$i<$nnn;$i++)
{
	if($fig2 eq $i){print ".";$fig2=$fig2+$fig;}
	my @grp1=grep{/$ortho_i[$i]/i}@B5;
	if($grp1[0] ne "")
	{
		my @C5=split(/\t/,$grp1[0]);
		my @grp2=grep{/$C5[0]/i}@pp;
		my @grp3=grep{/$C5[1]/i}@pp;
		if($grp2[0] ne "" && $grp3[0] ne "")
		{
			my @D5=split(/\t/,$grp2[0]);
			my @E5=split(/\t/,$grp3[0]);
			if($D5[0] ne $E5[0])
			{
				push(@ppis,"$D5[0]\t$E5[0]");
			}
		}
	}
}
@ppis=uniq(@ppis);
$nppis=@ppis;
open(OUT,">$dir/Interlog/IN_PPIs.txt");
for($i=0;$i<$nppis;$i++)
{
	print OUT "$ppis[$i]\n";
}
close OUT;
$ref_genome_t=$ref_genome;
$ref_genome_t=~s/.faa$/.csv/ig;
open(FILE6,$ref_genome_t);
my @A6=<FILE6>;
close FILE6;
$A6=join("",@A6);
$A6=~s/\r//ig;
$A6=~s/\"//ig;
my @B6=split(/\n/,$A6);
$n6=@B6;
for($i=0;$i<$n6;$i++)
{
	my @C6=split(/,/,$B6[$i]);
	if($C6[6] ne "")
	{
		push(@gene_symbols,"$C6[8]\t$C6[6]");
	}
	else
	{
		push(@gene_symbols,"$C6[8]\t$C6[7]");
	}
}
open(OUT,">$dir/Interlog/IN_PPIs_gs.txt");
$ngrp1=0;$ngrp2=0;
print OUT "Protein_1\tProtein_2\n";
for($i=0;$i<$nppis;$i++)
{
	my @C7=split(/\t/,$ppis[$i]);
	my @grp1=grep{/$C7[0]/i}@gene_symbols;
	$ngrp1=@grp1;
	my @grp2=grep{/$C7[1]/i}@gene_symbols;
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
