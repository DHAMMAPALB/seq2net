#!/usr/bin/perl
#use strict;
#use warnings;
use List::Util qw(uniq);

$file=$ARGV[0];
$ref_genome=$ARGV[1];
$e_value=$ARGV[2];
$coverage=$ARGV[3];
$dir=$ARGV[4];

print "Orthologs prediction ...\n";
open(FILE,$file);
my @A=<FILE>;
close FILE;
$A=join("",@A);
$A=~s/\r//ig;
my @B=split(/\n/,$A);
@B=sort(@B);
$n=@B;

mkdir "$dir/Orthologs";
for($i=0;$i<$n;$i++)
{
	`sed 's/"//ig' $B[$i]`;

	my @in=split(/\//,$B[$i]);$nin=@in;
	$in_id="$dir/Orthologs/$in[$nin-1]";
	$in_id=~s/.faa$/.txt/ig;
	push(@ids,$in_id);
	open(FILE2,"$B[$i]");my @A2=<FILE2>;close FILE2;$A2=join("",@A2);$A2=~s/\r//ig;my @B2=split(/>/,$A2);$n2=@B2;
	open(OUT,"> $in_id"); for($j=1;$j<$n2;$j++){ my @D2=split(/\[/,$B2[$j]); print OUT "$D2[0]\n";	} close OUT;
	`makeblastdb -in $B[$i] -dbtype prot`;
}
for($i=0;$i<$n;$i++)
{
	if($B[$i] eq $ref_genome)
	{
		$ref_n=$i;last;
	}
}

for($i=0;$i<$n;$i++)
{
	if($i!=$ref_n)
	{
		my @out_q=split(/\//,$B[$ref_n]);$nout_q=@out_q;
		my @out_s=split(/\//,$B[$i]);$nout_s=@out_s;
		$out_f="$dir/Orthologs/$out_q[$nout_q-1]-vs-$out_s[$nout_s-1]";
		push(@fblastg,$out_f);
		`blastp -query $B[$ref_n] -subject $B[$i] -out $out_f -evalue $e_value -outfmt \"6 qseqid sseqid pident evalue bitscore\"`;
		$out_r="$dir/Orthologs/$out_s[$nout_s-1]-vs-$out_q[$nout_q-1]";
		push(@bblastg,$out_r);
		`blastp -query $B[$i] -subject $B[$ref_n] -out $out_r -evalue $e_value -outfmt \"6 qseqid sseqid pident evalue bitscore\"`;
		
		open(FILE3,"$out_f");my @A3=<FILE3>;close FILE3;$nA3=@A3;
		open(FILE4,"$out_r");my @A4=<FILE4>;close FILE4;$nA4=@A4;
		$out_rbbh="$dir/Orthologs/RBBH.txt";
		$ngrp=0;
		for($j=0;$j<$nA3;$j++)
		{
			my @B3=split(/\t/,$A3[$j]);my @grp=grep{/$B3[1]\t$B3[0]\t/ig}@A4;$ngrp=@grp;
			if($ngrp!=0)
			{
				my @C3=split(/\t/,$grp[0]); if($B3[2]>=$coverage && $C3[2]>=$coverage){ open(OUT,">>$out_rbbh");print OUT "$B3[0]\t$B3[1]\n";close OUT; } $ngrp=0;
			}
		}
	}
}

open(FILE5,"$dir/Orthologs/RBBH.txt");
my @A5=<FILE5>;
close FILE5;
$A5=join("",@A5);
$A5=~s/\r//ig;
my @B5=split(/\n/,$A5);
$n5=@B5;
@B5=sort(@B5);
for($i=0;$i<$n5;$i++)
{
	@ss=split(/\t/,$B5[$i]);
	push(@RV,$ss[0]);
	push(@ortho,$ss[1]);
}
$nn=@RV;
$out_rbbh_t="$dir/Orthologs/RBBH_table.txt";
open(OUT,">$out_rbbh_t");
print OUT $RV[0];
for($l=0,$k=1;$l<$nn;$l++,$k++)
{
	if($RV[$l] eq $RV[$k]) {print OUT "\t$ortho[$l]";}
	else { print OUT "\t$ortho[$l]\n$RV[$k]";}
}
close OUT;


open(FILE6,"$out_rbbh_t");my @A6=<FILE6>;close FILE6;$A6=join("",@A6);$A6=~s/\r//ig;my @B6=split(/\n/,$A6);$n6=@B6;
$nids=@ids;
$out_ortho_t="$dir/Orthologs/Orthologs_table.txt";
open(OUT,">$out_ortho_t");
my @in2=split(/\//,$ids[$ref_n]);$nin2=@in2;
print OUT "$in2[$nin2-1]";
for($i=0;$i<$nids;$i++)
{
	my @in2=split(/\//,$ids[$i]);$nin2=@in2;
	if($i!=$ref_n)
	{
		print OUT "\t$in2[$nin2-1]";
	}
}
print OUT "\n";

$fig=int($n6/100);
$fig2=$fig;
for($i=0;$i<$n6;$i++)
{
	if($fig2 eq $i){print ".";$fig2=$fig2+$fig;}
	my @C6=split(/\t/,$B6[$i]);
	print OUT "$C6[0]";
	$base=shift(@C6);
	@C6=uniq(@C6);
#	$nC6=@C6;
	for($j=0;$j<$nids;$j++)
	{
		if($j!=$ref_n)
		{
			open(FILE7,"$ids[$j]");
			my @A7=<FILE7>;
			close FILE7;
			$A7=join("",@A7);
			$A7=~s/\r//ig;
			my @B7=split(/\n/,$A7);
			$nB7=@B7;
			$nlist=0;
			for($k=0;$k<$nB7;$k++)
			{

				my @C7=split(/ /,$B7[$k]);
				my @grp=grep(/$C7[0]/,@C6);
				if($grp[0] ne "")
				{
					push(@list,@grp[0]);
				}
			}
			$nlist=@list;
			if($nlist==0)
			{
				print OUT "\t0";
			}
			else
			{
				print OUT "\t$list[0]";
				for($k=1;$k<$nlist;$k++)
				{
					print OUT ",$list[$k]"
				}
			}
			@list=();
		}
	}
	print OUT "\n";
}
close OUT;
print "\nDone.\n";
