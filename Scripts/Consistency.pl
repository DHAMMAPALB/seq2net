#!usr/bin/perl
#use strict;
#use warnings;
use List::Util qw(uniq);

$number=$ARGV[0];
$dir=$ARGV[1];

print "PPIs consistently predicted by atleast $number methods ...\n";
open(FILE1,"$dir/Gene_Neighborhood/GN_PPIs_gs.txt");my @A1=<FILE1>;close FILE1;$A1=join("",@A1);$A1=~s/\r//ig;my @B1=split(/\n/,$A1);
$legend=shift(@B1);@B1=sort(@B1);@B1=uniq(@B1);$n1=@B1;
open(FILE2,"$dir/Phylogenetic_Profile/PP_PPIs_gs.txt");my @A2=<FILE2>;close FILE2;$A2=join("",@A2);$A2=~s/\r//ig;my @B2=split(/\n/,$A2);
$legend=shift(@B2);@B2=sort(@B2);@B2=uniq(@B2);$n2=@B2;
open(FILE3,"$dir/Gene_Fusion/GF_PPIs_gs.txt");my @A3=<FILE3>;close FILE3;$A3=join("",@A3);$A3=~s/\r//ig;my @B3=split(/\n/,$A3);
$legend=shift(@B3);@B3=sort(@B3);@B3=uniq(@B3);$n3=@B3;
open(FILE4,"$dir/Gene_Co-evolution/GC_PPIs_gs.txt");my @A4=<FILE4>;close FILE4;$A4=join("",@A4);$A4=~s/\r//ig;my @B4=split(/\n/,$A4);
$legend=shift(@B4);@B4=sort(@B4);@B4=uniq(@B4);$n4=@B4;
open(FILE5,"$dir/Interlog/IN_PPIs_gs.txt");my @A5=<FILE5>;close FILE5;$A5=join("",@A5);$A5=~s/\r//ig;my @B5=split(/\n/,$A5);
$legend=shift(@B5);@B5=sort(@B5);@B5=uniq(@B5);$n5=@B5;

push(@ppis,uniq(@B1,@B2,@B3,@B4,@B5));
@ppis=sort(@ppis);
$nppis=@ppis;

open(OUT,">$dir/Consistent_PPIs_table.txt");
print OUT "Protein_1\tProtein_2\tGN\tPP\tGF\tGC\tIN\n";
$count1=0;$count2=0;$count3=0;$count4=0;$count5=0;
$fig=int($nppis/100);
$fig2=$fig;
for($i=0;$i<$nppis;$i++)
{
	if($fig2 eq $i){print ".";$fig2=$fig2+$fig;}
	my @grp1=grep{/$ppis[$i]/i}@B1;my @grp2=grep{/$ppis[$i]/i}@B2;my @grp3=grep{/$ppis[$i]/i}@B3;my @grp4=grep{/$ppis[$i]/i}@B4;my @grp5=grep{/$ppis[$i]/i}@B5;
	print OUT "$ppis[$i]\t";
	if($grp1[0] ne ""){ print OUT "1\t"; $count1++;} else { print OUT "0\t"; }
	if($grp2[0] ne ""){ print OUT "1\t"; $count2++;} else { print OUT "0\t"; }
	if($grp3[0] ne ""){ print OUT "1\t"; $count3++;} else { print OUT "0\t"; }
	if($grp4[0] ne ""){ print OUT "1\t"; $count4++;} else { print OUT "0\t"; }
	if($grp5[0] ne ""){ print OUT "1\n"; $count5++;} else { print OUT "0\n"; }
	$count=$count1+$count2+$count3+$count4+$count5;
	if($count>=$number)
	{
		push(@table,$ppis[$i]);
	}
	$count1=0;$count2=0;$count3=0;$count4=0;$count5=0;
}
close OUT;
$ntable=@table;

$out_f="$dir/Consistent_$number"."_PPIs_gs.txt";
open(OUT,">$out_f");
print OUT "Protein_1\tProtein_2\n";
for($i=0;$i<$ntable;$i++)
{
	print OUT "$table[$i]\n";
}
close OUT;
print "\nDone.\n";
