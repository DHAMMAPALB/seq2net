#!usr/bin/perl
#use strict;
#use warnings;

$file=$ARGV[0];
$ref_genome=$ARGV[1];
$dir=$ARGV[2];

print "Gene Co-evolution ...\n";
open(FILE,$file);
my @A=<FILE>;
close FILE;
$A=join("",@A);
$A=~s/\r//ig;
my @files=split(/\n/,$A);
@files=sort(@files);

$o_table="$dir/OrthoParalogs/Orthologs_table.txt";
open(FILE,$o_table);
my @A=<FILE>;
close FILE;
$A=join("",@A);
my @B=split(/\n/,$A);
$n=@B;
my @names=split(/\t/,$B[0]);
$n_names=@names;
for($i=0;$i<$n_names;$i++)
{
	print "$i. $names[$i]\n";
	$names[$i]=~s/.txt$/.faa/ig;
	my @pids;
	for($j=1;$j<$n;$j++)
	{
		my @C=split(/\t/,$B[$j]);
		my @D=split(/\,/,$C[$i]);
		if($D[0] ne 0)
		{
			push(@pids,$D[0]);
		}
	}
	$npids=@pids;
	
	my @grpf=grep{/$names[$i]$/ig}@files;
	open(FILE,$grpf[0]);
	my @AA=<FILE>;
	close FILE;
	$AA=join("",@AA);
	my @BB=split(/\>/,$AA);
	$nn=@BB;
	my @pids2;
	for($j=0;$j<$nn;$j++)
	{
		my @CC=split(/ /,$BB[$j]);
		push(@pids2,$CC[0]);
	}
	$npids2=@pids2;
	
	for($j=0;$j<$npids;$j++)
	{
		for($k=0;$k<$npids2;$k++)
		{
			if($pids[$j] eq $pids2[$k])
			{
				push(@clust_input,$BB[$k]);
			}
		}
	}
}
my @BB=@clust_input;
mkdir "$dir/Gene_Co-evolution";
mkdir "$dir/Gene_Co-evolution/Mirror_Tree";
$fig=int($n/100);
$fig2=$fig;
for($i=1;$i<$n;$i++)
{
	if($fig2 eq $i){print ".";$fig2=$fig2+$fig;}
	my @C=split(/\t/,$B[$i]);
	for($j=$i+1;$j<$n;$j++)
	{
		my @D=split(/\t/,$B[$j]);
		$nD=@D;
		my @one; my @two; my @names2;	
		for($k=0;$k<$nD;$k++)
		{
			if($D[$k] ne 0 && $C[$k] ne 0)
			{
				my @D1=split(/\,/,$C[$k]);
				my @D2=split(/\,/,$D[$k]);
				push(@one,$D1[0]);
				push(@two,$D2[0]);
			}
		}
		$n_one=@one;
		$n_two=@two;
		if($n_one>=5)
		{
			print "$i.$j \t$n_one\t$n_two\n";
			$directory="$one[0]-"."$two[0]";
			mkdir "$dir/Gene_Co-evolution/Mirror_Tree/$C[0]/";
			mkdir "$dir/Gene_Co-evolution/Mirror_Tree/$C[0]/$directory";
			open(OUT,">$dir/Gene_Co-evolution/Mirror_Tree/$C[0]/$directory/$one[0].txt");
			for($m=0;$m<$n_one;$m++)
			{
				my @grp1=grep{/$one[$m]/ig}@BB;
				$grp1[0]=~s/\>//ig;
				print OUT ">$grp1[0]";
			}
			close OUT;
			open(OUT,">$dir/Gene_Co-evolution/Mirror_Tree/$C[0]/$directory/$two[0].txt");
			for($m=0;$m<$n_two;$m++)
			{
				my @grp2=grep{/$two[$m]/ig}@BB;
				$grp2[0]=~s/\>//ig;
				print OUT ">$grp2[0]";
			}
			close OUT;
			
			$dist_mat1="$one[0]".".mat";
			$dist_dnd1="$one[0]".".dnd";
			system("clustalo -i $dir/Gene_Co-evolution/Mirror_Tree/$C[0]/$directory/$one[0].txt --dealign --full --distmat-out=$dir/Gene_Co-evolution/Mirror_Tree/$C[0]/$directory/$dist_mat1 --guidetree-out=$dir/Gene_Co-evolution/Mirror_Tree/$C[0]/$directory/$dist_dnd1 --force");
			$dist_mat2="$two[0]".".mat";
			$dist_dnd2="$two[0]".".dnd";
			system("clustalo -i $dir/Gene_Co-evolution/Mirror_Tree/$C[0]/$directory/$two[0].txt --dealign --full --distmat-out=$dir/Gene_Co-evolution/Mirror_Tree/$C[0]/$directory/$dist_mat2 --guidetree-out=$dir/Gene_Co-evolution/Mirror_Tree/$C[0]/$directory/$dist_dnd2 --force");

			open(FILE,"$dir/Gene_Co-evolution/Mirror_Tree/$C[0]/$directory/$dist_mat1");
			my @AAAA=<FILE>;
			close FILE;
			$AAAA=join("",@AAAA);
			$AAAA=~s/ /\t/ig;
			my @BBBB=split(/\n/,$AAAA);
			$nnnn=@BBBB;
			$mat1="$dist_mat1".".txt";
			open(OUT,">$dir/Gene_Co-evolution/Mirror_Tree/$C[0]/$directory/$mat1");
			for($ii=1;$ii<$nnnn;$ii++)
			{
				my @CCCC=split(/\t/,$BBBB[$ii]);
				$nnnnn=@CCCC;
				for($jj=0;$jj<$nnnnn;$jj++)
				{
					if($CCCC[$jj] ne "")
					{
					print OUT "$CCCC[$jj]\t";
					}
				}
				print OUT "\n";
			}
			close OUT;

			open(FILE,"$dir/Gene_Co-evolution/Mirror_Tree/$C[0]/$directory/$dist_mat2");
			my @AAAA=<FILE>;
			close FILE;
			$AAAA=join("",@AAAA);
			$AAAA=~s/ /\t/ig;
			my @BBBB=split(/\n/,$AAAA);
			$nnnn=@BBBB;
			$mat2="$dist_mat2".".txt";	
			open(OUT,">$dir/Gene_Co-evolution/Mirror_Tree/$C[0]/$directory/$mat2");
			for($ii=1;$ii<$nnnn;$ii++)
			{
				my @CCCC=split(/\t/,$BBBB[$ii]);
				$nnnnn=@CCCC;
					for($jj=0;$jj<$nnnnn;$jj++)
				{
					if($CCCC[$jj] ne "")
					{
						print OUT "$CCCC[$jj]\t";
					}
				}
				print OUT "\n";
			}
			close OUT;
			
			open(FILE1,"$dir/Gene_Co-evolution/Mirror_Tree/$C[0]/$directory/$mat1");
			my @A1=<FILE1>;
			close FILE1;
			$A1=join("",@A1);
			my @B1=split(/\n/,$A1);
			$n1=@B1;

			$x_sum=0;my @xi;
			for($ii=0;$ii<$n1;$ii++)
			{
				my @C1=split(/\t/,$B1[$ii]);
				$nC1=@C1;
				for($jj=$ii+2;$jj<$nC1;$jj++)
				{
					$x_sum=$x_sum+$C1[$jj];
					push(@xi,$C1[$jj]);
				}
			}
			$nxi=@xi;
		
			$x_mean=$x_sum/$nxi;

			open(FILE2,"$dir/Gene_Co-evolution/Mirror_Tree/$C[0]/$directory/$mat2");
			my @A2=<FILE2>;
			close FILE2;
			$A2=join("",@A2);
			my @B2=split(/\n/,$A2);
			$n2=@B2;

			$y_sum=0; my @yi;
			for($ii=0;$ii<$n2;$ii++)
			{
				my @C2=split(/\t/,$B2[$ii]);
				$nC2=@C2;
				for($jj=$ii+2;$jj<$nC2;$jj++)
				{
					$y_sum=$y_sum+$C2[$jj];
					push(@yi,$C2[$jj]);
				}
			}
			$nyi=@yi;
			$y_mean=$y_sum/$nyi;

			$cov_sum=0;
			for($ii=0;$ii<$nxi;$ii++)
			{
				$cov=($xi[$ii]-$x_mean)*($yi[$ii]-$y_mean);
				$cov_sum=$cov_sum+$cov;
			}
			
			$var_x_sum=0;
			for($ii=0;$ii<$nxi;$ii++)
			{
				$var_x=($xi[$ii]-$x_mean)*($xi[$ii]-$x_mean);
				$var_x_sum=$var_x_sum+$var_x;
			}
			$std_x=sqrt($var_x_sum);
			
			$var_y_sum=0;
			for($ii=0;$ii<$nyi;$ii++)
			{
				$var_y=($yi[$ii]-$y_mean)*($yi[$ii]-$y_mean);
				$var_y_sum=$var_y_sum+$var_y;
			}
			$std_y=sqrt($var_y_sum);

			$correlation=$cov_sum/($std_x*$std_y);
			$line="$one[0]\t$two[0]";
			if($correlation>=0.8)
			{
				open(OUT,">>$dir/Gene_Co-evolution//GC_PPIs.txt");
				print OUT "$line\n";
				push(@ppis,$line);
			}

			$ddd="$one[0].txt";
			open(OUT,">>$dir/Gene_Co-evolution/Mirror_Tree/$C[0]/$ddd");
			print OUT "$line\t$correlation\n";
			close OUT;
						
		}
	}
}
$nppis=@ppis;

$ref_genome_t=$ref_genome;
$ref_genome_t=~s/.faa$/.csv/ig;
open(FILE3,$ref_genome_t);
my @A3=<FILE3>;
close FILE3;
$A3=join("",@A3);
$A3=~s/\r//ig;
$A3=~s/\"//ig;
my @B3=split(/\n/,$A3);
$n3=@B3;
for($i=0;$i<$n3;$i++)
{
	my @C3=split(/,/,$B3[$i]);
	if($C3[6] ne "")
	{
		push(@gene_symbols,"$C3[8]\t$C3[6]");
	}
	else
	{
		push(@gene_symbols,"$C3[8]\t$C3[7]");
	}
}
open(OUT,">$dir/Gene_Co-evolution/GC_PPIs_gs.txt");
$ngrp1=0;$ngrp2=0;
print OUT "Protein_1\tProtein_2\n";
for($i=0;$i<$nppis;$i++)
{
	my @C4=split(/\t/,$ppis[$i]);
	my @grp1=grep{/$C4[0]/i}@gene_symbols;
	$ngrp1=@grp1;
	my @grp2=grep{/$C4[1]/i}@gene_symbols;
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
