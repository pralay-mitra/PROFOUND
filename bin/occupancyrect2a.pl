#!/usr/bin/perl

my $id=$a=$b=$c=$c1=$d=$d1=$f=$f1=$g=$h=$i=$j=$k=$k1=$l=$concat=0;
my @t1;
my @t2;
my @t3;
my @t4;
my @t5;
my @t6;
my @t7;

$id=$ARGV[0];#pdb id to be rectified

@t1=`cat $id.pdb`;
open(FP,'>',"$id.rect");

@t6=`cat $id.pdb|grep "CA"`;
$h=0;
for($a=0;$a<$#t6;$a++){
	chomp($t6[$a]);
	chomp($t6[$a+1]);
	$b=substr $t6[$a],17,3;
	$c=substr $t6[$a],22,4;
	$d=substr $t6[$a],56,4;
	$e=substr $t6[$a+1],17,3;
	$f=substr $t6[$a+1],22,4;
	$g=substr $t6[$a+1],56,4;
	if(($f==$c)&&($e ne $b)){
		if($d>=$g){
			$t7[$h][0]=$f;
			$t7[$h][1]=$e;
			}else{
			$t7[$h][0]=$f;
			$t7[$h][1]=$b;
			}
		$h++;
		}
	}
for($a=0;$a<=$#t1;$a++){
	$t2[$a]=1;
	chomp($t1[$a]);
	$b=substr $t1[$a],17,3;
	$c=substr $t1[$a],22,4;
	for($d=0;$d<=$#t7;$d++){
		if(($t7[$d][0]==$c)&&($t7[$d][1]eq$b)){
			$t2[$a]=0;
			last;
			}
		}
	}

for($a=0;$a<=$#t1;$a++){
	chomp($t1[$a]);
	$b=substr $t1[$a],0,16;
	$k=substr $t1[$a],0,3;
	$c=substr $t1[$a],13,3;
	$d=substr $t1[$a],22,4;
	$e=substr $t1[$a],17,39;
	$f=substr $t1[$a],56,4;
	$g=substr $t1[$a],60;
	if($k eq 'TER'){
		print FP "$t1[$a]\n";
		next;
		}
	if($f==1.00){
		print FP "$b $e$f$g\n";
		}else{
		if($t2[$a]==0){
			next;
			}else{
			$d1=$d;
			$h=1;
			$i=$a+$h;
			while(($d1==$d)&&($i<=$#t1)){
				$i=$a+$h;
				if($t2[$i]==0){
					$h++;
					$i++;
					$d1=substr $t1[$i],22,4;
					next;
					}
				$c1=substr $t1[$i],13,3;
				$f1=substr $t1[$i],56,4;
				if(($c1 eq $c)&&($f1<=$f)){
					$t2[$i]=0;
					}
				if(($c1 eq $c)&&($f1>$f)){
					$t2[$a]=0;
					}
				$h++;
				$i++;
				$d1=substr $t1[$i],22,4;
				}
			if($t2[$a]==1){
				$concat=$b.' '.$e.'1.00'.$g;
				print FP "$concat\n";
				}
			}
		}
	}
undef(@t1);
close(FP);
@t1=`cat $id.rect`;
open(FP,'>',"$id.final");

for($a=0;$a<=$#t1;$a++){
	chomp($t1[$a]);
	$h=substr $t1[$a],0,3;
	$b=substr $t1[$a],11;
	$c=$a+1;
	if($c<10){
		if($a==$#t1){
			if($h eq 'TER'){
				print FP "TER       $c$b\n";
				}else{
				print FP "ATOM      $c$b\n";
				}
			last;
			}
		print FP "ATOM      $c$b\n";
		}
	if($c>9&&$c<100){
		if($a==$#t1){
			if($h eq 'TER'){
				print FP "TER      $c$b\n";
				}else{
				print FP "ATOM     $c$b\n";
				}
			last;
			}
		print FP "ATOM     $c$b\n";
		}
	if($c>99&&$c<1000){
		if($a==$#t1){
			if($h eq 'TER'){
				print FP "TER     $c$b\n";
				}else{
				print FP "ATOM    $c$b\n";
				}
			last;
			}
		print FP "ATOM    $c$b\n";
		}
	if($c>999&&$c<10000){
		if($a==$#t1){
			if($h eq 'TER'){
				print FP "TER    $c$b\n";
				}else{
				print FP "ATOM   $c$b\n";
				}
			last;
			}
		print FP "ATOM   $c$b\n";
		}
	if($c>9999&&$c<100000){
		if($a==$#t1){
			if($h eq 'TER'){
				print FP "TER   $c$b\n";
				}else{
				print FP "ATOM  $c$b\n";
				}
			last;
			}
		print FP "ATOM  $c$b\n";
		}
	}
`mv $id.final $id.rect`;
#`perl resnumchange2.pl $id`;	
#`perl /home/Anupam/PhD/library/occupancy_rect/resnumchange2.pl $id`;	
