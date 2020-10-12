#!/usr/bin/perl
#For new usage. Provide starting residue id and chain id to renumber and rename accordingly.

my $id=$ind=$a=$b=$bb=$c=$cc=$d=$e=$f=$chain=0;

my @t1;
my @t2;
my @t3;

$id=$ARGV[0];
@t1=`cat $id.pdb`;
open(FP,'>',"$id.rect");

chomp($t1[0]);
$cc=substr $t1[0],0,21;
$ind=$ARGV[1];
$chain=$ARGV[2];
$b=$cc.$chain;
$d=substr $t1[0],27;
$ind++;
$ind--;
if($ind>=0&&$ind<=9){
	print FP "$b   $ind $d\n";
	}elsif(($ind>=10&&$ind<=99)||($ind<0&&$ind>-10)){
	print FP "$b  $ind $d\n";
	}elsif($ind>=1000||($ind<-99&&$ind>-1000)){
	print FP "$b$ind $d\n";
	}else{
	print FP "$b $ind $d\n";
	}
for($a=1;$a<=$#t1;$a++){
	chomp($t1[$a]);
	$bb=substr $t1[$a],0,3;
	if($bb eq 'TER'){
		print FP "TER";
		}else{
		$cc=substr $t1[$a],0,21;
		$b=$cc.$chain;
		$c=substr $t1[$a],23,4;
		$d=substr $t1[$a],27;
		$e=substr $t1[$a-1],23,4;
		if($c ne $e){
			$ind++;
			if($ind>=0&&$ind<=9){
				print FP "$b   $ind $d\n";
				}elsif(($ind>=10&&$ind<=99)||($ind<0&&$ind>-10)){
				print FP "$b  $ind $d\n";
				}elsif($ind>=1000||($ind<-99&&$ind>-1000)){
				print FP "$b$ind $d\n";
				}else{
				print FP "$b $ind $d\n";
				}
			}else{
			if($ind>=0&&$ind<=9){
				print FP "$b   $ind $d\n";
				}elsif(($ind>=10&&$ind<=99)||($ind<0&&$ind>-10)){
				print FP "$b  $ind $d\n";
				}elsif($ind>=1000||($ind<-99&&$ind>-1000)){
				print FP "$b$ind $d\n";
				}else{
				print FP "$b $ind $d\n";
				}
			}
		}
	}	
