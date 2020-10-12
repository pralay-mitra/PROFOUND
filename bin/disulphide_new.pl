#!/usr/bin/perl

my $pdb=$a=$b=$c=$d=$e=$f=$g=$h=$sbridge=0;

my @t1;
my @t2;
my @t3;
my @t4;
my @t5;

$pdb=$ARGV[0];
#$resind=$ARGV[1];

@t1=`cat $pdb|grep 'ATOM'`;

open(FP1,'>',"temp_dsp.txt");
open(FP2,'>',"disulphide.txt");

while($a<=$#t1){
	chomp($t1[$a]);
	$aa=substr $t1[$a],17,3;
	$res_temp=$res=(substr $t1[$a],22,4)+0;
	$temp_x=$temp_y=$temp_z=$ncount=0;
	if($aa eq 'CYS'){
		while($res_temp==$res){
			$atom=substr $t1[$a],13,2;
			if($atom eq 'SG'){
				$x=substr $t1[$a],30,8;
				$y=substr $t1[$a],38,8;
				$z=substr $t1[$a],46,8;
				print FP1 "$res\tcys\t$x\t$y\t$z\n";
				}
			$a++;
			$res_temp=substr $t1[$a],22,4;
			}
		}	
	else{
		$a++;
		}
	}
undef @t1;
@t1=`cat temp_dsp.txt`;
for($a=0;$a<=$#t1;$a++){
	$t2[$a]=0;
	}
for($a=0;$a<$#t1;$a++){
	chomp($t1[$a]);
	@t3=split(/\s+/,$t1[$a]);
	chomp($t3[0]);
	chomp($t3[1]);
	chomp($t3[2]);
	chomp($t3[3]);
	chomp($t3[4]);
	#if($t3[0]!=$resind){
	#	next;
	#	}
	if($t2[$a]==0){
		for($b=0;$b<=$#t1;$b++){
			if($t2[$b]==0&&$t2[$a]==0){
				chomp($t1[$b]);
				@t4=split(/\s+/,$t1[$b]);
				chomp($t4[0]);
				chomp($t4[1]);
				chomp($t4[2]);
				chomp($t4[3]);
				chomp($t4[4]);
				if($t4[0]==$t3[0]){
					next;
					}
				$dist=sqrt(($t4[2]-$t3[2])**2+($t4[3]-$t3[3])**2+($t4[4]-$t3[4])**2);
				if(($dist<=4)&&($t3[1]eq$t4[1])){
					for($c=0;$c<=$#t1;$c++){
						if($t2[$c]==0){
							@t5=split(/\s+/,$t1[$c]);
							chomp($t5[0]);
							if($t5[0]==$t3[0]||$t5[0]==$t4[0]){
								#print FP2 "$c\n";
								$t2[$c]=1;
								}
							undef @t5;
							}
						}
					print FP2 "Disulphide bond between $t3[0] and $t4[0]\tdistance $dist\n";
					$sbridge++;
					}
				undef @t4;
				}
			}
		}
	undef @t3;
	}
print FP2 "Disulphide_bonds= $sbridge";
			
	
