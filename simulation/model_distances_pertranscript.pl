#!/usr/bin/env perl
use Modern::Perl;
use Math::Gauss ':all';

my %counts_5p;
my $ribo_footprint = 29;
my $cutdifficultyfactor = 100;
my $speed_factor = 20; #larger decreases height of 0 vs background
my $prob_of_capture = 0.01;
my $prob_to_cut = 0.6;
my $second_ribo_prob = 0.01;
my $second_ribo_distance = 15;

TRANSCRIPT: for (my $transcript = 0; $transcript < 100; $transcript++){

	warn $transcript."\n";
	my %loc_5p;
	my %loc_ribo;
	
	my $seq_length = 200 + int(rand(5000));
	
	my @prob_to_cut;
	for (my $i = 0; $i < $seq_length; $i++){
		$prob_to_cut[$i] = $prob_to_cut;
	}

	my @speed;
	for (my $i = 0; $i < $seq_length; $i++){
		$speed[$i] = 1;
	}
	
	for (my $i = 0; $i < $seq_length; $i+=3){
		$speed[$i+1] = 1.3*$speed[$i];
		$speed[$i+2] = 1.3*$speed[$i];
	}
	
	my $max_bumps = int($seq_length/200);
	for (my $i = 0; $i < $max_bumps; $i++){
		my $width_of_stall = 6;
		
		### get bump start so that midpoint of slowdown is in frame zero
		my $bump_start = $width_of_stall + int(rand($seq_length-(2*$width_of_stall)));
		my $bmid = $bump_start + ($width_of_stall/2);
		$bmid = $bmid - ($bmid % 3);
		$bump_start = $bmid - ($width_of_stall/2);
		
		for (my $pos = 0; $pos < $width_of_stall; $pos++){
			my $x = abs((($width_of_stall-1)/2)-$pos) * (3/($width_of_stall-1));
			my $slowdownfactor = 10 * pdf( $x, 0, 1 );
			my $abs_pos = $bump_start + $pos;
			
			if ($abs_pos % 3 == 0){
				$speed[$abs_pos] = $speed[$abs_pos] / $slowdownfactor;
			}
			if ($abs_pos % 3 != 0){
				$speed[$abs_pos] = $speed[$abs_pos] / ($slowdownfactor * 0.8);
			}
			if ($speed[$abs_pos] < 0.01){$speed[$abs_pos] = 0.01;}
		}
		for (my $pos = -1; $pos <= 1; $pos++){
			my $center_of_slowdown = $bump_start + int($width_of_stall/2);
			$prob_to_cut[$center_of_slowdown + $pos] = 0; #symmetric around center of slowdown
		}
	}
 
	my $transcript_level = 1 + int(rand(1000));
	ITERATION: for (my $perm = 0; $perm < $transcript_level; $perm++){
		for (my $i = 0; $i < $seq_length-$ribo_footprint; $i++){
			my $local_speed = $speed[$i+$ribo_footprint]; ### counting speed at the front of the ribosome
			my $time_on_nt = int($speed_factor/$local_speed);
			
			for (my $time = 0; $time < $time_on_nt; $time++){
				if (rand() < $prob_of_capture){
					$loc_ribo{$i}++;
				}
				if (rand() < $prob_to_cut[$i]){
					$loc_5p{$i}++;
					if (rand() < $second_ribo_prob){ #second ribosome stalls
						my $pos_second_ribo = int($i-$second_ribo_distance - rand(5) + rand(5));
						$loc_ribo{$pos_second_ribo}++;
					}
				}
			}
		}
	}

	foreach my $pos1 (keys %loc_5p){
		foreach my $pos2 (keys %loc_ribo){
			if($pos2 < 80 or $pos2 > ($seq_length-80)){next;}
			my $dist = $pos1 - $pos2;
			if (abs($dist) > 40){next;}
			else {
				$counts_5p{$dist} += $loc_5p{$pos1} * $loc_ribo{$pos2};
			}
		}
	}
}
print join("\t", ("from","distance","count"))."\n";

for (my $dist = -40; $dist <= 40; $dist++){
	if (!exists $counts_5p{$dist}){$counts_5p{$dist} = 0;}
	print join("\t", '5p', $dist, $counts_5p{$dist})."\n";
}

exit;
