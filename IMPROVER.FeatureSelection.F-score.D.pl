#!usr/bin/perl -w 
use strict;
use warnings;

my $line;
my %symbol = ();
my %SampleLable = ();

die "perl $0 <Annotation File><ExpressionMatrix><Group File>" unless $#ARGV== 2;

&main();

my $cmd = "perl IMPROVER.TestExtractor.pl $ARGV[0] GDS1209.NormalizedExpression.xls $ARGV[1].$0.Features.csv GDS1209.Groups.3.table";
$cmd = `$cmd`;

$cmd = "perl IMPROVER.TestExtractor.pl $ARGV[0] GDS2736.NormalizedExpression.xls $ARGV[1].$0.Features.csv GDS2736.Groups.3.table";
$cmd = `$cmd`;

sub main () {
	my $annotation = &load_annotation($ARGV[0]);
	my $expression_matrix = &load_Expression($ARGV[1]);
	my $groups = &load_groups($ARGV[2]);
	my $Zscore = &getZscore($expression_matrix);
	my $features = &FeaturesProcess($Zscore,$annotation,$groups,$expression_matrix);
	
	open (OUT,">".$ARGV[1].".".$0.".Train.csv") or die "$!";
	print OUT ("Label,Class,",join(",", sort {$a cmp $b} keys %{$features}),"\n");
	foreach my $iclass (grep {$SampleLable{$_} ne "Test"} keys %SampleLable) {
		print OUT ($iclass,",",$SampleLable{$iclass});
		foreach my $ifeature (sort {$a cmp $b} keys %{$features}) {
			my $matrix = &getPAC($features->{$ifeature},$Zscore);
			print OUT (",",$matrix->{$iclass});
		}
		print OUT ("\n");
	}
	close OUT or die "$!";
	
	open (OUT,">".$ARGV[1].".".$0.".Features.csv") or die "$!";
	foreach my $ifeature (sort {$a cmp $b} keys %{$features}) {
		print OUT ($ifeature,",",join(",",@{$features->{$ifeature}}),"\n");
	}
	close OUT or die "$!";
}


sub FeaturesProcess () {
	my $Zscore = shift;
	my $annotation = shift;
	my $groups = shift;
	my $expression_matrix = shift;
	my %results = ();my %resultsDS = ();
	my @achor = grep {$_ ne "Test"} keys %$groups;
	for (my $i = 0; $i <= $#achor; $i++) {
		my %features = ();
		my %features_DS = ();
		warn "Selecting :".$achor[$i]." specific probes.";
		my @class_1 = keys %{$groups->{$achor[$i]}};
		my @class_2 = grep {$SampleLable{$_} ne "Test"} keys %SampleLable;
		foreach my $domain (keys %{$annotation}) {
			#warn $domain;
			my @subset = keys %{$annotation->{$domain}};
			@subset = grep {exists $expression_matrix->{$_}} @subset;
			my $tmp = &pickCFS(\@subset,$Zscore,\@class_1,\@class_2,$domain);
			my $posiMarker = 0;
			$posiMarker = abs(&getDS(&getPAC($tmp -> {'Positive'},$Zscore),\@class_1,\@class_2)) unless (scalar(@{$tmp->{'Positive'}}) == 0);
			my $negaMarker = 0;
			$negaMarker = abs(&getDS(&getPAC($tmp -> {'Negative'},$Zscore),\@class_1,\@class_2)) unless (scalar(@{$tmp->{'Negative'}}) == 0);
			do {
				if (not (exists $features{$domain."_Positive"})) {
					@{$features{$domain."_Positive"}} = @{$tmp -> {'Positive'}};
					$features_DS{$domain."_Positive"} = $posiMarker;
				}else{
					warn "CrossValidation::".$domain."_Positive"." is already exists.\n";
				}
			} unless (scalar(@{$tmp->{'Positive'}}) == 0 || $posiMarker < 1);
			do {
				if (not (exists $features{$domain."_Negative"})) {
					@{$features{$domain."_Negative"}} = @{$tmp -> {'Negative'}};
					$features_DS{$domain."_Negative"} = $negaMarker;
				}else{
					warn "CrossValidation::".$domain."_Negative"." is already exists.\n";
				}
			} unless (scalar(@{$tmp->{'Negative'}}) == 0 || $negaMarker < 1);
		}
		#die "FeaturesProcess::PickTop5 Failed.\n" unless (&pickTop5(\%output,\%features_DS,$achor[$i]));
		if (not (exists $results{$achor[$i]})) {
			foreach my $ifeature (keys %features) {
				@{$results{$achor[$i]}{$ifeature}} = @{$features{$ifeature}};
				$resultsDS{$achor[$i]}{$ifeature} = $features_DS{$ifeature};
			}
		}else{
			die "FeaturesProcess:: ".$achor[$i]." is already exists in the HASH Features.\n";
		}
	}
	my %output = ();
	die "FeaturesProcess:: FeatureSelection Failed." unless (&FeatureSelection(\%results,\%resultsDS,\%output));
	#die "FeaturesProcess::PickTop5 Failed.\n" unless (&pickTop5(\%output,\%features_DS,$achor[$i]));
	return \%output;
}


sub FeatureSelection () {
	my $features = shift;
	my $featuresDS = shift;
	my $neo_Features = shift;
	my %median = ();my %temp = ();my %libsH = ();my @libs = ();
	warn "FeatureSelection::Preformed.";
	foreach my $ilab (keys %{$features}) {
		foreach my $ifeature (keys %{$features->{$ilab}}) {
			if (not (exists $median{$ifeature}{$ilab})) {
				@{$median{$ifeature}{$ilab}} = @{$features->{$ilab}{$ifeature}} unless (scalar(@{$features->{$ilab}{$ifeature}}) <= 0);
				#warn "FeatureSelection::".$ifeature." ".$ilab." inject ".scalar(@{$features->{$ilab}{$ifeature}})." probes.";
			}else{
				warn "FeatureSelection::".$ifeature." ".$ilab." is already exists in the HASH median.";
			}
		}
	}
	foreach my $ifeature (keys %median) {
		@libs = keys %{$median{$ifeature}};
		my %deletes = ();
		for (my $i = 0; $i < $#libs; $i++) {
			for (my $j = $i + 1; $j <= $#libs; $j++) {
				if (&ArrayCompare($median{$ifeature}{$libs[$i]},$median{$ifeature}{$libs[$j]})) {
					if ($featuresDS->{$libs[$i]}{$ifeature} > $featuresDS->{$libs[$j]}{$ifeature} ) {
						$deletes{$libs[$j]} = 1;
					}
				}
			}
		}
		@libs = grep {not exists $deletes{$_}} @libs;
		foreach my $ilib (@libs) {
			@{$temp{$ilib}{$ifeature}} = @{$median{$ifeature}{$ilib}};
		}
		@libs = keys %{$features};
	}
	map {die "FeatureSelection::PickTop5 Failed.\n" unless (&pickTop5(\%temp,$featuresDS,$_))} @libs;
	foreach my $ilib (@libs) {
		foreach my $ifeature (keys %{$temp{$ilib}}) {
			if (not (exists $neo_Features -> {$ifeature."_".$ilib})) {
				@{$neo_Features -> {$ifeature."_".$ilib}} = @{$temp{$ilib}{$ifeature}} ;
			}else{
				die "ERRO::FeatureSelection::Redundency Found";
			}
		}
	}
	return 1;
}

sub pickTop5 () {
	my $features = shift;
	my $DS = shift;
	my $type = shift;
	my @deletes = ();
	my @totalFeatures = keys %{$features->{$type}};
	#########deredundcy######
	for (my $i = 0; $i < $#totalFeatures; $i++) {
		for (my $j = $i+1; $j <= $#totalFeatures; $j++) {
			die "empty array 1" unless (scalar(@{$features->{$type}{$totalFeatures[$i]}}) > 0);
			die "empty array 2" unless (scalar(@{$features->{$type}{$totalFeatures[$j]}}) > 0);
			if (&ArrayCompare($features->{$type}{$totalFeatures[$i]},$features->{$type}{$totalFeatures[$j]})) {
				push @deletes, $totalFeatures[$i];
			}
		}
	}
	#map {delete($features->{$type}{$_});warn"Warning::pickTop5:Ã»ÓÐÉ¾¸É¾»:$type $_"} @deletes;
	@totalFeatures = keys %{$features->{$type}};
	####################ranking###########
	@totalFeatures = keys %{$features->{$type}};
	for (my $i = 0; $i < $#totalFeatures; $i++) {
		for (my $j = $i+1; $j <= $#totalFeatures; $j++) {
			if ($DS->{$type}{$totalFeatures[$i]} < $DS->{$type}{$totalFeatures[$j]}) {
				my $temp = $totalFeatures[$i];
				$totalFeatures[$i] = $totalFeatures[$j];
				$totalFeatures[$j] = $temp;
			}
		}
	}
	##################
	@deletes = ();
	if (scalar(@totalFeatures) > 5) {
		@deletes = @totalFeatures[5..$#totalFeatures];
	}else{
		warn "<= 5 Feature Found.Picking Top 5 Features Proccess Skipped.\n";
	}
	open (OUT,">".$0.".".$ARGV[2].".".$type.".xls") or die "$!";
	foreach my $i (@totalFeatures) {
		print OUT ($i,"\t",$DS->{$type}{$i},"\n");
	}
	close OUT or die "$!";
	map {
		delete($features->{$type}{$_})
		} @deletes;
	return 1;
}

sub ArrayCompare () {
	my $array1 = shift;
	my $array2 = shift;
	return 0 unless (scalar(@{$array1}) == scalar(@{$array2}));
	my %tmp = ();
	map {$tmp{$_} = 1} @{$array2};
	foreach my $i (@{$array1}) {
		if (not (exists $tmp{$i})) {
			return 0;
		}
	}
	return 1;
}

sub pickCFS () {
	#warn "Picking Correlated Subset...\n";
	my $addre = shift;
	my $Zscore = shift;
	my $class_1 = shift;
	my $class_2 = shift;
	my $tag = shift;
	my %out = ();
	#warn "pickCFS:: Class1 =".scalar(@{$class_1})." Class2 =".scalar(@{$class_2});
	#warn "pickCFS:: candidate probes: ".scalar(@{$addre});
	my $f_Matrix = &getFscoreMatrix($addre,$Zscore,$class_1,$class_2);
	my @tmp = @{$addre}; ##### Gene List
	do {
		warn "No gene included,skipped.\n";
	} unless (scalar(@tmp) >= 1);
	my @temp = ();
	my %positiveData = ();
	my %negativeData = ();
	map {
		if ($f_Matrix->{$_} >= 0) {
			$positiveData{$_} = $f_Matrix->{$_};
		}else{
			$negativeData{$_} = $f_Matrix->{$_};
		}
	} @tmp;
	my @ranked_P_Tscore_index = @{&rank(\%positiveData,1)};
	my @ranked_N_Tscore_index = @{&rank(\%negativeData,2)};
	my $index_max = &larger(scalar(@ranked_P_Tscore_index),scalar(@ranked_N_Tscore_index));
	#do {
	#	#warn "ranked_P_Tscore_index ".scalar(@ranked_P_Tscore_index)." != ". scalar(@ranked_N_Tscore_index)." ranked_N_Tscore_index.\nProcess...";
	#	if (scalar(@ranked_P_Tscore_index) > scalar(@ranked_N_Tscore_index)) {
	#		@ranked_N_Tscore_index[scalar(@ranked_N_Tscore_index)..$#ranked_P_Tscore_index] = (-1)x($#ranked_P_Tscore_index-scalar(@ranked_N_Tscore_index) + 1);
	#	}else{
	#		@ranked_P_Tscore_index[scalar(@ranked_P_Tscore_index)..$#ranked_N_Tscore_index] = (-1)x($#ranked_N_Tscore_index-scalar(@ranked_P_Tscore_index) + 1);
	#	}
	#	die "Failed.\n" unless (scalar(@ranked_P_Tscore_index) == scalar(@ranked_N_Tscore_index));
	#	#warn "Filling Ok.\n";
	#} unless (scalar(@ranked_P_Tscore_index) == scalar(@ranked_N_Tscore_index));
	my @posiSet = ();my @negaSet = ();
	my $positive_marker = 0;
	my $negative_marker = 0;
	if (scalar(@ranked_P_Tscore_index) > 2) {
		@posiSet = ($ranked_P_Tscore_index[0]) ;
		#warn "////////Positive Analyzing\n";
		die "Picking Positive Data Failed.\n" unless (&selectHits(\@posiSet,$Zscore,\@ranked_P_Tscore_index,$class_1,$class_2,\$positive_marker,\%positiveData));
	}
	if (scalar(@ranked_N_Tscore_index) > 2 ) {
		@negaSet = ($ranked_N_Tscore_index[0]);
		#warn "////////Negative Analyzing\n";
		die "Picking Negative Data Failed.\n" unless (&selectHits(\@negaSet,$Zscore,\@ranked_N_Tscore_index,$class_1,$class_2,\$negative_marker,\%negativeData));
	}
	warn "Positive: ".scalar(@posiSet)." Negative: ".scalar(@negaSet)." of ".scalar(@tmp)." Subset Selected.\n";
	@{$out{"Positive"}} = @posiSet;
	@{$out{"Negative"}} = @negaSet;
	return \%out;
}

sub selectHits () {
	my $selectSet = shift;
#	my $selectData = shift;
	my $Zscore = shift;
	my $selectIndex = shift;
	my $class_1 = shift;
	my $class_2 = shift;
	my @classA = (@{$class_1},@{$class_2});
	my $marker = shift;
	my $index_max = scalar(@{$selectIndex});
	my $fscore = shift;
	#my $ifPositive = shift;
	my $PAC_mark = abs(&getDS(&getPAC($selectSet,$Zscore),$class_1,$class_2));
	#warn "initial DS value $PAC_mark.\n";
	#warn "Positive TScore: ".$selectData->{$selectIndex->[0]}."\n";
	my $PAC_current = 99999999;
	my $iter = 0;
	while (1) {
		#warn "selectHits:: iter: ".$iter;
		$iter ++;
		last unless (scalar(@{$selectSet}) < scalar(@{$selectIndex}));
		last unless ($iter <= scalar(@{$selectIndex}));
		#warn $iter." ".join("\t",@{$selectSet});
		my %mCORGsScore = ();
		my @remain = @{$selectIndex};
		foreach my $key (@{$selectSet}) {
			@remain = grep {($_ ne $key) && ($symbol{$_} ne $symbol{$key})} @remain;
		}
		last unless (scalar(@remain) > 0);
		foreach my $ilab (@remain) {
			###$selectSet,$genei,$Zscore,$fscore,$class
			$mCORGsScore{$ilab} = &getmCORGsScoreD($selectSet,$ilab,$Zscore,$fscore,\@classA);  #######################################Different#############################
		}
		my @RankedmCORGsScore = @{&rank(\%mCORGsScore,1)};
		#warn "Remains:".scalar(@remain)." ".scalar(@RankedmCORGsScore);
		#warn $RankedmCORGsScore[0]."[".$mCORGsScore{$RankedmCORGsScore[0]}."] ".$RankedmCORGsScore[1]."[".$mCORGsScore{$RankedmCORGsScore[1]};
		#warn "Going to iterating...$iter\n";
		#my $flag = 0;
		#map {if ($symbol{$_} eq $symbol{$RankedmCORGsScore[0]}) {$flag =1}} @{$selectSet};
		#if ($flag != 1) {
		#	push @{$selectSet},$RankedmCORGsScore[0];
		#}else{
		#	next;
		#	#warn "selectHits:: redundency Probe ".$RankedmCORGsScore[0]."\n";
		#}
		#warn "Select TScore: ".$selectData->{$selectIndex->[$iter]}."\n";
		push @{$selectSet},$RankedmCORGsScore[0];
		$PAC_current = &getDS(&getPAC($selectSet,$Zscore),$class_1,$class_2);
		#warn "PAC_current: $PAC_current <??> PAC_mark: $PAC_mark .\n";
		if (abs($PAC_current) >= abs($PAC_mark)) {
			$PAC_mark = abs($PAC_current);
		}else{
			#warn "Current PAC_Marker small than last one out of itertion.";
			if (scalar(@{$selectSet}) > 1) {
				pop(@{$selectSet});
			}
			last;
		}
		#sleep(2);
	}
	$$marker = $PAC_mark;
	return 1;
}

sub getPAC () {
	my $subset = shift;
	my $Zscore = shift;
	my %PACmatrix = ();
	do {
		warn "getPAC:: empty subset.\n";
		return \%PACmatrix} unless (scalar(@$subset) > 0);
	my @samples = keys %{$Zscore->{$subset->[0]}};
	for (my $i = 0; $i <= $#samples; $i++) {
		#warn "getPACi::i:$i.\n";
		#warn "pacCalc: ".&pacCalc($plus,$i)." ".&pacCalc($minu,$i)."\n";
		if (not (exists $PACmatrix{$samples[$i]})) {
			$PACmatrix{$samples[$i]} = &pacCalc($Zscore,$subset,$samples[$i]);
		}else{
			warn "Warn:: $0 :: getPAC : ".$samples[$i]." is already exists.\n";
		}
	}
	return \%PACmatrix;
}

sub output () {
	my $Zscore = shift;
	my $probeSet = shift;
	my $tag = shift;
	my $handle2 = shift;
	my $class_1 = shift;
	my $class_2 = shift;
	foreach my $iprobe (@{$probeSet}) {
		print $handle2 ($tag,"\t",$iprobe);
		foreach my $isample (@{$class_1},@{$class_2}) {
			print $handle2 ("\t",$Zscore->{$iprobe}{$isample});
		}
		print $handle2 ("\n");
	}
	return 1;
}

sub pacCalc () { ## the row 1 and the col 2
	my $Zscore = shift;
	my $subset = shift;
	my $sample = shift;
	my $num = scalar(@{$subset});
#	warn "pacCalc::Num: ".join("\t",@{$subset})." Sample: ".$sample;
	my @temp = ();
	for (my $i = 0; $i < scalar(@{$subset}); $i++) {
		push @temp,$Zscore->{$subset -> [$i]}{$sample}/sqrt($num);
	}
	return &getSum(@temp);
}

sub getZscore () {
	my $Raw = shift;
	my %Zscore = ();
	warn "Calculating Zscore...";
	my @temp = ();
	foreach my $ilab (keys %{$Raw}) {
		foreach my $isample (keys %{$Raw->{$ilab}}) {
			push @temp,$Raw->{$ilab}{$isample};
			#warn "PUSH: ".$ilab." ".$isample." ".$Raw->{$ilab}{$isample};
		}
	}
	my $mean = &getMean(@temp);
	my $sd = &getSD(@temp);
	#die "getZscore:: debug: SD = 0: ".join("_",@temp) unless ($sd != 0);
	foreach my $ilab (keys %{$Raw}) {
		foreach my $isample (keys %{$Raw->{$ilab}}) {
			$Zscore{$ilab}{$isample} = ($Raw->{$ilab}{$isample} - $mean) / $sd;
		}
	}
	warn "Calculating Zscore completed.\nZscore: ".scalar(keys %Zscore)."\n";
	return \%Zscore;
}

sub getDS () {
	my $addre = shift;
	my $class_1 = shift;
	my $class_2 = shift;
	my @class1 = ();
	my @class2 = ();
	map {push @class1, $addre -> {$_}} @{$class_1};
	map {push @class2, $addre -> {$_}} @{$class_2};
	return &getFscore(\@class1,\@class2);
}

sub larger () {
	if ($_[0] > $_[1]) {
		return $_[0];
	}else{
		return $_[1];
	}
}

sub getFscore () {
	my ($class1,$class2) = @_;
	my $mean1 = &getMean(@{$class1});
	my $mean2 = &getMean(@{$class2});
#	warn "Mean1: ".$mean1."\tMean2: ".$mean2."\n";
	my $num1 = scalar(@{$class1});
	my $num2 = scalar(@{$class2});
#	warn "Num1: ".$num1."\tNum2: ".$num2."\n";
	my $value = $num1*(($mean1 - $mean2)**2);
	my $value2 = 0;
	map {$value2 += ($_-$mean1)**2} @{$class1};
	if ($value2 == 0) {$value2 = 1}
	if ($mean1 - $mean2 < 0) {
		return -($value/$value2);
	}else{
		return $value/$value2;
	}
}

sub getmCORGsScoreD () {
	my ($selectSet,$genei,$Zscore,$fscore,$class) = @_;
	my @classi = ();my $sum = 0;
	map {
		if (exists $Zscore->{$genei}{$_}) {
			push @classi,$Zscore->{$genei}{$_};
		}else{
			die "ERRO::getmCORGsScoreD::$genei $_ is not exists in the Zscore matrix";
		}
		
		} @{$class};
	foreach my $ilab (@{$selectSet}) { 
		my @classA = ();
		map {
			if (exists $Zscore->{$ilab}{$_}) {
				push @classA,$Zscore->{$ilab}{$_};
			}else{
				die "ERRO::getmCORGsScoreD::$ilab $_ is not exists in the Zscore matrix";
			}
			
			} @{$class};
		$sum += abs(&pearson(\@classA,\@classi));
	}
	die "ERRO::getmCORGsScoreD" unless (exists $fscore->{$genei});
	return abs($fscore->{$genei}) - ($sum/scalar(@{$selectSet}));
}


sub getmCORGsScoreQ () {
	my ($selectSet,$genei,$Zscore,$fscore,$class) = @_;
	my @classi = ();my $sum = 0;
	map {
		if (exists $Zscore->{$genei}{$_}) {
			push @classi,$Zscore->{$genei}{$_};
		}else{
			die "ERRO::getmCORGsScoreQ::$genei $_ is not exists in the Zscore matrix";
		}
		
		} @{$class};
	foreach my $ilab (@{$selectSet}) { 
		my @classA = ();
		map {
			if (exists $Zscore->{$ilab}{$_}) {
				push @classA,$Zscore->{$ilab}{$_};
			}else{
				die "ERRO::getmCORGsScoreQ::$ilab $_ is not exists in the Zscore matrix";
			}
			
			} @{$class};
		$sum += abs(&pearson(\@classA,\@classi));
	}
	die "ERRO::getmCORGsScoreQ" unless (exists $fscore->{$genei});
	return abs($fscore->{$genei})/($sum/scalar(@{$selectSet}));
}

sub rank () { ### 1 in descending 2 in asending
	my $addre = shift;
	my $dir = shift;
	my %tmp = ();
	my @ranked = ();
	foreach my $key (keys %{$addre}) {
		push @{$tmp{$addre->{$key}}},$key;
	}
	if ($dir == 2) {
		foreach my $key (sort {$a <=> $b} keys %tmp) {
			push @ranked,@{$tmp{$key}};
		}
	}elsif ($dir == 1) {
		foreach my $key (sort {$b <=> $a} keys %tmp) {
			push @ranked,@{$tmp{$key}};
		}
	}else{
		die "rank::erro dir $dir.\n";
	}
	return \@ranked;
}

sub getMean () {
	my $num = scalar(@_);
	die "getMean:: no element found.\n" unless ($num > 0);
	my $sum = &getSum(@_);
	return $sum/$num;
}

sub getSum () {
	my $total = 0;
	map {$total += $_} @_;
#	warn "getSum:: Sum = 0: ".join("_",@_) unless ($total != 0);
	return $total;
}

sub getSD () {
	my $mean = &getMean(@_);
	my $num = scalar(@_);
	my @temp = ();
	map {push @temp,($_-$mean)**2} @_;
#	warn "getSD: Temp size: ".scalar(@temp)." [0]: ".$temp[0]." [-1]: ".$temp[-1];
	my $tmp = &getSum(@temp);
	if ($tmp != 0) {
		return sqrt((1/($num-1))*$tmp);
	}else{
		return 0;
	}
}

sub pearson () {
	my ($ad1,$ad2) = @_;
	my @x = @$ad1;
	my @y = @$ad2;
	my $x2=0;
	my $y2=0;
	my $xy=0;
	my $s_x=0;
	my $s_y=0;
	####################
		if ($#x != $#y) {
		die "cout x != cout y\n";
	}
	my $n =$#x + 1;
	######################
	for(my $i=0;$i<$n;$i++) {
		$x2 += $x[$i]*$x[$i];
		$y2 += $y[$i]*$y[$i];
		$xy += $x[$i]*$y[$i];
		$s_x += $x[$i];
		$s_y += $y[$i];
	}
	my $low = sqrt($n*$x2-$s_x*$s_x)*sqrt($n*$y2-$s_y*$s_y);
	if ($low != 0) {
		return ($n*$xy-$s_x*$s_y)/$low;
	}else{
		return 0;
	}
}

sub getFscoreMatrix() {
	my $subset = shift;
	my $Zscore = shift;
	my $class_1 = shift;
	my $class_2 = shift;
	my %fMatrix = ();
	foreach my $ilab (@$subset) {
		my @class1 = ();
		my @class2 = ();
		if (exists $Zscore->{$ilab}) {
			map {
				if (exists $Zscore->{$ilab}{$_}) {
					push @class1,$Zscore->{$ilab}{$_};
				}else{
					warn "Warn:: getTtestMatrix: Class 1: <".$ilab."> <".$_."> do not exists.\n";
				}
				} @{$class_1};
			map {
				if (exists $Zscore->{$ilab}{$_}) {
					push @class2,$Zscore->{$ilab}{$_};
				}else{
					warn "Warn:: getTtestMatrix: Class 2: <".$ilab."> <".$_."> do not exists.\n";
				}
				} @{$class_2};
		}else{
			warn "Warn:: getTtestMatrix: Labels: <".$ilab."> do not exists.\n";
			die "getTtestMatrix:: inject ".$ilab." into Zscore Failed.\n" unless (&repair($Zscore,$ilab,$class_1,$class_2));
			warn "getTtestMatrix:: inject ".$ilab." into Zscore Success.\n";
			return &getFscoreMatrix($subset,$Zscore,$class_1,$class_2);
		}
		if (not (exists $fMatrix{$ilab})) {
			$fMatrix{$ilab} = &getFscore(\@class1,\@class2);
		}else{
			warn "Warn:: getTtestMatrix: ".$ilab." is already exists in the tMatrix.\n";
		}
	}
	return \%fMatrix;
}

sub repair () {
	my $Zscore = shift;
	my $ilab = shift;
	my $class_1 = shift;
	my $class_2 = shift;
	map {
		if (not (exists $Zscore -> {$ilab}{$_})) {
			$Zscore -> {$ilab}{$_} = 0;
		}else{
			warn "repair:: ".$ilab." ".$_." is already exists.\n";
		}
	} @$class_1;
	map {
		if (not (exists $Zscore -> {$ilab}{$_})) {
			$Zscore -> {$ilab}{$_} = 0;
		}else{
			warn "repair:: ".$ilab." ".$_." is already exists.\n";
		}
	} @$class_2;
	return 1;
}


sub load_annotation () {
	my $ifile = shift;
	my %annotation = ();
	warn "load_annotation:: loading ".$ifile." ...";
	open (IN,$ifile) or die "$!";
	while ($line = <IN>) {
		$line=~s/\r//g;$line=~s/\b//g;$line=~s/\n//g;
		next unless ($line !~/^\#/);
		next unless ($line !~/^\"Probe/);
		my @ele = split ("\",\"",$line);
		map {$_=~s/\"//g} @ele;
		if (not (exists $symbol{$ele[0]})) {
			$symbol{$ele[0]} = $ele[14];
		}else{
			warn "load_annotation:: ".$ele[0]." is already exists in the HASH:symbol.\n";
		}
		#die "load_annotation:: annotation injection failed.\n" unless (&parse_interpro($ele[34],$ele[0],\%annotation));
		do {
			$ele[33] = "Others";
			} unless ($ele[33] =~ /\w/);
		die "load_annotation:: annotation injection failed.\n" unless (&parse_pathway($ele[33],$ele[0],\%annotation));
	}
	close IN or die "$!";
	warn "load_annotation:: loading ".scalar(keys %annotation)." domain cats.\n";
	return \%annotation;
}

sub parse_interpro () {
	my $tag = shift;
	my $id = shift;
	my $anno = shift;
	my @ele = split (/\/\/\//,$tag);
	foreach my $key (@ele) {
		my @temp = split(/\/\//,$key);
		next unless (scalar(@temp) == 3);
		map {$_=~s/ //g} @temp;
		$anno -> {$temp[0]}{$id} = $temp[1];
	}
	return 1;
}

sub parse_pathway () {
	my $tag = shift;
	my $id = shift;
	my $anno = shift;
	my @ele = split (/\/\/\//,$tag);
	foreach my $key (@ele) {
		my @temp = split(/\/\//,$key);
		@temp = grep {$_ !~/ GenMAPP/} @temp;
		map {$_ =~s/_GenMAPP//g} @temp;
		map {$_ =~s/_KEGG//g} @temp;
		map {$_ =~s/_Reactome//g} @temp;
		map {$_ =~s/ //g} @temp;	
		map {
			#warn "Found <".$_."> Pathway.\n";sleep(1);
			$anno -> {$_}{$id} = $_} @temp;
	}
	return 1;
}

sub load_groups () {
	my $ifile = shift;
	my %groups = ();
	my $iter = 0;
	open (IN,$ifile) or die "$!";
	while ($line = <IN>) {
		$line =~s/\r//g;$line =~s/\b//g;$line =~s/\n//g;
		my @ele = split("\t",$line);
		$groups{$ele[1]}{$ele[0]} = 1;
		$SampleLable{$ele[0]} = $ele[1];
	}
	close IN or die "$!";
	return \%groups;
}

sub load_Expression () {
	my $ifile = shift;
	my %head = ();
	my %Raw = ();
	warn "loading $ifile...";
	open (IN,$ifile) or die "$!";
	$line = <IN>;
	chomp($line);
	my @ele = split("\t",$line);
	for (my $i = 0;$i <= $#ele; $i++) {
		$head{$i} = $ele[$i];
		$head{$ele[$i]} = $i;
		#warn "Found identifier: ".$ele[$i]." ".$i;
	}
	while ($line = <IN>) {
		chomp($line);
		@ele = split("\t",$line);
		do{
			warn "Warn: ".$ele[0]." identifier do not found.Skipped.\n";
			next; 
		}unless ($ele[0]=~/\w/);
		for (my $i = 1;$i <= $#ele;$i++) {
			if (not (exists $Raw{$ele[0]}{$head{$i}})) {
				$Raw{$ele[0]}{$head{$i}} = $ele[$i];
			}else{
				warn "Warn: ".$ele[0]." ".$head{$i}." is already exists.\n";
			}
		}
	}
	close IN or die "$!";

	warn "Loading Completed.\nData: ".scalar(keys %Raw)." signals identified.\n";
	return \%Raw;
}