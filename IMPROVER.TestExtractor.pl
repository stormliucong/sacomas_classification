#!usr/bin/perl -w 
use strict;
use warnings;

my $line;
my %symbol = ();
my %SampleLable = ();

die "perl $0 <Annotation File><ExpressionMatrix><Features><Group File>" unless $#ARGV== 3;

&main();

sub main () {
	my $annotation = &load_annotation($ARGV[0]);
	my $expression_matrix = &load_Expression($ARGV[1]);
	my $groups = &load_groups($ARGV[3]);
	my $Zscore = &getZscore($expression_matrix);
	my $features = &load_features($ARGV[2]);
	open (OUT,">".$ARGV[1].".".$0.".Test.csv") or die "$!";
	print OUT ("Label,Class,",join(",",sort {$a cmp $b} keys %{$features}),"\n");
	foreach my $iclass (keys %SampleLable) {
		print OUT ($iclass,",",$SampleLable{$iclass});
		foreach my $ifeature (sort {$a cmp $b} keys %{$features}) {
			my $matrix = &getPAC($features->{$ifeature},$Zscore);
			print OUT (",",$matrix->{$iclass});
		}
		print OUT ("\n");
	}
	close OUT or die "$!";
}


sub load_features () {
	my $ifile = shift;
	my %features = ();
	open (IN,$ifile) or die "$!";
	while ($line = <IN>) {
		$line =~s/\r//g;$line =~s/\b//g;$line =~s/\n//g;
		my @ele = split(",",$line);
		if ($ele[0] =~/\w/) {
			@{$features{$ele[0]}} = @ele[1..$#ele];
		}
	}
	close IN or die "$!";
	return \%features;
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
	return &getTtest(\@class1,\@class2);
}

sub larger () {
	if ($_[0] > $_[1]) {
		return $_[0];
	}else{
		return $_[1];
	}
}

sub getTtest () {
	my ($class1,$class2) = @_;
	my $mean1 = &getMean(@{$class1});
	my $mean2 = &getMean(@{$class2});
#	warn "Mean1: ".$mean1."\tMean2: ".$mean2."\n";
	my $num1 = scalar(@{$class1});
	my $num2 = scalar(@{$class2});
#	warn "Num1: ".$num1."\tNum2: ".$num2."\n";
	my $S12 = sqrt(((($num1-1)*(&getSD(@{$class1})**2))+(($num2-1)*(&getSD(@{$class2})**2)))/($num1+$num2-2));
#	warn "SD1: ".&getSD(@{$class1})."\tSD2: ".&getSD(@{$class2})."\n";
#	warn "S12: ".$S12."\n";
	my $tmp = ($S12*sqrt((1/$num1)+(1/$num2)));
	if ($tmp == 0) {
		return 0;
	}else{
		return (($mean1-$mean2)/$tmp);
	}
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

sub getTtestMatrix() {
	my $subset = shift;
	my $Zscore = shift;
	my $class_1 = shift;
	my $class_2 = shift;
	my %tMatrix = ();
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
			return &getTtestMatrix($subset,$Zscore,$class_1,$class_2);
		}
		if (not (exists $tMatrix{$ilab})) {
			$tMatrix{$ilab} = &getTtest(\@class1,\@class2);
		}else{
			warn "Warn:: getTtestMatrix: ".$ilab." is already exists in the tMatrix.\n";
		}
	}
	return \%tMatrix;
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