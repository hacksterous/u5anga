open(INPUT, "<tithi-input.txt") or die "Couldn't open file tithi-input.txt, $!";

while($line = <INPUT>) {
	$line =~ s/^\s+|\s+$//g;
	#print "$line\n";
	if ($line =~ /^(\d+)\/(\d+)\/(\d+).*\s(\d+):(\d+)\s([A|P])M\s+\d+\.\s+(\S+)\s.*$/) {
		$hr = $4;
		if ($6 eq "P" && $4 < 12) {
			$hr += 12;
		} elsif ($6 eq "A" && $4 == 12) {
			$hr = "0";
		}
		$yr = 2000 + $3;
		print "$2-$1-$yr,$hr:$5,$7\n";
	}
}
