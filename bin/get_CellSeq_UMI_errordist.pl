#!/usr/bin/perl -w
use strict;
# collect UMI mismatch-frequency for each genomic location (with at least two different UMIs), so that we can later decide whether correction for sequencing/PCR errors is needed
die "Usage: $0   \"output_filename\" \"SAM_file\"  \"strand==plus_or_minus_or_both\"  \"barcode_pos==3\"  \"UMI_pos==1\"   \(debug==0\)" if (@ARGV < 2);
my $fileout=$ARGV[0];
my $filein=$ARGV[1];
my $strand="no";      # will NOT consider strand info
if (scalar(@ARGV) > 2) {$strand=$ARGV[2]; }
if (($strand ne "no") and ($strand ne "yes") ) {
    die "strand info can only be one of the following two: 1)yes 2)no \n";
}
my $barcode_pos=3;  # 0-based cell-barcodes which are used to mark different cells, by default pos2==sequence and pos3==number
if (scalar(@ARGV) > 3) {$barcode_pos=$ARGV[3]; }
my $UMI_pos=1;  # 0-based
if (scalar(@ARGV) > 4) {$UMI_pos=$ARGV[4]; }
my $debug=0;
if (scalar(@ARGV) > 5) {$debug=$ARGV[5]; }
my %uniq;
my %sample;
my %Redun;

sub countN($){
    my @A=split //, lc shift; # lower case
    my %letter;
    my $Nr=scalar(@A);
    for(@A) {$letter{$_}++;}
    my $top=0;
    my $sec=0;
    foreach my $NT (sort{$letter{$b} <=> $letter{$a}} keys %letter) {
        if ($top eq 0) {$top=$letter{$NT};}
        else {$sec=$letter{$NT}; last;}
    }
    if ((($top / $Nr) > 0.6) or (($top + $sec) / $Nr > 0.8)) {
        return 1;   # yes, simple repeat, should not consider
    }
    else {
        return 2;   # seems all right
    }
}

sub pFLAG {
    my $str=shift;
    my $seq=shift;
    my $cFlag=dec2bin($str);
    my @aFlag=split('',$cFlag);
    my $strand="+"; # genomic strandness, by default
    if ($aFlag[3] eq 1) { $strand="*"; return $strand; }    # ignore "second fragment"
    if ($aFlag[6] eq 1) { $strand="-"; }
    my $simple=2;
    # only check if reported as "secondary alignment"
    if ($aFlag[2] eq 1) { $simple=countN($seq); }
    if ($simple eq 2) { return $strand;}
    else { $strand="*"; return $strand;}
}
sub bin2dec {
    return unpack("N", pack("B32", substr("0" x 32 . shift, -32)));
}
sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    #$str =~ s/^0+(?=\d)//; # otherwise you'll get leading zeros return $str;
    my $decode = substr($str,-11);
    return $decode;
}

open(IN, $filein) or die "Cannot open SAM file : ".$filein."!\n";
while (<IN>) {
    chomp;
    my @a=split("\t",$_);
    # SN541:340:H3YYHADXX:2:1214:4767:23437_TATATT_GAGAGA_94	272	1	3041505	0	51M	*	0	0	TTTGGTTTTTGGTTTTTGGTTTTTGGTTTTTGGTTTTTGGTTTTTGGTTTT	D;EEA2=8..(.8/8??8)D:B1:)1@;?<27FGBG;IFC<3:=DDDA81=	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:51	YT:Z:UU	NH:i:20	CC:Z:=	CP:i:4751264	HI:i:0
    my $seq_strand=pFLAG($a[1],$a[9]);
    if ($seq_strand eq "*") {next;}     # discard simple-composite reads
    my @b=split(/\_/,$a[0]);
    my $gpos=join("_",$a[2],$a[3]);
    if ($strand eq "yes") {
        $gpos=join("_",$a[2],$a[3],$seq_strand);
    }
    my $weight=1;
    for (my $i=11; $i<scalar(@a); $i++) { if($a[$i]=~/NH/){ my @tmp=split(/\:/,$a[$i]); $weight=1/$tmp[-1]; } }
    $uniq{$b[$barcode_pos]}{$gpos}{$b[$UMI_pos]}+=$weight;
}
close IN;

sub levenshtein($$){
    my @A=split //, lc shift; # lower case
    my @B=split //, lc shift;
    my @W=(0..@B);
    my ($i, $j, $cur, $next);
    for $i (0..$#A){
	$cur=$i+1;
	for $j (0..$#B){
	    $next=min( $W[$j+1]+1, $cur+1, ($A[$i] ne $B[$j])+$W[$j] );
            $W[$j]=$cur;
            $cur=$next;
	}
	$W[@B]=$next;
    }
    return $next;
}

sub min($$$){
    if ($_[0] < $_[2]){ pop @_; } else { shift @_; }
    return $_[0] < $_[1]? $_[0]:$_[1];
} 

open(OUT, ">".$fileout) or die "Cannot wite to file : ".$fileout."!\n";
foreach my $sample (sort keys %uniq) {
    my $sample_pos=0;
    my $sample_pos_singleumi=0;
    my $sample_reads=0;
    my $sample_reads_singleumi=0;
    foreach my $gpos (sort keys %{$uniq{$sample}}) {
        $sample_pos++;
        my %Dist;
        my %cnt;
        my $f=0;            # number of pairs among all non-self umi
        my $s=0;            # sum of pair-wise distance among all non-self umi
        my $uniqumi=0;      # number of unique umi
        my $totalumi=0;     # total number of umi
        my $maumi="";       # most abundant umi
        my $umidist=0;      # sum of distance of sub-umi to the maumi
        foreach my $umi (sort{$uniq{$sample}{$gpos}{$b} <=> $uniq{$sample}{$gpos}{$a}} keys %{$uniq{$sample}{$gpos}}) {
            $uniqumi++;
            $totalumi+=$uniq{$sample}{$gpos}{$umi};
            $cnt{$umi}=$uniq{$sample}{$gpos}{$umi};
            if ($maumi eq "") { $maumi=$umi; }
            foreach my $mykey (sort{$cnt{$b} <=> $cnt{$a}} keys %cnt) {
                if ($mykey ne $umi) {
                    my $pairs=$cnt{$umi} * $cnt{$mykey} ;
                    #if (($cnt{$umi} eq 1) or ($cnt{$mykey} eq 1)) { $pairs=$cnt{$umi} * $cnt{$mykey}; } 
                    $f += $pairs;
                    $s += $pairs * levenshtein($umi, $mykey);
                }
            }
            # assess the distance of all (sub) umi to the most abundant umi
            if ($umi ne $maumi) {
                $umidist+=levenshtein($umi, $maumi);
            }
            
        }
        if ($uniqumi > 1) { print OUT join("\t",$sample,$uniqumi,$totalumi,$umidist,$f,sprintf("%.3f",$s/$f)),"\n"; $sample_reads+=$totalumi; }
        else { $sample_pos_singleumi++; $sample_reads_singleumi+=$cnt{$maumi}; }
    }
    print OUT join("\t","#".$sample,$sample_pos,$sample_pos_singleumi,$sample_reads,$sample_reads_singleumi),"\n";
}