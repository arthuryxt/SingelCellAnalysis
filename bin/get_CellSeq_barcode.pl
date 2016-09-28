#!/usr/bin/perl -w
use strict;
die "Usage: $0   \"output_filename\"  \"Barcode\"  \"input_fasta\"  \"len_barcode\"  \"len_UMI\"  \"error==0\"   \(UMI_first==0_or_1\)" if (@ARGV < 3);
my $fileout=$ARGV[0];
my $Barcode=$ARGV[1];
my $filein=$ARGV[2];
my $lenBC=6;
if (scalar(@ARGV) > 3) {$lenBC=$ARGV[3];}
my $lenUMI=6;
if (scalar(@ARGV) > 4) {$lenUMI=$ARGV[4];}
my $MM=1;
if (scalar(@ARGV) > 5) {$MM=$ARGV[5]; }
if (($MM < 0) or ($MM > 6)) {
    die "The number of errors tolerated in the barcode sequence should between 0 and 6.\n";
}
my $UMI_first=0;    #
if (scalar(@ARGV) > 6) {$UMI_first=$ARGV[6]; }
if (($UMI_first ne 0) and ($UMI_first ne 1)) {
    die "the parameter UMI_first should be either 0 or 1.\n";
}

my %uniq;
open(IN, $Barcode) or die "Cannot open Barcode file.\n";
my $cnt=1;
while (<IN>) {
    chomp;
    $uniq{$_}=$cnt;
    $cnt++;
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

open(IN1, $filein) or die "Cannot open input_fasta file.\n";
open(OUT, ">".$fileout) or die "Cannot open output file.\n";
while (<IN1>) {
    chomp;
    s/^@//;
    my @t1=split(" ",$_);
    my $id=$t1[0];
    my $seq=<IN1>;
    chomp $seq;
    my $UMI=substr($seq,0,$lenUMI);
    my $CB=substr($seq,$lenUMI,$lenBC);
    if ($UMI_first eq 0) {
        $CB=substr($seq,0,$lenBC);
        $UMI=substr($seq,$lenBC,$lenUMI);
    }
    if (exists $uniq{$CB}) {
        # find a hit to exact CellBarcode
        print OUT join("\t",$id,$UMI,$CB,$uniq{$CB}),"\n";
    }
    elsif ($MM > 0) {
        my $tmpCB="";
        my $randomMax=999999;
        my $tmpdist=$randomMax;
        foreach my $b (keys %uniq) {
            if (length($b) < 1) { next;}
            my $dist=levenshtein($b, $CB);
            if ($dist <= $MM) {
                if ($MM eq 1) {
                    $tmpdist = $dist;
                    $tmpCB=$b;
                    last;
                }
                else {
                    if ($tmpdist > $dist) {
                        $tmpdist = $dist;
                        $tmpCB=$b;
                    }
                }
            }
        }
        if ($tmpdist ne $randomMax){
            # found a match
            print OUT join("\t",$id,$UMI,$CB,$uniq{$tmpCB}),"\n";
        }
        else {
            print OUT join("\t",$id,$UMI,$CB,"NA"),"\n";
        }
    }
    $seq=<IN1>;
    $seq=<IN1>;
}
close IN1;
close OUT;


