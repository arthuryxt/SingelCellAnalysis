#!/usr/bin/perl -w
use strict;
die "Usage: $0   \"output_filename\" \"ERCC_gtf\" \"ERCC_info\"  \"strand==plus_or_minus_or_both\"   \(debug==0\)" if (@ARGV < 3);
my $fileout=$ARGV[0];
my $gtf=$ARGV[1];
# samtools view accepted_hits.bam | cut -f1-4 > accepted_hits.info
# cut -d\_ -f2-4 accepted_hits.info | sort | uniq -c | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}' > accepted_hits.info1
my $filein=$ARGV[2]; # 5-columns: UMI_Barcode FLAG ERCC_id pos redun
my $strand="plus";
if (scalar(@ARGV) > 3) {$strand=$ARGV[3]; }
if (($strand ne "plus") and ($strand ne "minus") and ($strand ne "both")) {
    die "strand info can only be one of the following three: 1)plus 2)minus 3)both\n";
}
my $debug=0;
if (scalar(@ARGV) > 4) {$debug=$ARGV[4]; }
open(IN0, $gtf) or die "Cannot open ERCC_gtf file.\n";
my %Len;
while (<IN0>) {
    chomp;
    # ERCC-00002	ERCC	exon	1	1061	0.000000	+	.	gene_id "ERCC-00002"; transcript_id "DQ459430";
    my @a=split("\t",$_);
    $Len{$a[0]}=$a[4];
}
close IN0;

sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    #$str =~ s/^0+(?=\d)//; # otherwise you'll get leading zeros return $str;
    my $decode = substr($str,-11);
    return $decode;
}

sub bin2dec {
    return unpack("N", pack("B32", substr("0" x 32 . shift, -32)));
}

open(IN1, $filein) or die "Cannot open ERCC_info file.\n";
my $sum_plus=0;
my $sum_minus=0;
my %uniq;
my %ERCC;
while (<IN1>) {
    chomp;
    # AAAAAA_11	0	ERCC-00096	1045	3
    my @a=split(" ",$_);
    my $cFlag=dec2bin($a[1]);
    my @aFlag=split('',$cFlag);
    if ($aFlag[6] eq 1){ $sum_minus++; }
    else { $sum_plus++; }
    $ERCC{$a[2]}++;
    if (($strand eq "minus") and ($aFlag[6] eq 1)) {
        $uniq{$a[2]}{$Len{$a[2]}-$a[3]}++;
    }
    elsif (($strand eq "plus") and ($aFlag[6] eq 0)) {
        $uniq{$a[2]}{$Len{$a[2]}-$a[3]}++;
    }
    elsif ($strand eq "both") {
        $uniq{$a[2]}{$Len{$a[2]}-$a[3]}++;
    }
}
close IN1;

print "Plus Strand = ",$sum_plus,"\nMinus Strand = ",$sum_minus,"\n";
open(OUT, ">".$fileout) or die "Cannot open output file.\n";
foreach my $gene (sort keys %uniq) {
    my $header=$gene;
    for(my $i=0; $i<$Len{$gene}; $i++) {
        if (exists $uniq{$gene}{$i}) {
            $header=$header."\t".$uniq{$gene}{$i};
        }
        else {
            $header=$header."\t0";
        }
    }
    print OUT $header,"\n";
}
close OUT;


