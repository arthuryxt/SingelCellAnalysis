#!/usr/bin/perl -w
use strict;
die "Usage: $0   \"output_filename\" \"ERCC_info\"  \"strand==plus_or_minus_or_both\"   \(debug==0\)" if (@ARGV < 2);
my $fileout=$ARGV[0];
# samtools view accepted_hits.bam | cut -f1-4 > accepted_hits.info
# cut -d\_ -f2-4 accepted_hits.info | sort | uniq -c | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}' > accepted_hits.info1
my $filein=$ARGV[1]; # 5-columns: UMI_Barcode FLAG ERCC_id pos redun
my $strand="plus";
if (scalar(@ARGV) > 2) {$strand=$ARGV[2]; }
if (($strand ne "plus") and ($strand ne "minus") and ($strand ne "both")) {
    die "strand info can only be one of the following three: 1)plus 2)minus 3)both\n";
}
my $debug=0;
if (scalar(@ARGV) > 3) {$debug=$ARGV[3]; }

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
my %Sample;
while (<IN1>) {
    chomp;
    # AAAAAA_11	0	ERCC-00096	1045	3
    my @a=split("\t",$_);
    my @b=split("_",$a[0]);
    my $cFlag=dec2bin($a[1]);
    my @aFlag=split('',$cFlag);
    if ($aFlag[6] eq 1){ $sum_minus++; }
    else { $sum_plus++; }
    $Sample{$b[1]}++;
    $ERCC{$a[2]}++;
    if (($strand eq "minus") and ($aFlag[6] eq 1)) {
        $uniq{$a[2]}{$b[1]}++;
    }
    elsif (($strand eq "plus") and ($aFlag[6] eq 0)) {
        $uniq{$a[2]}{$b[1]}++;
    }
    elsif ($strand eq "both") {
        $uniq{$a[2]}{$b[1]}++;
    }
}
close IN1;

print "Plus Strand = ",$sum_plus,"\nMinus Strand = ",$sum_minus,"\n";
open(OUT, ">".$fileout) or die "Cannot open output file.\n";
my $header="geneID";
foreach my $id (sort keys %Sample) { $header=$header."\t".$id; }
print OUT $header,"\n";
foreach my $gene (sort keys %ERCC) {
    $header=$gene;
    foreach my $id (sort keys %Sample) {
        if (exists $uniq{$gene}{$id}) {
            $header=$header."\t".$uniq{$gene}{$id};
        }
        else {
            $header=$header."\t0";
        }
    }
    print OUT $header,"\n";
}
close OUT;


