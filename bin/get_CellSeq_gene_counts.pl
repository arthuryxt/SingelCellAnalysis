#!/usr/bin/perl -w
use strict;
die "Usage: $0   \"output_filename\" \"alignment_log\"  \"strand==plus_or_minus_or_both\"  \"barcode_pos==3\"  \"use_UMI==0_or_1\"  \"UMI_pos==1\"   \(debug==0\)" if (@ARGV < 2);
my $fileout=$ARGV[0];   # accepted_hits.ref.exon.cnt
my $filein=$ARGV[1];    # accepted_hits.ref.exon.log
my $strand="plus";
if (scalar(@ARGV) > 2) {$strand=$ARGV[2]; }
if (($strand ne "plus") and ($strand ne "minus") and ($strand ne "both")) {
    die "strand info can only be one of the following three: 1)plus 2)minus 3)both\n";
}
my $barcode_pos=3;  # 0-based cell-barcodes which are used to mark different cells
if (scalar(@ARGV) > 3) {$barcode_pos=$ARGV[3]; }
my $useUMI=0;
if (scalar(@ARGV) > 4) {$useUMI=$ARGV[4]; }
if ($useUMI eq 0) { print "count all reads and ignore the UMIs even if they exist.\n"; }
else {print "count only non-redundant reads using the UMI info.\n"}
my $UMI_pos=1;  # 0-based
if (scalar(@ARGV) > 5) {$UMI_pos=$ARGV[5]; }
my $debug=0;
if (scalar(@ARGV) > 6) {$debug=$ARGV[6]; }
if ($useUMI eq 0) { $UMI_pos=0; }
my %uniq;
my %sample;
my %Redun;
my %Gene;
open(IN, $filein) or die "Cannot open alignment_log file : ".$filein."!\n";
while (<IN>) {
    chomp;
    my @a=split("\t",$_);
    # SN541:333:HTJFJADXX:1:1211:10044:6134_TGCATG_CTAGAC_35	1	-	391268	50M	NM_001099460___2___6	391018	391300	+
    my @b=split(/\_/,$a[0]);
    my @c=split(/\_\_\_/,$a[5]);
    if (scalar(@b) eq 3) {
        # (there is no UMI)
        if (($strand eq "plus") and ($a[2] eq $a[8])) {
            my $aln_info=join("\t",$b[0],$a[1],$a[2],$a[3]);
            #if (exists $uniq{$b[3]}{$c[0]}{$aln_info}) {
            #    # reduncancy, ignore
            #}
            if (exists $Redun{$b[0]}) {
                $Redun{$b[0]}++;
                $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}=join("\t",$b[0],$a[5],$a[6],$a[7]);
                $sample{$b[$barcode_pos]}++;
                $Gene{$c[0]}++;
            }
            else {
                $Redun{$b[0]}=1;
                $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}=join("\t",$b[0],$a[5],$a[6],$a[7]);
                $sample{$b[$barcode_pos]}++;
                $Gene{$c[0]}++;
            }
        }
        elsif (($strand eq "minus") and ($a[2] ne $a[8])) {
            my $aln_info=join("\t",$b[0],$a[1],$a[2],$a[3]);
            #if (exists $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}) {
            #    # reduncancy, ignore
            #}
            if (exists $Redun{$b[0]}) {
                $Redun{$b[0]}++;
                $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}=join("\t",$b[0],$a[5],$a[6],$a[7]);
                $sample{$b[$barcode_pos]}++;
                $Gene{$c[0]}++;
            }
            else {
                $Redun{$b[0]}=1;
                $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}=join("\t",$b[0],$a[5],$a[6],$a[7]);
                $sample{$b[$barcode_pos]}++;
                $Gene{$c[0]}++;
            }
        }
        elsif ($strand eq "both") {
            my $aln_info=join("\t",$b[0],$a[1],$a[2],$a[3]);
            #if (exists $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}) {
            #    # reduncancy, ignore
            #}
            if (exists $Redun{$b[0]}) {
                $Redun{$b[0]}++;
                $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}=join("\t",$b[0],$a[5],$a[6],$a[7]);
                $sample{$b[$barcode_pos]}++;
                $Gene{$c[0]}++;
            }
            else {
                $Redun{$b[0]}=1;
                $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}=join("\t",$b[0],$a[5],$a[6],$a[7]);
                $sample{$b[$barcode_pos]}++;
                $Gene{$c[0]}++;
            }
        }
    }
    elsif ($useUMI eq 0) {
        # (don't use UMI)
        if (($strand eq "plus") and ($a[2] eq $a[8])) {
            my $aln_info=join("\t",$b[$UMI_pos],$a[1],$a[2],$a[3]);
            #if (exists $uniq{$b[3]}{$c[0]}{$aln_info}) {
            #    # reduncancy, ignore
            #}
            if (exists $Redun{$b[0]}) {
                $Redun{$b[0]}++;
                $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}=join("\t",$b[0],$a[5],$a[6],$a[7]);
                $sample{$b[$barcode_pos]}++;
                $Gene{$c[0]}++;
            }
            else {
                $Redun{$b[0]}=1;
                $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}=join("\t",$b[0],$a[5],$a[6],$a[7]);
                $sample{$b[$barcode_pos]}++;
                $Gene{$c[0]}++;
            }
        }
        elsif (($strand eq "minus") and ($a[2] ne $a[8])) {
            my $aln_info=join("\t",$b[$UMI_pos],$a[1],$a[2],$a[3]);
            #if (exists $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}) {
            #    # reduncancy, ignore
            #}
            if (exists $Redun{$b[0]}) {
                $Redun{$b[0]}++;
                $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}=join("\t",$b[0],$a[5],$a[6],$a[7]);
                $sample{$b[$barcode_pos]}++;
                $Gene{$c[0]}++;
            }
            else {
                $Redun{$b[0]}=1;
                $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}=join("\t",$b[0],$a[5],$a[6],$a[7]);
                $sample{$b[$barcode_pos]}++;
                $Gene{$c[0]}++;
            }
        }
        elsif ($strand eq "both") {
            my $aln_info=join("\t",$b[$UMI_pos],$a[1],$a[2],$a[3]);
            #if (exists $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}) {
            #    # reduncancy, ignore
            #}
            if (exists $Redun{$b[0]}) {
                $Redun{$b[0]}++;
                $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}=join("\t",$b[0],$a[5],$a[6],$a[7]);
                $sample{$b[$barcode_pos]}++;
                $Gene{$c[0]}++;
            }
            else {
                $Redun{$b[0]}=1;
                $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}=join("\t",$b[0],$a[5],$a[6],$a[7]);
                $sample{$b[$barcode_pos]}++;
                $Gene{$c[0]}++;
            }
        }
    }
    else{
        if (($strand eq "plus") and ($a[2] eq $a[8])) {
            my $aln_info=join("\t",$b[$UMI_pos],$a[1],$a[2],$a[3]);
            if (exists $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}) {
                # reduncancy, ignore
            }
            elsif (exists $Redun{$b[0]}) {
                $Redun{$b[0]}++;
                $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}=join("\t",$b[0],$a[5],$a[6],$a[7]);
                $sample{$b[$barcode_pos]}++;
                $Gene{$c[0]}++;
            }
            else {
                $Redun{$b[0]}=1;
                $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}=join("\t",$b[0],$a[5],$a[6],$a[7]);
                $sample{$b[$barcode_pos]}++;
                $Gene{$c[0]}++;
            }
        }
        elsif (($strand eq "minus") and ($a[2] ne $a[8])) {
            my $aln_info=join("\t",$b[$UMI_pos],$a[1],$a[2],$a[3]);
            if (exists $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}) {
                # reduncancy, ignore
            }
            elsif (exists $Redun{$b[0]}) {
                $Redun{$b[0]}++;
                $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}=join("\t",$b[0],$a[5],$a[6],$a[7]);
                $sample{$b[$barcode_pos]}++;
                $Gene{$c[0]}++;
            }
            else {
                $Redun{$b[0]}=1;
                $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}=join("\t",$b[0],$a[5],$a[6],$a[7]);
                $sample{$b[$barcode_pos]}++;
                $Gene{$c[0]}++;
            }
        }
        elsif ($strand eq "both") {
            my $aln_info=join("\t",$b[$UMI_pos],$a[1],$a[2],$a[3]);
            if (exists $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}) {
                # reduncancy, ignore
            }
            elsif (exists $Redun{$b[0]}) {
                $Redun{$b[0]}++;
                $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}=join("\t",$b[0],$a[5],$a[6],$a[7]);
                $sample{$b[$barcode_pos]}++;
                $Gene{$c[0]}++;
            }
            else {
                $Redun{$b[0]}=1;
                $uniq{$b[$barcode_pos]}{$c[0]}{$aln_info}=join("\t",$b[0],$a[5],$a[6],$a[7]);
                $sample{$b[$barcode_pos]}++;
                $Gene{$c[0]}++;
            }
        }
    }
}
close IN;
open(OUT, ">".$fileout) or die "Cannot open output file!\n";
my $result="geneID";
foreach my $s (sort keys %sample) { $result=$result."\t".$s; }
print OUT $result,"\n";
foreach my $g (sort keys %Gene) {
    $result=$g;
    foreach my $s (sort keys %sample) {
        my $cnt=0;
        foreach my $aln (keys %{$uniq{$s}{$g}}) {
            my @a=split("\t",$uniq{$s}{$g}{$aln});
            $cnt+=1/$Redun{$a[0]};
        }
        $result=$result."\t".sprintf("%.3f",$cnt);
    }
    print OUT $result,"\n";
}
