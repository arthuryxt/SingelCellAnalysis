#!/usr/bin/perl -w
use strict;
die "Usage: $0   \"output_filename\"  \"Barcode_info\"  \"input_fastq\"  \"prefix-len\"  \(split_file==0_or_1\)" if (@ARGV < 3);
my $fileout=$ARGV[0];
my $Barcode=$ARGV[1];
my $filein=$ARGV[2];
my $preflen=1000; # output the first $preflen nt for each read
if (scalar(@ARGV) > 3) {$preflen=$ARGV[3]; }
my $split_file=0;
if (scalar(@ARGV) > 4) {$split_file=$ARGV[4]; }

if ($split_file eq 0) {
    my %uniq;
    open(IN, $Barcode) or die "Cannot open Barcode file.\n";
    while (<IN>) {
        chomp;
        my @a=split("\t",$_);
        # SN541:333:HTJFJADXX:1:1107:1155:2118	NCTTTT	ACGTAC	32
        # OR
        # SRR567993.7	GACACCGC	19
        if (scalar(@a) eq 4) {
            $uniq{$a[0]}=join("_",$a[1],$a[2],$a[3]);
        }
        else {
            $uniq{$a[0]}=join("_",$a[0],$a[1],$a[2]);
        }
    }
    close IN;
    print "reading Barcode file done.\n";
    open(IN1, $filein) or die "Cannot open input_fastq file.\n";
    open(OUT, ">".$fileout) or die "Cannot open output file.\n";
    while (<IN1>) {
        chomp;
        s/^@//;
        my @t1=split(" ",$_);
        my $id="@".$t1[0]."_".$uniq{$t1[0]};
        my $seq=<IN1>;
        chomp $seq;
        my $qua=<IN1>;
        $qua=<IN1>;
        print OUT $id,"\n",substr($seq,0,$preflen),"\n+\n",substr($qua,0,$preflen),"\n";
    }
    close IN1;
    close OUT;
}
else {
    my %uniq;
    open(IN, $Barcode) or die "Cannot open Barcode file.\n";
    while (<IN>) {
        chomp;
        my @a=split("\t",$_);
        # SN541:333:HTJFJADXX:1:1107:1155:2118	NCTTTT	ACGTAC	32
        # OR
        # SRR567993.7	GACACCGC	19
        if (scalar(@a) eq 4) {
            $uniq{$a[3]}{$a[0]}=join("_",$a[1],$a[2],$a[3]);
        }
        else {
            $uniq{$a[2]}{$a[0]}=join("_",$a[0],$a[1],$a[2]);
        }
    }
    close IN;
    print "reading Barcode file done.\n";
    foreach my $sample (sort keys %uniq){
        open(IN1, $filein) or die "Cannot open input_fastq file.\n";
        open(OUT, ">".$fileout."_".$sample) or die "Cannot open output file.\n";
        while (<IN1>) {
            chomp;
            s/^@//;
            my @t1=split(" ",$_);
            if (exists $uniq{$sample}{$t1[0]}) {
                my $id="@".$t1[0]."_".$uniq{$sample}{$t1[0]};
                my $seq=<IN1>;
                chomp $seq;
                my $qua=<IN1>;
                $qua=<IN1>;
                print OUT $id,"\n",substr($seq,0,$preflen),"\n+\n",substr($qua,0,$preflen),"\n";
            }
            else {
                my $seq=<IN1>;
                $seq=<IN1>;
                $seq=<IN1>;
            }
        }
        close IN1;
        close OUT;
    }    
}



