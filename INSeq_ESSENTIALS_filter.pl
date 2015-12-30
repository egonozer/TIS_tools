#!/usr/bin/perl

my $version = 0.1;

use strict;
use warnings;

$|++;

my $usage = "
INSeq_filter.pl [options] -m <local minumum value from essentialgenes_densityplot.png> -c essentialgenes_alloutputmerged.tsv -t gene_alloutputmerged.tsv 

Using INSeq data analyzed by ESSENTIALS, filters genes likely to be essential
for growth in media and outputs list of genes above and below a given log2 fold-
change.

Outputs three files:

<prefix>_essentialgenes_filtered.tsv: A list of all genes considered to be
                                      essential for growth in media
<prefix>_genes_filtered.tsv: A list of all genes, but with counts for all genes
                             essential for growth in media replaced with \"NA\"
<prefix>_genes_treatment_related.tsv: A list of only non-essential genes meeting
                                      the fold-change and p-value cutoffs.

Required:
  -m    local minimum value identified from essentialgenes_densityplot.png file
        from ESSENTIALS analysis of only control pool(s)
  -c    essentialgenes_alloutputmerged.tsv file from ESSENTIALS analysis of only
        control pool(s)
  -t    gene_alloutputmerged.tsv file from ESSENTIALS analysis of both control
        and target pools

Options:
  -f    lower log2 fold-change cutoff
        (default: -1)
  -F    upper log2 fold-change cutoff
        (default: 1)
  -p    maximum p-value
        (default: 0.05)
  -o    output file prefix
        (default \"out\")
        
";

use Getopt::Std;
use vars qw( $opt_m $opt_c $opt_t $opt_f $opt_F $opt_p $opt_o );
getopts('m:c:t:f:F:p:o:');

die $usage unless ($opt_m and $opt_c and $opt_t);

my $localmin    = $opt_m;
my $ctsv        = $opt_c;
my $ttsv        = $opt_t;
my $lo          = $opt_f ? $opt_f : -1;
my $hi          = $opt_F ? $opt_F : 1;
my $maxpval     = $opt_p ? $opt_p : 0.05;
my $pref        = $opt_o ? $opt_o : "out";

# get list of essential genes
open (my $cin, "<$ctsv") or die "ERROR: Can't open $ctsv: $!\n";
my %essen;
my ($tcount, $ecount) = (0) x 2;
my $fccol;
open (my $eout, "> $pref\_essentialgenes_filtered.tsv");
while (my $line = <$cin>){
    chomp $line;
    my @tmp = split("\t", $line);
    unless ($line =~ m/^\d/){
        print $eout "$line\n";
        for my $i (0 .. $#tmp){
            if ($tmp[$i] eq "logFC"){
                $fccol = $i;
                last;
            }
        }
        next;
    }
    #	Row.names	Input_1_1_control.x	Input_2_2_control.x	Input_3_3_control.x	unique_flanking_sequences	Input_1_1_control.y	Input_2_2_control.y	Input_3_3_control.y	expected_reads	logConc	logFC	P.Value	adj.P.Val	Start	Stop	Strand	Length	PID	Gene	Product
    $tcount++;
    my ($id, $fc) = ($tmp[1], $tmp[$fccol]);
    if ($fc eq "NA"){
        $essen{$id}++;
        $ecount++;
        print $eout "$line\n";
        next;
    }
    if ($fc < $localmin){
        $essen{$id}++;
        $ecount++;
        print $eout "$line\n";
    }
}
close ($cin);
close ($eout);
print STDERR "$tcount total genes in $opt_c, $ecount genes are essential for growth in media (log2 fold-change less than $localmin)\n";

# filter output to remove essential genes and output list of genes with changes above or below fold-change threshold
open (my $tin, "<$ttsv") or die "ERROR: Can't open $ttsv: $!\n";
my @list;
my $ccount = 0;
my ($fcol, $pcol);
my $head;
open (my $fout, "> $pref\_genes_filtered.tsv");
while (my $line = <$tin>){
    chomp $line;
    my @tmp = split("\t", $line);
    unless ($line =~ m/^\d/){
        print $fout "$line\n";
        $head = $line;
        for my $i (0 .. $#tmp){
            if ($tmp[$i] eq "logFC"){
                $fcol = $i;
            }
            if ($tmp[$i] eq "adj.P.Val"){
                $pcol = $i;
            }
        }
        next;
    }
    my ($id, $fc, $pval) = ($tmp[1], $tmp[$fcol], $tmp[$pcol]);
    if ($essen{$id}){
        for my $i (2 .. $pcol){
            $tmp[$i] = "NA";
        }
        print $fout join("\t", @tmp), "\n";
        next;
    }
    print $fout "$line\n";
    next if ($fc eq "NA");
    if (($fc < $lo or $fc > $hi) and $pval < $maxpval){
        push @list, [@tmp];
        $ccount++;
    }
}
close ($tin);
close ($fout);
print STDERR "$ccount genes with log2 fold-change less than $lo or greater than $hi AND with adjusted p-value less than $maxpval\n";
@list = sort{$a->[$fcol] <=> $b->[$fcol]} @list;
open (my $cout, "> $pref\_genes_treatment_related.tsv");
print $cout "$head\n";
foreach my $line (@list){
    my @tmp = @{$line};
    print $cout join("\t", @tmp), "\n";
}
close ($cout);
