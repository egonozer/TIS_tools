#!/usr/bin/perl

use strict;
use warnings;

my $usage = "
INSeq_ESSENTIALS_filter_compare.pl <file_of_files>

Compares two or more ESSENTIALS experiments to determine what is the same / different

Required:
  <file_of_files>: Text file conttaining paths to treatment-related genes files
                   output by INSeq_ESSENTIALS_filter.pl
                   Format:
                   /path/to/expt1_genes_treatment_related.tsv<tab>experiment_1
                   /path/to/expt2_genes_treatment_related.tsv<tab>experiment_2
                   etc.

";

die $usage unless (@ARGV);

my $fof = $ARGV[0];
my @order;
my %seen;
my %results;
open (my $fin, "<$fof") or die "ERROR: Can't open $fof: $!\n";
while (my $fline = <$fin>){
    chomp $fline;
    my ($path, $id) = split("\t", $fline);
    die "ERROR: IDs must be unique ($id seen more than once)\n" if $seen{$id};
    $seen{$id}++;
    push @order, $id;
    open (my $in, "<$path") or die "ERROR: Can't open $path: $!\n";
    my $logfc_col;
    while (my $line = <$in>){
        chomp $line;
        my @tmp = split("\t", $line);
        unless ($logfc_col){
            for my $i (0 .. $#tmp){
                my $val = $tmp[$i];
                if ($val eq "logFC"){
                    $logfc_col = $i;
                    last;
                }
            }
            next;
        }
        my ($idx, $lid, $logfc, $logCPM, $pval, $fdr, $gene, $prod) = @tmp[0,1,$logfc_col,$logfc_col+1,$logfc_col+2,$logfc_col+3,$logfc_col+9,$logfc_col+10];
        @{$results{$lid}{$id}} = ($logfc, $logCPM, $pval, $fdr);
        @{$results{$lid}{'gene_and_product'}} = ($gene, $prod);
    }
    close ($in);
}
close ($fin);

my @array;
my $avgpos;
foreach my $lid (sort keys %results){
    my ($gene, $prod) = @{$results{$lid}{'gene_and_product'}};
    my @tmp;
    push @tmp, $lid;
    my @logfcs;
    foreach my $id (@order){
        my ($logfc, $logCPM, $pval, $fdr) = ("-") x 4;
        if ($results{$lid}{$id}){
            ($logfc, $logCPM, $pval, $fdr) = @{$results{$lid}{$id}};
            push @logfcs, $logfc;
        }
        push @tmp, $logfc;
    }
    my $num = scalar(@logfcs);
    my $sum = 0;
    foreach (@logfcs){
        $sum += $_;
        $num-- if $_ eq "-";
    }
    my $avg = sprintf("%.5f", $sum / $num);
    push @tmp, $avg;
    $avgpos = $#tmp;
    push @tmp, ($gene, $prod);
    push @array, [@tmp];
}

@array = sort{$a->[$avgpos] <=> $b->[$avgpos]}@array;
print "locus_id";
foreach my $id (@order){
    print "\t$id\_logFC";
}
print "\tavg_logFC\tgene\tproduct\n";
foreach my $slice (@array){
    my @tmp = @{$slice};
    print join("\t", @tmp), "\n";
}
