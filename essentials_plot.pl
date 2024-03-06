#!/usr/bin/perl

my $version = 0.2;

use strict;
use warnings;

#Tracks:
#0) All genes (and directions) - from genes_filtered.tsv
#1) Insertion sites (with at least ___ reads) - from wiggle files, preprocess? OR ta_alloutputmerged.tsv from full run
#2) Essential Genes - from essentialgenes_filtered.tsv
#3) Treatment-related genes - from genes_treatment_related.tsv
#4) Fold-changes (frac of min/max) - from genes_filtered.tsv, logFC

#Inputs:
#1) files listed above:
#   -genes_filtered.tsv
#   -essentialgenes_filtered.tsv
#   -genes_treatment_related.tsv
#   -wiggle file summary
#2) Genome size, in bp


my $usage = "
essentials_plot.pl

Generates xml-formatted file for plotting with cgview.jar

Required:
  -g    genes_filtered.tsv
  -e    essentialgenes_filtered.tsv
  -t    genes_treatment_related.tsv
  -i    insertion_sites.txt (from INSeq_read_preprocess.pl)
  -s    reference genome size (in bp)
  
Optional:
  -f    Display positive and negative log fold-changes relative to the maximum
        positive and negative log fold-changes, respectively
        (default: will display positive and negative changes relative to the
        maximum absolute positive or negative log fold-changes)

";

use Getopt::Std;
use vars qw( $opt_g $opt_e $opt_t $opt_i $opt_s $opt_f );
getopts('g:e:t:i:s:f');
die $usage unless ($opt_g and $opt_e and $opt_t and $opt_i and $opt_s);

my $g_file      = $opt_g;
my $e_file      = $opt_e;
my $t_file      = $opt_t;
my $i_file      = $opt_i;
my $gensize     = $opt_s;

## read in the genes file;
open (my $gin, "<$g_file") or die "ERROR: Can't open $g_file: $!\n";
my @cols;
my @genes_plus;
my @genes_minus;
my @fc_all;
my ($minfc, $maxfc) = (0) x 2;
my $absfc = 0;
while (my $line = <$gin>){
    chomp $line;
    my @tmp = split("\t", $line);
    if ($line =~ m/^\s*Row.names\s/){
        for my $i (0 .. $#tmp){
            if ($tmp[$i] eq "Row.names"){
                $cols[0] = $i;
            }
            if ($tmp[$i] eq "logFC"){
                $cols[1] = $i;
            }
            if ($tmp[$i] eq "Start"){
                $cols[2] = $i;
            }
            if ($tmp[$i] eq "Stop"){
                $cols[3] = $i;
            }
            if ($tmp[$i] eq "Strand"){
                $cols[4] = $i;
            }
        }
        next;
    }
    my ($id, $fc, $start, $stop, $strand) = @tmp[@cols];
    push @genes_plus, ([$start, $stop]) if $strand eq "+";
    push @genes_minus, ([$start, $stop]) if $strand eq "-";
    unless ($fc eq "NA"){
        $maxfc = $fc if $fc > $maxfc;
        $minfc = $fc if $fc < $minfc;
        $absfc = abs($fc) if abs($fc) > $absfc;
        push @fc_all, ([$start, $stop, $fc]);
    }
}
close ($gin);
@genes_plus = sort{$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @genes_plus;
@genes_minus = sort{$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @genes_minus;
@fc_all = sort{$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @fc_all;

my %gene_types;

## read in treatment-related genes
open (my $tin, "<$t_file") or die "ERROR: Can't open $t_file: $!\n";
my @t_genes;
while (my $line = <$tin>){
    chomp $line;
    my @tmp = split("\t", $line);
    next if ($line =~ m/^\s*Row.names\s/);
    my ($id, $fc, $start, $stop, $strand) = @tmp[@cols];
    push @t_genes, ([$start, $stop]);
    $gene_types{"$start..$stop"} = "t";
}
close ($tin);
@t_genes = sort{$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @t_genes;

#read in essential genes;
open (my $ein, "<$e_file") or die "ERROR: Can't open $e_file: $!\n";
@cols = ();
my @e_genes;
while (my $line = <$ein>){
    chomp $line;
    my @tmp = split("\t", $line);
    if ($line =~ m/^\s*Row.names\s/){
        for my $i (0 .. $#tmp){
            if ($tmp[$i] eq "Row.names"){
                $cols[0] = $i;
            }
            if ($tmp[$i] eq "logFC"){
                $cols[1] = $i;
            }
            if ($tmp[$i] eq "Start"){
                $cols[2] = $i;
            }
            if ($tmp[$i] eq "Stop"){
                $cols[3] = $i;
            }
            if ($tmp[$i] eq "Strand"){
                $cols[4] = $i;
            }
        }
        next;
    }
    my ($id, $fc, $start, $stop, $strand) = @tmp[@cols];
    push @e_genes, ([$start, $stop]);
    $gene_types{"$start..$stop"} = "e";
}
close ($ein);
@e_genes = sort{$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @e_genes;

#read in insertion sites
open (my $iin, "<$i_file") or die "ERROR: Can't open $i_file: $!\n";
my @ins;
while (my $line = <$iin>){
    chomp $line;
    $line =~ s/\s//g;
    push @ins, $line;
}
close ($iin);
@ins = sort{$a <=> $b}@ins;

my $minimumFeatureLength = "x-small"; #default: xxx-small
my $backboneRadius = 1600;
my $rulerFont = "SansSerif,plain,40";
my $tickDensity = 0.3;
my $tickThickness = 8;
my $shortTickThickness = 8;
my $rulerPadding = 40;
my $height = 4000;
my $width = 4000;
my $featureThickness = 100;

my $ecolor = "green";
my $tcolor = "orange";
my $pluscolor = "blue";
my $minuscolor = "red";

my $legendFont = "SansSerif,plain,42";

print qq{<?xml version="1.0" encoding="ISO-8859-1"?>
<cgview sequenceLength="$gensize" showBorder="false" showShading="false" minimumFeatureLength="$minimumFeatureLength" backboneRadius="$backboneRadius" rulerFont="$rulerFont" tickDensity="$tickDensity" tickThickness="$tickThickness" shortTickThickness="$shortTickThickness" rulerPadding="$rulerPadding" height="$height" width="$width" featureThickness="$featureThickness">
};

## Gene locations
## forward
print qq{<featureSlot strand="direct">
    <feature color="red" decoration="arc">
};
my $last_stop = 0;
foreach (@genes_plus){
    my ($start, $stop) = @{$_};
    my $color = "black";
    if (my $val = $gene_types{"$start..$stop"}){
        $color = $tcolor;
        $color = $ecolor if $val eq "e";
    }
    #if ($start <= $last_stop){
    #    $start = $last_stop + 1;
    #}
    #$last_stop = $stop;
    print qq{       <featureRange start="$start" stop="$stop" color="$color"/>};
    print "\n";
}
print qq{   </feature>
</featureSlot>
};
## reverse
print qq{<featureSlot strand="reverse">
    <feature color="blue" decoration="arc">
};
$last_stop = 0;
foreach (@genes_minus){
    my ($start, $stop) = @{$_};
    my $color = "black";
    if (my $val = $gene_types{"$start..$stop"}){
        $color = $tcolor;
        $color = $ecolor if $val eq "e";
    }
    #if ($start <= $last_stop){
    #    $start = $last_stop + 1;
    #}
    #$last_stop = $stop;
    print qq{       <featureRange start="$start" stop="$stop" color="$color" />};
    print "\n";
}
print qq{   </feature>
</featureSlot>
};

## Insertion sites
print qq{<featureSlot strand="reverse" featureThickness="50">
    <feature color="black" decoration="arc">
};
foreach my $val (@ins){
    print qq{       <featureRange start="$val" stop="$val" />};
    print "\n";
}
print qq{   </feature>
</featureSlot>
};

## Fold change
print qq{<featureSlot strand="reverse" featureThickness="600">
    <feature decoration="arc">
};
$last_stop = 0;
foreach (@fc_all){
    my ($start, $stop, $fc) = @{$_};
    #if ($start <= $last_stop){
    #    $start = $last_stop + 1;
    #}
    #$last_stop = $stop;
    print STDERR "$start - $stop\n" if ($stop < $start);
    my $color = "black";
    my $barHeight = 0.02;
    my $radiusShift = 0.5;
    if ($fc > 0){
        $barHeight = abs($fc) * 0.5 / $absfc;
        $barHeight = abs($fc) * 0.5 / abs($maxfc) if $opt_f;
        $radiusShift = 0.5 - $barHeight/ 2;
        $color = $pluscolor if $fc >= 1;
    } elsif ($fc < 0){
        $barHeight = abs($fc) * 0.5 / $absfc;
        $barHeight = abs($fc) * 0.5 / abs($minfc) if $opt_f;
        $radiusShift = 0.5 + $barHeight / 2;
        $color = $minuscolor if $fc <= -1;
    } 
    print qq{       <featureRange color="$color" start="$start" stop="$stop" proportionOfThickness="$barHeight" radiusAdjustment="$radiusShift" />};
    print "\n";
}
print qq{   </feature>
</featureSlot>
};

print qq{
<legend position="upper-left" font="$legendFont">
    <legendItem text="Genes" />
    <legendItem text="Essential" drawSwatch="true" swatchColor="$ecolor" />
    <legendItem text="Treatment" drawSwatch="true" swatchColor="$tcolor" />
    <legendItem text="Neutral" drawSwatch="true" swatchColor="black" />
    <legendItem text="Log Fold Change" />
    <legendItem text="Greater than 1" drawSwatch="true" swatchColor="$pluscolor" />
    <legendItem text="Less than -1" drawSwatch="true" swatchColor="$minuscolor" />
    <legendItem text="Between -1 -and 1" drawSwatch="true" swatchColor="black" />
</legend>
};



print "</cgview>\n";






