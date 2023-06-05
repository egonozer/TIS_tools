#!/usr/bin/perl

my $version = 0.3;

## Changes from v0.2:
## Calculate number of shared and total insertion sites. Will separately calculate for input pool(s) containing the word "input" or as marked by the user in the input file
## Outputs stats to a file by default instead of just to STDOUT. 
## Removed File::Which dependency. Added subroutine to test for whether executable is in PATH that uses only core Perl modules
## Minor code cleanup

## Changes from v0.1:
## Changed name from essentials_read_preprocess.pl to INSeq_read_preprocess.pl (probably will rename again at some point)
## Add code to allow non-Unix style line endings (i.e. if somebody uses a .tsv produced by Excel)
## Allow user to set path to bowtie and bowtie-build and check to see if it's present (uses File::Which)
## Added option to run bowtie alignments multithreaded

use warnings;
use strict;
#use File::Which;
use File::Spec::Functions qw ( catfile path );

$|++;

## Inspired in part by Andrew Goodman's INSeq Analysis Pipeline

my $usage = "
INSeq_read_preprocess.pl (version: $version)

For INSeq experiments, peforms the following processing steps on raw sequencing
reads:
1) identify reads with appropriate MmeI sequence and remove all upstream
   sequence --> will keep all resulting 23 - 25 bp reads 
   (expected if 51 bp reads produced)
2) align reads to the reference genome
3) filter alignments

Prerequisites:
bowtie (version 0.12.8 or higher)

Required:
  -r    File of demultiplexed read files corresponding to separate pools,
        conditions, and/or replicates. Read files must be in fastq format and 
        can be gzipped. 
        OPTIONAL: Mark input pool(s) with 'i' in a third column. Otherwise will 
        guess input are files that have 'input' (case insensitive) at the 
        beginning of the pool_ID name.
        
        File should have the following format and be tab delimited with one 
        pool / read file per line:
        
        path/to/read_file.fastq(.gz) <tab> pool_ID 
  
  -g    Reference genome sequence file, in fasta format. For best results, this
        file must include ALL sequences present in the organism (i.e. chromsomes
        and plasmids)

Options:
  -t    Trim sequence
        (default: 'GGGGACTTATCATCCAACCTGTTA')
  -m    Number of trim sequence mismatches allowed. Maximum is 3
        (default: 1)
  -x    Maximum number of genome alignment mismatches allowed
        (default: 1)
  -o    Output format for count files. Possible options include:
            'wiggl': Files will take the wiggle format that can be used as
                     input to Essentials
            'inseq': Files will take the format of bowtiemap_processed.txt files
                     produced by Goodman's INSeq_analysis software
            'essen': Files will take the format of allta_split_*.counts.txt
                     files produced by Essentials
        (default: 'wiggl')
  -s    Minimum number of total reads at a site (left+right) required for output
        (default: 3)
  -b    Path to directory containing bowtie and bowtie-build. For example:
        /Users/myname/applications/bowtie_folder
        (default: assumes this directory is in your PATH)
  -p    Number of threads
        (default: 1)
  -P    Output stats files prefix. Output files will be titled 
        '<prefix>.inseq_read_preprocess_stats.txt'
        '<prefix>.inseq_read_preprocess_site_counts.txt'
        (default: 'out')

";

use Getopt::Std;
use vars qw( $opt_r $opt_g $opt_t $opt_m $opt_x $opt_o $opt_s $opt_h $opt_d $opt_b $opt_p $opt_P );
getopts('r:g:t:m:x:o:s:hd:b:p:P:');

die $usage unless ($opt_r and $opt_g);

my $reads   = $opt_r if $opt_r;
my $reffile = $opt_g if $opt_g;

my $tn      = $opt_t ? $opt_t : "GGGGACTTATCATCCAACCTGTTA";
my $maxmm   = defined $opt_m ? $opt_m : 1; #added 'defined' so users can enter 0
$maxmm = 3 if $maxmm > 3;
my $maxrmm  = defined $opt_x ? $opt_x : 1;
my $outform = $opt_o ? $opt_o : "wiggl";
my $minrd   = $opt_s ? $opt_s : 3;
my $maxdiff = $opt_d ? $opt_d : 1;
if ($maxdiff < 1){
    $maxdiff = 1;
}
my $bowpath = $opt_b if $opt_b;
my $threads = $opt_p ? $opt_p : 1;
my $outpref = $opt_P ? $opt_P : "out";

$outform = lc($outform);
unless ($outform =~ m/inseq|essen|wiggl/){
    print STDERR "WARNING: Output format '$outform' does not match one of the possible options (inseq/essen/wiggl). Defaulting to 'wiggl'\n";
    $outform = "wiggl";
}

## Count TA sites in the genome;
my %tasite_count;
open (my $gin, "<$reffile") or die "ERROR: Can't open $reffile: $!\n";
my ($seq, $id);
while (my $line = <$gin>){
    chomp $line;
    $line =~ s/\s//g;
    if ($line =~ m/^>/){
        if ($seq){
            $seq = uc($seq);
            my @array = split(//, $seq);
            for my $i (0 .. $#array - 1){
                my $bi = $array[$i].$array[$i+1];
                if ($bi eq "TA"){
                    $tasite_count{$id} += 1;
                }
            }
            $tasite_count{$id} += 1 if $array[$#array].$array[0] eq "TA"; ## Assumes the sequence is circular
            $seq = "";
        }
        $id = substr($line, 1);
        next;
    }
    $seq .= $line;
}
close $gin;
if ($seq){
    $seq = uc($seq);
    my @array = split(//, $seq);
    for my $i (0 .. $#array - 1){
        my $bi = $array[$i].$array[$i+1];
        if ($bi eq "TA"){
            $tasite_count{$id} += 1;
        }
    }
    $tasite_count{$id} += 1 if $array[$#array].$array[0] eq "TA"; ## Assumes the sequence is circular
    $seq = "";
}

## process the trim sequence
$tn = uc($tn); 
my @tnarray = split(//, $tn);
my $tnleng = length($tn);
my $tnmm = $tn;
if ($maxmm > 0){
    my @mmarray;
    for my $i (0 .. ($tnleng - $maxmm)){
        my $temp = $tn;
        substr($temp, $i, 1, ".");
        if ($maxmm > 1){
            for my $j ($i+1 .. $tnleng - ($maxmm-1)){
                my $temp1 = $temp;
                substr($temp1, $j, 1, ".");
                if ($maxmm == 3){
                    for my $k ($j+1 .. $tnleng - ($maxmm-2)){
                        my $temp2 = $temp1;
                        substr($temp2, $k, 1, ".");
                        push @mmarray, $temp2;
                    }
                } else {
                    push @mmarray, $temp1;
                }
            }
        } else {
            push @mmarray, $temp;
        }
    }
    $tnmm = join("|", @mmarray);
}


## check if bowtie and bowtie-build are available
my $bow_loc = is_path("bowtie");
my $bb_loc = is_path("bowtie-build");
if ($bowpath){
    if (-e "$bowpath/bowtie"){
        $bow_loc = "$bowpath/bowtie";
    } else {
        print STDERR "WARNING: Could not find 'bowtie' at $bowpath. Searching PATH...\n";
    }
    if (-e "$bowpath/bowtie-build"){
        $bb_loc = "$bowpath/bowtie-build";
    }
    unless (-e "$bowpath/bowtie-build"){
        print STDERR "WARNING: Could not find 'bowtie-build' at $bowpath. Searching PATH...\n";
    }
}
die "ERROR: Could not find bowtie in PATH. Make sure bowtie is installed and executable or use the '-b' option to give the path to the bowtie directory.\n" if !$bow_loc;
die "ERROR: Could not find bowtie-build in PATH. Make sure bowtie is installed and executable or use the '-b' option to give the path to the bowtie directory.\n" if !$bb_loc;

my %results;

## Process reads 
print STDERR "Starting read processing...\n";
my %sorted;
my %read_serials;
my @pools;
my $proc_reads = 1;
my %inpool_guess;
my %inpool;
if ($proc_reads){
    my $read_serial = 0;
    open (my $in, "<$reads");
    while (my $line = <$in>){
        chomp $line;
        next if $line =~ m/^\s*#/; #Skip commented lines
        my ($path, $pool_id, $is_in) = split("\t", $line);
        push @pools, $pool_id;
        if ($is_in and $is_in =~ m/^i/i){
            $inpool{$pool_id} = 1;
        }
        if ($pool_id =~ m/^input/i){
            $inpool_guess{$pool_id} = 1;
        }
        my $rin;
        if ($path =~ m/\.gz$/){
            open ($rin, "gzip -cd $path | ") or die "ERROR: Can't open $path: $!\n";
        } else {
            open ($rin, "<", $path) or die "ERROR: Can't open $path: $!\n";
        }
        my $total_reads = 0;
        my %reads_w_bc;
        my %reads_w_tn;
        my ($tot_tn, $tot_ln) = (0) x 2;
        my %tn_mms;
        print STDERR "Processing $pool_id read $total_reads ... ";
        my %test;
        while (<$rin>){
            my $x1 = $_;
            chomp (my $seq = <$rin>);
            my $x2 = <$rin>;
            my $x3 = <$rin>;
            $total_reads++;
            print STDERR "\rProcessing $pool_id read $total_reads ... " if ($total_reads % 10000 == 0);
            my $tn_match;
            if ($seq =~ m/$tn/){ #try without mismatches first. Should go much faster
                $tn_match = $+[0] - 2; #returns the end position of the sequence match, but includes the TA
                $tot_tn++;
                $tn_mms{0}++;
            } elsif ($seq =~ m/$tnmm/) {
                $tn_match = $+[0] - 2; #returns the end position of the sequence match, but includes the TA
                $tot_tn++;
                $tn_mms{1}++;
            }
            next unless $tn_match;
            my $finseq = substr($seq, $tn_match);
            my $seqleng = length($finseq);
            #$test{$seqleng}++;
            if ($seqleng < 23 or $seqleng > 25) {  # Eventually make this user-controllable
                #$test{$seqleng}++;
                next;
            }
            $tot_ln++;
            unless ($sorted{$finseq}){
                $read_serial++;
                $read_serials{$read_serial} = $finseq;
                $sorted{$finseq}{'serialnumber'} = $read_serial;
            }
            $sorted{$finseq}{$pool_id}{'count'}++;
            $sorted{$finseq}{$pool_id}{'leng'} = $seqleng; ## Vast majority should be 24 bp including the TA with this protocol (if read trimming is off)
        }
        close ($rin);
        print STDERR "\rProcessing $pool_id read $total_reads ... Done!\n";
        print STDERR "\tTotal reads with trim sequence: $tot_tn\n";
        
        #my @tarray;
        #foreach my $lleng (sort{$a <=> $b} keys %test){
        #    push @tarray, "$lleng:$test{$lleng}";
        #}
        #print STDERR join(",", @tarray), "\n\n";
        
        $results{$pool_id}{'reads_total'} = $total_reads;
        $results{$pool_id}{'reads_tn'} = $tot_tn;
        $results{$pool_id}{'reads_tn_perfect'} = $tn_mms{0};
        $results{$pool_id}{'reads_right_length'} = $tot_ln;
    }
    close ($in);
}
unless (%inpool){
    if (%inpool_guess){
        print STDERR "\nGuessing input pools:\n";
        print STDERR join("\n", sort(keys %inpool_guess)), "\n";
        %inpool = %inpool_guess;
    }
}

#create bowtie index files
print STDERR "\nBuilding bowtie index ...";
my $command = "$bb_loc $reffile temp_ispp >/dev/null 2>&1";
my $status = system($command);
die "ERROR: command '$command' exited with status $status\n" if $status;
print STDERR " Done\n";

my %refs;
my %align_stats;
my %all_st;
my %all_corr;

#multithread bowtie alignments
#output reads, run bowtie, process results
open (my $out, ">temp_reads.fasta") or die "ERROR: Can't open output file 'temp_reads.fasta': $!\n";
foreach my $seq (sort keys %sorted){
    #my $count = $sorted{$seq}{'count'};
    #my $leng = $sorted{$seq}{'leng'};
    my $serial = $sorted{$seq}{'serialnumber'};
    #print $out ">R:$leng:$count\n$seq\n";
    print $out ">$serial\n$seq\n";
}
close ($out);
my $aout = "temp.bowtie";
print STDERR "\tRunning bowtie  ...\n";
#could just pipe results of bowtie directly into processing, but for the purposes of keeping a trail of steps, will output alignments, then process
$command = "$bow_loc -m 1 --best --strata -a --fullref -n $maxrmm -l 25 temp_ispp -f temp_reads.fasta -p $threads > $aout 2>/dev/null";
$status = system($command);
die "ERROR: command '$command' exited with status $status\n" if $status;

## Delete bowtie index files
unlink glob "temp_ispp*ebwt";

#process the bowtie alignment output
open (my $in, "<$aout") or die "ERROR: Can't open bowtie alignment file '$aout': $!\n";
my %sites;
my %read_seqs; #This will be to store actual aligning read sequences at each site, if comparison to an assembly sequence is done.

#following are for end-of-run stats
my $num_aligned = 0;
my $num_aligned_to_TA = 0;
my $num_aligned_perfect = 0;
my $num_aligned_mismatch = 0;
my $num_aligned_shifted = 0;
my $num_aligned_filtered = 0;
my $num_not_TA = 0;

#parse the bowtie alignment and assign read counts to positions and alignment types.
print STDERR "\tParsing bowtie alignments ...\n";
while (my $line = <$in>){
    chomp $line;
    my ($id, $dir, $ref, $pos, $seq, $qual, $x1, $mm) = split("\t", $line);
    my $serial = $id;
    my $readseq = $read_serials{$serial};
    my ($skip, $side, $site, $type) = (0) x 4;
    if ($dir eq "+"){
        $side = "right";
        $site = $pos + 1;
        my $head = substr($seq, 0, 2);
        if ($head eq "TA"){
            if ($mm){
                my @mms = split(",", $mm); ### This seems dumb. Bowtie should be filtering hits with too many mismatches, I think?
                if (scalar @mms > $maxrmm){
                    $skip = "max_mm";
                    next;
                }
                $type = "imperf";
            } else {
                $type = "perf";
            }
        } else {
            $skip = "not_ta"; ## Just going to filter all these because they seem to be mostly garbage 
        }
    } else {
        $side = "left";
        $site = $pos + length($seq) - 1;  #bowtie is 0-based
        my $tail = substr($seq, -2);
        if ($tail eq "TA"){
            if ($mm){
                my @mms = split(",", $mm);
                if (scalar @mms > $maxrmm){
                    $skip = "max_mm";
                    next;
                }
                $type = "imperf";
            } else {
                $type = "perf";
            }
        } else {
            $skip = "not_ta"; ## Just going to filter all these because they seem to be mostly garbage 
        }
    }
    foreach my $pool_id (@pools){
        next unless $sorted{$readseq}{$pool_id};
        my $count = $sorted{$readseq}{$pool_id}{'count'};
        $results{$pool_id}{'reads_aligned'} += $count;
        if ($skip eq "0"){
            $sites{$pool_id}{$ref}{$site}{$side}{$type} += $count;
            $results{$pool_id}{$type} += $count;
        } else {
            #not_ta = did not align to TA site
            #max_mm = number of mismatches exceeded maximum set by user
            $results{$pool_id}{$skip} += $count;
        }
    }
}
close ($in);


#print "Total reads aligned: $num_aligned\n";
#print "Total reads aligned to TA sites: $num_aligned_to_TA\n";
#print "Total reads aligned to non-TA sites (removed): $num_not_TA\n";
#print "Total reads perfectly aligned to TA sites: $num_aligned_perfect\n";
#print "Total reads aligned to TA sites with mismatch: $num_aligned_mismatch\n";
#print "Total reads aligned to TA sites with >$maxrmm mismatches (removed): $num_aligned_filtered\n";
#print "\n";


#print "Total reads shifted to TA sites: $num_aligned_shifted\n\n";

#$align_stats{$bar} = ([$num_aligned, $num_aligned_to_TA, $num_aligned_perfect, $num_aligned_mismatch, $num_aligned_shifted]);

#process and correct transposon insertion site counts
my %site_counts;
my %in_site_counts;
my %site_counts_gte;
my %in_site_counts_gte;
foreach my $pool_id (@pools){
    print STDERR "Outputting $pool_id alignments\n";
    foreach my $ref (sort keys %{$sites{$pool_id}}){
        $refs{$ref}++;
        my $outfile = "$pool_id.$ref.wiggle";
        $outfile = "allta_split_$pool_id\_$ref.counts.txt" if $outform eq "essen";
        $outfile = "INSEQ_experiment.scarf.bowtiemapped_processed.txt_$pool_id\_$ref" if $outform eq "inseq";
        open (my $wig_out, "> $outfile") or die "ERROR: Can't open output file '$outfile': $!\n";
        
        #for status update
        my $num_sites = keys %{$sites{$pool_id}{$ref}};
        my $site_count = 0;
        my $site_gte_count = 0;
        print STDERR "\tProcessing site #$site_count (of $num_sites)";    
        
        my @site_type_counts;
        #set each base value to 0;
        for my $i (1 .. 8){
            $site_type_counts[$i] = 0;
        }
        #my $beg_mm_correct_count = 0; #for keeping track of how many sites are corrected for a beginning mismatch
        
        foreach my $site (sort {$a <=> $b} keys %{$sites{$pool_id}{$ref}}){
            $site_count++;
            print STDERR "\r\tProcessing site #$site_count (of $num_sites)" if $site_count % 100 == 0;
            
            
            #my ($l_perf, $l_imperf, $r_perf, $r_imperf) = (0) x 4;
            my $l_perf = defined $sites{$pool_id}{$ref}{$site}{'left'}{'perf'} ? $sites{$pool_id}{$ref}{$site}{'left'}{'perf'} : 0;
            my $l_imperf = defined $sites{$pool_id}{$ref}{$site}{'left'}{'imperf'} ? $sites{$pool_id}{$ref}{$site}{'left'}{'imperf'} : 0;
            my $r_perf = defined $sites{$pool_id}{$ref}{$site}{'right'}{'perf'} ? $sites{$pool_id}{$ref}{$site}{'right'}{'perf'} : 0;
            my $r_imperf = defined $sites{$pool_id}{$ref}{$site}{'right'}{'imperf'} ? $sites{$pool_id}{$ref}{$site}{'right'}{'imperf'} : 0;
            my $tot = ($l_perf + $l_imperf + $r_perf + $r_imperf);

            $site_counts{$ref}{$site}++;
            $in_site_counts{$ref}{$site}++ if (%inpool and defined $inpool{$pool_id});

            next if $tot < $minrd; #skip any sites with less than 3 (or user input) reads
            $site_gte_count++;
            
            $site_counts_gte{$ref}{$site}++;
            $in_site_counts_gte{$ref}{$site}++ if (%inpool and defined $inpool{$pool_id});

            #assign site to site type and correct counts if possible
            my $assigned;
            #my ($l_out, $r_out) = (0,0);
            my $l_out = $l_perf + $l_imperf;
            my $r_out = $r_perf + $r_imperf;
            
            #output read counts in format requested
            my $tot_out = $l_out + $r_out;
            my $outline = "-$site\t$l_out\n$site\t$r_out\n";
            $outline = "$l_out\t-$site\n$r_out\t$site\n" if $outform eq "essen";
            $outline = ">$ref\t$site\t$l_out\t$r_out\t$tot_out\n" if $outform eq "inseq";
            print $wig_out "$outline";
        }
        close ($wig_out);
        print STDERR "\r\tProcessing site #$site_count (of $num_sites) ... Done!\n";
        print STDERR "\t\t$ref: # sites total: $site_count, # sites with at least $minrd reads: $site_gte_count\n";
        $results{$pool_id}{$ref}{'site_count'} = $site_count;
        $results{$pool_id}{$ref}{'site_gte_count'} = $site_gte_count;
        
        #print "# of sites of each type (# of sites modified/filtered):\n";
        #print "\t#1 both perfect, no imperfect (mismatch or shifted): $site_type_counts[1] ($corrected_site_counts[1])\n";
        #print "\t#2 both perfect, minor imperfect on one or both sides: $site_type_counts[2] ($corrected_site_counts[2])\n";
        #print "\t#3 both perfect, major imperfect on one or both sides: $site_type_counts[3] ($corrected_site_counts[3])\n";
        #print "\t#4 perfect on one side, imperfect on the other: $site_type_counts[4] ($corrected_site_counts[4])\n";
        #print "\t#5 imperfect only on both sides: $site_type_counts[5] ($corrected_site_counts[5])\n";
        #print "\t#6 perfect on one side, 0 reads on other: $site_type_counts[6] ($corrected_site_counts[6])\n";
        #print "\t#7 imperfect on one side, 0 reads on other: $site_type_counts[7] ($corrected_site_counts[7])\n";
        #print "\t#8 unassigned: $site_type_counts[8] ($corrected_site_counts[8])\n\n";
        #$site_type_counts[0] = $site_gte_count;
        #$all_st{$ref}{$bar} = [@site_type_counts];
        #$all_corr{$ref}{$bar} = [@corrected_site_counts];
    }
}
print STDERR "\n";

open (my $sout, ">$outpref.inseq_read_preprocess_stats.txt") or die "ERROR: Can't open stats output file: $!\n";
my $string = "pool\tinput?\ttotal_reads\treads_w_tn\treads_w_tn_perfect\treads_correct_length\treads_aligned_to_refs\treads_aligned_TA\treads_aligned_TA_perfect";
my @refids;
foreach my $ref (sort keys %refs){
    push @refids, $ref;
    $string .= "\t$ref\_sites\t$ref\_sites_gte_$minrd";
}
$string .= "\n";
###                        0          1           2                 3                4         5     6      7
my @result_types = qw(reads_total reads_tn reads_tn_perfect reads_right_length reads_aligned perf not_ta max_mm);
foreach my $pool_id (@pools){
    my @results;
    foreach my $type (@result_types){
        my $val = defined $results{$pool_id}{$type} ? $results{$pool_id}{$type} : 0;
        push @results, $val;
    }
    foreach my $ref (@refids){
        my $site_count = defined $results{$pool_id}{$ref}{'site_count'} ? $results{$pool_id}{$ref}{'site_count'} : 0;
        my $site_gte_count = defined $results{$pool_id}{$ref}{'site_gte_count'} ? $results{$pool_id}{$ref}{'site_gte_count'} : 0;
        push @results, ($site_count, $site_gte_count);
    }
    my $reads_aligned_ta = $results[4] - $results[6];
    my $is_input = "N";
    $is_input = "Y" if (%inpool and $inpool{$pool_id});
    $string .= "$pool_id\t$is_input\t" . join("\t", @results[0..4]) . "\t$reads_aligned_ta\t$results[5]\t" . join("\t", @results[8 .. $#results]) . "\n" ;
}
output($string);
close ($sout);
$string = "";

print "\n";

open ($sout, ">$outpref.inseq_read_preprocess_site_counts.txt") or die "ERROR: Can't open site_counts output file: $!\n";
$string .= "ID\tTA_sites\tpools\tsites_w_ins\tsites_w_ins_all_pools\tsites_w_ins_gte${minrd}rds\tsites_w_ins_all_pools_gte${minrd}rds\n";
foreach my $ref (sort keys %refs){
    my $total_sites_any = scalar(%{$site_counts{$ref}});
    my $total_sites_any_gte = scalar(%{$site_counts_gte{$ref}});
    my $in_sites_any = scalar(%{$in_site_counts{$ref}});
    my $in_sites_any_gte = scalar(%{$in_site_counts_gte{$ref}});
    my ($total_sites_all, $total_sites_all_gte, $in_sites_all, $in_sites_all_gte) = (0) x 4;
    foreach my $site (keys %{$site_counts{$ref}}){
        $total_sites_all++ if $site_counts{$ref}{$site} == scalar(@pools);
    }
    foreach my $site (keys %{$site_counts_gte{$ref}}){
        $total_sites_all_gte++ if $site_counts_gte{$ref}{$site} == scalar(@pools);
    }
    foreach my $site (keys %{$in_site_counts{$ref}}){
        $in_sites_all++ if $in_site_counts{$ref}{$site} == scalar(%inpool);
    }
    foreach my $site (keys %{$in_site_counts_gte{$ref}}){
        $in_sites_all_gte++ if $in_site_counts_gte{$ref}{$site} == scalar(%inpool);
    }
    #$string = "\nTotal TA sites: $tasite_count{$ref}\n";
    my $ta_count = $tasite_count{$ref};
    $string .= "$ref\t$ta_count\tall\t$total_sites_any\t$total_sites_all\t$total_sites_any_gte\t$total_sites_all_gte\n";
    if (%inpool){
        $string .= "$ref\t$ta_count\tinput\t$in_sites_any\t$in_sites_all\t$in_sites_any_gte\t$in_sites_all_gte\n";
    }
    output($string);
}
close($sout);

        #$results{$pool_id}{'reads_total'} = $total_reads;
        #$results{$pool_id}{'reads_tn'} = $tot_tn;
        #$results{$pool_id}{'reads_tn_perfect'} = $tn_mms{0};
        #$results{$pool_id}{'reads_right_length'} = $tot_ln;
        #$results{$pool_id}{'reads_aligned'} += $count;
        #$results{$pool_id}{$type} += $count; #perf / imperf
        ##not_ta = did not align to TA site
        #    #max_mm = number of mismatches exceeded maximum set by user
        #    $results{$pool_id}{$skip} += $count;
        #$results{$pool_id}{$ref}{'site_gte_count'}++;
        #$results{$pool_id}{$ref}{'site_count'} = $site_count;



#output summary
#print "ID\tbarcode\treads_aligned\treads_align_TA\tperfect\tmismatch\tshifted\n";
#foreach my $slice (@bc_list){
#    my ($id, $bar) = @{$slice};
#    if ($align_stats{$bar}){
#        my @array = @{$align_stats{$bar}};
#        print "$id\t$bar\t", join("\t", @array), "\n";
#    }
#}
#print "\nreference\tID\tbarcode\ttotal_sites\ttype1\t1_corrected\ttype2\t2_corrected\ttype3\t3_corrected\ttype4\t4_corrected\ttype5\t5_corrected\ttype6\t6_corrected\ttype7\t7_corrected\ttype8\t8_corrected\n";
#foreach my $ref (sort keys %refs){
#    foreach my $slice (@bc_list){
#        my ($id, $bar) = @{$slice};
#        print "$ref\t$id\t$bar";
#        my @st = @{$all_st{$ref}{$bar}};
#        my @corr = @{$all_corr{$ref}{$bar}};
#        for my $i (0 .. 8){
#            print "\t$st[$i]";
#            next if $i == 0;
#            print "\t$corr[$i]";
#        }
#        print "\n";
#    }
#}


#---------------SUBROUTINES----------------
sub is_path {
    ## Subroutine based on StackOverflow post by Sinan Unur (https://stackoverflow.com/a/8243770)
    my $exe = shift;
    my @path = path;
    my @pathext = ( q{} );
    if ($^O eq 'MSWin32'){
        push @pathext, map { lc } split /;/, $ENV{PATHEXT};
    }
    for my $dir (@path){
        for my $ext (@pathext){
            my $f = catfile $dir, "$exe$ext";
            return ($f) if -x $f;
        }
    }
    return();
}

sub output {
    my $ostring = shift;
    print "$ostring";
    print $sout "$ostring";
    return;
}
