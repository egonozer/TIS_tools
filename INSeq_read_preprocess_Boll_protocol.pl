#!/usr/bin/perl

my $version = 0.2;

## Changes from v0.1:
## Changed name from essentials_read_preprocess.pl to INSeq_read_preprocess.pl (probably will rename again at some point)
## Add code to allow non-Unix style line endings (i.e. if somebody uses a .tsv produced by Excel)
## Allow user to set path to bowtie and bowtie-build and check to see if it's present (uses File::Which)
## Added option to run bowtie alignments multithreaded

use warnings;
use strict;
use File::Which;

$|++;

## Inspired in part by Andrew Goodman's INSeq Analysis Pipeline

my $usage = "
INSeq_read_preprocess.pl (version: $version)

For INSeq experiments, peforms the following processing steps on raw sequencing
reads:
2) identify reads with appropriate MmeI sequence and remove all upstream
   sequence --> will keep all resulting 23 - 25 bp reads (expected if 51 bp reads produced)
3) align reads to the reference genome
4) filter alignments

Prerequisites:
bowtie (version 0.12.8 or higher)

Required:
  -r    File of demultiplexed read files corresponding to seprate pools,
        conditions, and/or replicates. Read files must be in fastq format
        and can be gzipped. File should have the following format and be tab
        delimited with one pool / read file per line:
        
        path/to/read_file.fastq(.gz) <tab> pool_ID
  
  -g    Reference genome sequence file, in fasta format. For best results, this
        file must include ALL sequences present in the organism (i.e. chromsomes
        and plasmids)

Encouraged, but optional:
  -w    assembled sequence of parent strain that was used to produce transposon
        mutagenesis library, in fasta format

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
  -s    Minimum number of total reads at a site required for output
        (default: 3)
  -b    Path to directory containing bowtie and bowtie-build. For example:
        /Users/myname/applications/bowtie_folder
        (default: assumes this directory is in your PATH)
  -p    Number of threads
        (default: 1)

";

use Getopt::Std;
use vars qw( $opt_r $opt_g $opt_t $opt_m $opt_x $opt_o $opt_w $opt_s $opt_h $opt_d $opt_b $opt_p );
getopts('r:g:t:m:x:o:w:s:hd:b:p:');

die $usage unless ($opt_r and $opt_g);

my $reads   = $opt_r if $opt_r;
my $reffile = $opt_g if $opt_g;

my $tn      = $opt_t ? $opt_t : "GGGGACTTATCATCCAACCTGTTA";
my $maxmm   = defined $opt_m ? $opt_m : 1; #added 'defined' so users can enter 0
$maxmm = 3 if $maxmm > 3;
my $maxrmm  = defined $opt_x ? $opt_x : 1;
my $outform = $opt_o ? $opt_o : "wiggl";
my $assem   = $opt_w if $opt_w;
my $minrd   = $opt_s ? $opt_s : 3;
my $maxdiff = $opt_d ? $opt_d : 1;
if ($maxdiff < 1){
    $maxdiff = 1;
}
my $bowpath = $opt_b if $opt_b;
my $threads = $opt_p ? $opt_p : 1;

$outform = lc($outform);
unless ($outform =~ m/inseq|essen|wiggl/){
    print STDERR "WARNING: Output format '$outform' does not match one of the possible options (inseq/essen/wiggl). Defaulting to 'wiggl'\n";
    $outform = "wiggl";
}

#process the trim sequence
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


#check if bowtie and bowtie-build are available
my $bow_loc = which("bowtie");
my $bb_loc = which("bowtie-build");
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
if ($proc_reads){
    my $read_serial = 0;
    open (my $in, "<$reads");
    while (my $line = <$in>) {
        chomp $line;
        next if $line =~ m/^\s*#/; #Skip commented lines
        my ($path, $pool_id) = split("\t", $line);
        print STDERR "\tProcessing $pool_id\n";
        push @pools, $pool_id;
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
        print STDERR "\tProcessing read $total_reads ... ";
        my %test;
        while (<$rin>){
            my $x1 = $_;
            chomp (my $seq = <$rin>);
            my $x2 = <$rin>;
            my $x3 = <$rin>;
            $total_reads++;
            print STDERR "\r\tProcessing read $total_reads ... " if ($total_reads % 10000 == 0);
            my $tn_match;
            if ($seq =~ m/$tn/){ #try without mismatches first. Should go much faster
                $tn_match = $+[0] - 2; #returns the end position of the sequence match, but includds the TA
                $tot_tn++;
                $tn_mms{0}++;
            } elsif ($seq =~ m/$tnmm/) {
                $tn_match = $+[0] - 2; #returns the end position of the sequence match, but includds the TA
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
        print STDERR "\r\tProcessing read $total_reads ... Done!\n";
        print STDERR "\nTotal reads: $total_reads\n\tTotal reads with trim sequence: $tot_tn\n";
        
        
        
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

#create bowtie index files
print STDERR "\nBuilding bowtie index ...";
my $command = "$bb_loc $reffile temp_ispp >/dev/null 2>&1";
my $status = system($command);
die "ERROR: command '$command' exited with status $status\n" if $status;
print STDERR " Done\n";

#read in the reference sequence(s)
my %fullseq;
open (my $in, "<", $reffile) or die "ERROR: Can't open $reffile: $!\n";
my $sid;
while (my $line = <$in>){
    chomp $line;
    $line =~ s/\s.*$//;
    if ($line =~ m/^>/){
        $sid = substr($line, 1);
        next;
    }
    $fullseq{$sid} .= $line;
}
close ($in);

##if parent assembly file was given, read it in as a string
#my $assem_string = "X";
my %assem_flanks;
#if ($assem){
#    $assem_string = "";
#    open (my $in, "<", $assem) or die "Can't open $assem: $!\n";
#    while (my $line = <$in>){
#        chomp $line;
#        $line =~ s/\s.*$//;
#        if ($line =~ m/^>/){
#            $assem_string .= "X"; ## Separate sequences with X to avoid identifying TA sequences at ends
#            next;
#        }
#        $assem_string .= uc($line);
#    }
#    $assem_string .= "X";
#    close ($in);
#    
#    #Find all TA sites in the assembly sequence file, load hash of TA site flanking sequences
#    my @tas = $assem_string =~ m/(?=([ACTG]{22}TA[ACTG]{22}))/g; #including negative lookahead (?=) to find overlapping matches
#    my $num = scalar @tas;
#    print "Potential transposon insertion sites (TA) in parent strain sequence: $num\n";
#    foreach (@tas){
#        my $full = $_;
#        my $left = substr($full, 0, 24);
#        my $right = substr($full, 24);
#        $right = reverse($right);
#        $right =~ tr/ACTG/TGAC/;
#        push @{$assem_flanks{$left}}, $right;
#        push @{$assem_flanks{$right}}, $left;
#    }
#}

my %refs;
my %align_stats;
my %all_st;
my %all_corr;

#multithread bowtie alignments
#output reads, run bowtie, process results
my ($aout);
open (my $out, ">temp_reads.fasta") or die "ERROR: Can't open output file 'temp_reads.fasta': $!\n";
foreach my $seq (sort keys %sorted){
    #my $count = $sorted{$seq}{'count'};
    #my $leng = $sorted{$seq}{'leng'};
    my $serial = $sorted{$seq}{'serialnumber'};
    #print $out ">R:$leng:$count\n$seq\n";
    print $out ">$serial\n$seq\n";
}
close ($out);
$aout = "temp.bowtie";
print STDERR "\tRunning bowtie  ...\n";
#could just pipe results of bowtie directly into processing, but for the purposes of keeping a trail of steps, will output alignments, then process
$command = "$bow_loc -m 1 --best --strata -a --fullref -n $maxrmm -l 25 temp_ispp -f temp_reads.fasta -p $threads > $aout 2>/dev/null";
$status = system($command);
die "ERROR: command '$command' exited with status $status\n" if $status;
#process the bowtie alignment output
open ($in, "<$aout") or die "ERROR: Can't open bowtie alignment file '$aout': $!\n";
my %sites;
my %potential_mm; #This will be to track and correct those sites where a side has perfect 16bp reads as well as 17bp reads where the first base is mismatched.
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
        my $seq_rc;
        if ($assem){
            $seq_rc = reverse($seq);
            $seq_rc =~ tr/ACTG/TGAC/;
        }
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
        #my @corrected_site_counts; #for keeping track of how many sites end up with modified output counts
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
            
            next if $tot < $minrd; #skip any sites with less than 3 (or user input) reads
            $site_gte_count++;
            
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
print "pool\ttotal_reads\treads_w_tn\treads_w_tn_perfect\treads_correct_length\treads_aligned_to_refs\treads_aligned_TA\treads_aligned_TA_perfect";
my @refids;
foreach my $ref (sort keys %refs){
    push @refids, $ref;
    print "\t$ref\_sites\t$ref\_sites_gte_$minrd";
}
print "\n";
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
    print "$pool_id\t" . join("\t", @results[0..4]) . "\t$reads_aligned_ta\t$results[5]\t" . join("\t", @results[8 .. $#results]) . "\n" ;
    
}

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
sub shift_assess{
    my $side1_perf = shift;
    my $side2_perf = shift;
    my $side2_shift = shift;
    if (abs($side2_perf - $side1_perf) > abs(($side2_perf + $side2_shift) - $side1_perf)){
        return ($side1_perf, ($side2_perf + $side2_shift));
    } else {
        return ($side1_perf, $side2_perf);
    }
}

sub flank_assess{
    my $maj_count = shift;
    my $maj_seq = shift;
    my $min_count = shift;
    my $min_seq = shift;
    my $bar = shift;
    my $ref = shift;
    my $type = shift;
    my $debug = shift;
    
    #my $shortmaj = substr($maj_seq, 0, -2);
    #my @r_flanks = $assem_string =~ /$shortmaj(TA[ACGT]{15})/g;
    
    #my $rev_shortmaj = reverse($shortmaj);
    #$rev_shortmaj =~ tr/ACTG/TGAC/;
    #my @flanks = $assem_string =~ /([ACTG]{15}TA)$rev_shortmaj/g;
    
    #foreach (@r_flanks){
    #    my $rc = reverse($_);
    #    $rc =~ tr/ACTG/TGAC/;
    #    push @flanks, $rc;
    #}
    
    my @flanks;
    @flanks = @{$assem_flanks{$maj_seq}} if $assem_flanks{$maj_seq};
    
    my ($flank_seq, $flank_seq_16);
    if (scalar @flanks <= 1){ #major sequence found only once in contigs (or not at all likely due to undersequencing, misassembly, or being at contig ends, or no assembly sequence given...)
        print STDERR "\t0 - 1 flanks found\n" if $debug;
        if (scalar @flanks == 1){
            $flank_seq = $flanks[0];
            $flank_seq_16 = substr($flank_seq, 1);
        } else {
            $flank_seq = $min_seq;
            $flank_seq_16 = substr($flank_seq, 1);
            $flank_seq_16 = $min_seq if (length($min_seq) == 16);
        }
    } elsif (scalar @flanks > 1) { #major sequence found more than once in contigs
        print STDERR "\t> 1 flanks found\n" if $debug;
        my %all_flanks;
        foreach (@flanks){
            $all_flanks{$_}++;
        }
        if (scalar (keys %all_flanks) == 1){ #IF there is only 1 flanking sequence at all sites:
            print STDERR "\tonly one flanking sequence at all sites, continuing\n" if $debug;
            $flank_seq = $flanks[0];
            $flank_seq_16 = substr($flank_seq, 1);
        } else { #IF there is more than 1 flanking sequence:
            #Output 2x the read count of the lower (non-majority) perfect side
            print STDERR "\tmore than one flanking sequence at all sites: return($min_count, $min_count)\n" if $debug;
            return($min_count, $min_count);
        }
    }
    my $flank_test = $flank_seq;
    $flank_test = $flank_seq_16 if length($min_seq) == 16;
    #Does the sequence of the flank of this site match aligned reads?
    if ($flank_test eq $min_seq){ #if yes...
        print STDERR "\tflank ($flank_test) = min_seq ($min_seq)\n" if $debug;
        my $align_count = 0;
        $align_count += $sorted{$bar}{$flank_seq}{'count'} if ($sorted{$bar}{$flank_seq} and $sorted{$bar}{$flank_seq}{'count'});
        $align_count += $sorted{$bar}{$flank_seq_16}{'count'} if ($sorted{$bar}{$flank_seq_16} and $sorted{$bar}{$flank_seq_16}{'count'});
        #Do the the number of 16bp and 17bp reads aligning to this flank match the total perfect reads for the side?
        if (abs($min_count - $align_count) < 2){ #if yes (with 1 read slop)...
            print STDERR "\talign_count ($align_count) = min_count ($min_count)\n" if $debug;
            print STDERR "\treturn($maj_count, $min_count)\n" if $debug;
            return($maj_count, $min_count);
        } else { #if no...
             print STDERR "\talign_count ($align_count) != min_count ($min_count)\n" if $debug;
            #Look for the number of occurences of the 16bp flank sequence in assembly
            #Get flank sequences at each of these sites that don't match the reference (majority side) sequence
            my $alt_maj_count = 0;
            my $alt_maj_sites = 0; #debug
            
            #my $fshort = substr($flank_seq_16, 0, -2);
            #my @r_flank_flanks = $assem_string =~ /$fshort(TA[ACGT]{15})/g;
            #@r_flank_flanks = $fullseq{$ref} =~ /$fshort(TA[ACGT]{15})/g if !$assem;
            #foreach (@r_flank_flanks){
            #    my $seq = reverse($_);
            #    $seq =~ tr/ACTG/TGAC/;
            #    my $seq16 = substr($seq, 1);
            #    next if $seq16 eq $maj_seq;
            #    next if $seq eq $maj_seq;
            #    $alt_maj_sites ++;
            #    $alt_maj_count += $sorted{$bar}{$seq}{'count'} if $sorted{$bar}{$seq}{'count'};
            #    $alt_maj_count += $sorted{$bar}{$seq16}{'count'} if $sorted{$bar}{$seq16}{'count'};
            #}
            #
            #my $rev_fshort = reverse($fshort);
            #$rev_fshort =~ tr/ACTG/TGAC/;
            #my @flank_flanks = $assem_string =~ /([ACGT]{15}TA)$rev_fshort/g;
            #@flank_flanks = $fullseq{$ref} =~ /([ACGT]{15}TA)$rev_fshort/g if !$assem;
            #foreach (@flank_flanks){
            #    my $seq = $_;
            #    my $seq16 = substr($seq, 1);
            #    next if $seq16 eq $maj_seq;
            #    next if $seq eq $maj_seq;
            #    $alt_maj_sites ++;
            #    $alt_maj_count += $sorted{$bar}{$seq}{'count'} if $sorted{$bar}{$seq}{'count'};
            #    $alt_maj_count += $sorted{$bar}{$seq16}{'count'} if $sorted{$bar}{$seq16}{'count'};
            #}
            
            my @flank_flanks;
            @flank_flanks = @{$assem_flanks{$flank_seq_16}} if $assem_flanks{$flank_seq_16};
            foreach (@flank_flanks){
                my $seq = $_;
                my $seq16 = substr($seq, 1);
                next if $seq16 eq $maj_seq;
                next if $seq eq $maj_seq;
                $alt_maj_sites ++;
                $alt_maj_count += $sorted{$bar}{$seq}{'count'} if ($sorted{$bar}{$seq} and $sorted{$bar}{$seq}{'count'});
                $alt_maj_count += $sorted{$bar}{$seq16}{'count'} if ($sorted{$bar}{$seq16} and $sorted{$bar}{$seq16}{'count'});
            }
            
            print STDERR "\t\talt_maj_sites:$alt_maj_sites, alt_maj_count:$alt_maj_count\n" if $debug;
            if ($alt_maj_count < 2) { #IF no reads aligning to alternative flank sequences
                #Add 16bp flank read count to perfect on that side and output perfect for both sides:
                print STDERR "\t\t1:return($maj_count, $align_count)\n" if $debug;
                return($maj_count, $align_count);
            } else { #IF (more than 1) reads align to alternative flank sequences:
                #Output twice the majority side perfect count
                print STDERR "\t\t2:return($maj_count, $maj_count)\n" if $debug;
                return($maj_count, $maj_count);
            }
        }
    } else { #IF the sequence of the flank of this site does NOT match aligned reads...
        print STDERR "\tflank ($flank_test) != min_seq ($min_seq)\n" if $debug;
        return($maj_count, $min_count) if $type == 7; #for type 7 (complete one-sided imperfect) sites, we want at least one side to be perfect. Otherwise we can't be confident we're aligning reads to the same site as the reference.
        #Get the true flank sequence from assembly
        #Count the number of 16bp and 17bp reads with that sequence
        my $true_flank_count = 0;
        $true_flank_count += $sorted{$bar}{$flank_seq}{'count'} if ($sorted{$bar}{$flank_seq} and $sorted{$bar}{$flank_seq}{'count'});
        $true_flank_count += $sorted{$bar}{$flank_seq_16}{'count'} if ($sorted{$bar}{$flank_seq_16} and $sorted{$bar}{$flank_seq_16}{'count'});
        print STDERR "\treturn($maj_count, $true_flank_count)\n" if $debug;
        return($maj_count, $true_flank_count);
        #### should I go through the process above and make sure the flank sequences are not found in multiple locations? Juice worth the squeeze?
    } 
}
