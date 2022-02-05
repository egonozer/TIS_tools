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
1) sort by barcode and remove barcode sequences
2) identify reads with appropriate transposon sequence and remove transposon
   sequence
3) align reads to the reference genome
4) filter alignments

Prerequisites:
bowtie (version 0.12.8 or higher)

Required:
  -r    Read sequence file, in fastq format. Can be gzipped.
        AND
  -c    Barcode file in following format: 
        pool_id<tab>barcode
            To skip processing of raw reads +/- bowtie alignment, you can also
            include paths to either the processed read output files +/- paths to
            the bowtie alignment files in this file.
            To include the sorted reads file (as previously output by this
            program):
            pool_id<tab>barcode<tab>/path/to/reads.fasta
            To also include the bowtie alignment (as previously output by this
            program):
            pool_id<tab>barcode<tab>/path/to/reads.fasta<tab>/path/to/alignment.bowtie
        AND
  -g    Reference genome sequence file, in fasta format. For best results, this
        file must include ALL sequences present in the organism (i.e. chromsomes
        and plasmids)

Encouraged, but optional:
  -w    assembled sequence of parent strain that was used to produce transposon
        mutagenesis library, in fasta format

Options:
  -t    Transposon sequence
        (default: 'ACAGGTTG')
  -m    Number of transposon sequence mismatches allowed
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
  -d    Minimum difference in read counts between sides to trigger
        investigation of site flanks. A higher number may result in faster
        processing. Minimum value is 1.
        (default: 1)
  -h    Ignore shifted reads
        (default: will add shifted reads to a side's total read count if
        addition of the reads balances the left and right sides of a site)
  -b    Path to directory containing bowtie and bowtie-build. For example:
        /Users/myname/applications/bowtie_folder
        (default: assumes this directory is in your PATH)
  -p    Number of threads
        (default: 1)

";

use Getopt::Std;
use vars qw( $opt_r $opt_c $opt_g $opt_t $opt_m $opt_o $opt_w $opt_s $opt_h $opt_d $opt_b $opt_p );
getopts('r:c:g:t:m:o:w:s:hd:b:p:');

die $usage unless ($opt_r and $opt_c and $opt_g);

my $reads   = $opt_r if $opt_r;
my $bclist  = $opt_c if $opt_c;
my $reffile = $opt_g if $opt_g;

my $tn      = $opt_t ? $opt_t : "ACAGGTTG";
my $maxmm   = defined $opt_m ? $opt_m : 1; #added 'defined' so users can enter 0
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
    print STDERR "WARNING: Output format '$outform' does not match one of the possible options (inseq/essen/wiggl). Defaulting to 'inseq'\n";
    $outform = "inseq";
}

#process the transposon sequence
$tn = uc($tn); 
my @tnarray = split(//, $tn);
my $tnleng = length($tn);

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

## Read in barcodes (+/- reads +/- alingments)
open (my $bin, "<", $bclist) or die "ERROR: Can't open $bclist: $!\n";
my %bc;
my @bc_list;
my $proc_reads;
my %has_reads;
while (my $fline = <$bin>){
    $fline =~ s/\R/\012/g; #converts to UNIX-style line endings
    my @lines = split("\n", $fline); #need to split lines by line-ending character in the case of Mac-formatted files which only have CR line terminators, not both CR and LF like DOS
    foreach my $line (@lines){
        next if $line =~ m/^\s*$/; #skip blank lines
        my @larray = split("\t", $line);
        my ($id, $bar) = @larray;
        $bc{$bar} = $id;
        push @bc_list, [@larray];
        $proc_reads = 1 unless $larray[2];
        $has_reads{$bar}++ if $larray[2];
    }
}
close ($bin);

## Process reads (unless read file was given for every pool)
print STDERR "Starting read processing...\n";
my %sorted;
if ($proc_reads){
    my $rin;
    if ($reads =~ m/\.gz$/){
        open ($rin, "gzip -cd $reads | ") or die "ERROR: Can't open $reads: $!\n";
    } else {
        open ($rin, "<", $reads) or die "ERROR: Can't open $reads: $!\n";
    }
    my $total_reads = 0;
    my %reads_w_bc;
    my %reads_w_tn;
    my ($tot_bc, $tot_tn) = (0) x 2;
    my %tn_mms;
    print STDERR "\tProcessing read $total_reads ... ";
    while (<$rin>){
        my $x1 = $_;
        chomp (my $seq = <$rin>);
        my $x2 = <$rin>;
        my $x3 = <$rin>;
        $total_reads++;
        print STDERR "\r\tProcessing read $total_reads ... " if ($total_reads % 100000 == 0);
        my $first4 = substr($seq, 0, 4, ""); #automatically strip the first 4 bases so I don't have to do another substr below (substr benchmarked faster than regex).
        next if $has_reads{$first4}; #skip saving these reads if a read file was given for this barcode
        if ($bc{$first4}){
            $reads_w_bc{$first4}++;
            $tot_bc++;
            #my $subseq = substr($seq, 4);
            my $tn_match;
            #look for perfect transposon sequence matches. This is fastest
            #for my $i (16 .. 17){
            #    #my $tseq = substr($subseq, $i, $tnleng);
            #    my $tseq = substr($seq, $i, $tnleng);
            #    if ($tseq eq $tn){
            #        $reads_w_tn{$first4}++;
            #        $tot_tn++;
            #        $tn_mms{0}++;
            #        $tn_match = $i;
            #        last;
            #    }
            #}
            #maybe this is faster since it tests for both 16 and 17 simultaneously?
            if ($seq =~ m/$tn/){
                if ($-[0] == 16 or $-[0] == 17){
                    $tn_match = $-[0];
                    $reads_w_tn{$first4}++;
                    $tot_tn++;
                    $tn_mms{0}++;
                }
            }
            
            #if a perfect transposon match was not found, look for imperfect transposon mismatch. This is a slower process and usually less common, so will try not to go through this every time.
            if (!$tn_match and $maxmm > 0){
                for my $i (16 .. 17){
                    #my $tseq = substr($subseq, $i, $tnleng);
                    my $tseq = substr($seq, $i, $tnleng);
                    
                    ### old way
                    #my $num_mm = 0;
                    #my @tseqa = split(//, $tseq);
                    #for my $j (0 .. $#tnarray){
                    #    $num_mm++ unless ($tnarray[$j] eq $tseqa[$j]);
                    #}
                    
                    ## new way (xor, should be much faster)
                    my $num_mm = ($tseq ^ $tn) =~ tr/\0//c;
                    
                    next unless $num_mm <= $maxmm;
                    $reads_w_tn{$first4}++;
                    $tot_tn++;
                    $tn_mms{$num_mm}++;
                    $tn_match = $i;
                    last;
                }
            }
            next unless $tn_match;
            #my $finseq = substr($subseq, 0, $tn_match);
            my $finseq = substr($seq, 0, $tn_match);
            $sorted{$first4}{$finseq}{'count'}++;
            $sorted{$first4}{$finseq}{'leng'} = $tn_match;
        }
    }
    close ($rin);
    print STDERR "\r\tProcessing read $total_reads ... Done!\n";
    print "\nTotal reads: $total_reads\nTotal reads with barcode: $tot_bc\nTotal reads with transposon: $tot_tn\n";
    for my $i (0 .. $maxmm){
        my $val = 0;
        $val = $tn_mms{$i} if $tn_mms{$i};
        print "\tTransposons with $i mismatch(es): $val\n" if $maxmm > 0;
    }
    print "pool\tbc\t#_w_bc\t%_w_bc\t#_w_tn\t%_w_tn\n";
    foreach my $slice (@bc_list){
        my ($id, $bar) = @{$slice};
        next if $has_reads{$bar};
        print "$id\t$bar";
        my ($num_bc, $pct_bc, $num_tn, $pct_tn) = (0) x 4;
        if ($reads_w_bc{$bar}){
            $num_bc = $reads_w_bc{$bar};
            $pct_bc = sprintf("%.1f", 100*($num_bc / $tot_bc));
        }
        if ($reads_w_tn{$bar}){
            $num_tn = $reads_w_tn{$bar};
            $pct_tn = sprintf("%.1f", 100*($num_tn / $num_bc));
        }
        print "\t$num_bc\t$pct_bc\t$num_tn\t$pct_tn\n";
    }
}

## input read files if any were given
foreach (@bc_list){
    my ($id, $bar, $rfile) = @{$_};
    next unless $rfile;
    open (my $in, "<$rfile") or die "ERROR: Can't open read file '$rfile': $!\n";
    my $count;
    my $tot_count = 0;
    while (my $line = <$in>){
        chomp $line;
        if ($line =~ m/^>/){
            $count = (split(":", $line))[2];
            next;
        }
        $line =~ s/\s//g;
        $sorted{$bar}{$line}{'count'} = $count;
        $sorted{$bar}{$line}{'leng'} = length($line);
        $tot_count += $count;
    }
    my $tot_seqs = scalar keys %{$sorted{$bar}};
    print "$id, $bar: Read in $tot_seqs unique read sequences ($tot_count total reads) from file.\n";
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

#if parent assembly file was given, read it in as a string
my $assem_string = "X";
my %assem_flanks;
if ($assem){
    $assem_string = "";
    open (my $in, "<", $assem) or die "Can't open $assem: $!\n";
    while (my $line = <$in>){
        chomp $line;
        $line =~ s/\s.*$//;
        if ($line =~ m/^>/){
            $assem_string .= "X";
            next;
        }
        $assem_string .= uc($line);
    }
    $assem_string .= "X";
    close ($in);
    
    #Find all TA sites in the assembly sequence file, load hash of TA site flanking sequences
    my @tas = $assem_string =~ m/(?=([ACTG]{15}TA[ACTG]{15}))/g; #including negative lookahead (?=) to find overlapping matches
    #my @tas = $assem_string =~ m/([ACTGNX]{15}TA[ACTGNX]{15})/g;
    my $num = scalar @tas;
    print "Potential transposon insertion sites (TA) in parent strain sequence: $num\n";
    foreach (@tas){
        my $full = $_;
        my $left = substr($full, 0, 17);
        my $left_16 = substr($left, 1);
        my $right = substr($full, 15);
        $right = reverse($right);
        $right =~ tr/ACTG/TGAC/;
        my $right_16 = substr($right, 1);
        push @{$assem_flanks{$left}}, $right;
        push @{$assem_flanks{$left_16}}, $right;
        push @{$assem_flanks{$right}}, $left;
        push @{$assem_flanks{$right_16}}, $left;
    }
}

my %refs;
my %align_stats;
my %all_st;
my %all_corr;

#multithread bowtie alignments

#output reads, run bowtie, process results
foreach my $slice (@bc_list){
    my ($id, $bar, $readf, $alignf) = @{$slice};
    if ($sorted{$bar}){
        my ($rout, $aout);
        if ($readf){
            $rout = $readf;
        } else {
            $rout = "$id\_$bar.reads.fasta";
            open (my $out, ">$rout") or die "ERROR: Can't open output file '$rout': $!\n";
            foreach my $seq (sort keys %{$sorted{$bar}}){
                my $count = $sorted{$bar}{$seq}{'count'};
                my $leng = $sorted{$bar}{$seq}{'leng'};
                print $out ">$bar:$leng:$count\n$seq\n";
            }
            close ($out);
        }
        if ($alignf){
            $aout = $alignf;
        } else {
            $aout = "$id\_$bar.bowtie";
            print STDERR "\tRunning bowtie on $id ($bar) ...\n";
            #could just pipe results of bowtie directly into processing, but for the purposes of keeping a trail of steps, will output alignments, then process
            my $command = "$bow_loc -m 1 --best --strata -a --fullref -n 1 -l 17 temp_ispp -f $rout -p $threads > $aout 2>/dev/null";
            my $status = system($command);
            die "ERROR: command '$command' exited with status $status\n" if $status;
        }
        #process the bowtie alignment output
        open (my $in, "<$aout") or die "ERROR: Can't open bowtie alignment file '$aout': $!\n";
        my %sites;
        my %potential_mm; #This will be to track and correct those sites where a side has perfect 16bp reads as well as 17bp reads where the first base is mismatched.
        my %read_seqs; #This will be to store actual aligning read sequences at each site, if comparison to an assembly sequence is done.
        
        #following are for end-of-run stats
        my $num_aligned = 0;
        my $num_aligned_to_TA = 0;
        my $num_aligned_perfect = 0;
        my $num_aligned_mismatch = 0;
        my $num_aligned_shifted = 0;
        
        #parse the bowtie alignment and assign read counts to positions and alignment types.
        print STDERR "\tParsing bowtie alignments ...\n";
        while (my $line = <$in>){
            chomp $line;
            my ($id, $dir, $ref, $pos, $seq, $qual, $x1, $mm) = split("\t", $line);
            my ($bar, $leng, $count) = split(":", $id);
            $num_aligned += $count;
            if ($dir eq "+"){
                my $side = "left";
                my $site = $pos + $leng - 1;
                my $tail = substr($seq, -2);
                if ($tail eq "TA"){
                    $num_aligned_to_TA += $count;
                    if ($mm){
                        $sites{$ref}{$site}{$side}{'imperf'}{'unshift'} += $count;
                        $num_aligned_mismatch += $count;
                        my $mpos;
                        if ($leng == 17){
                            $mpos = (split(":", $mm))[0];
                            $potential_mm{$ref}{$site}{$side}{17} += $count if $mpos == 0;
                        }
                        if ($assem){
                            $mpos = (split(":", $mm))[0] if !$mpos;
                            push @{$read_seqs{$ref}{$site}{$side}}, ([$count, $leng, $seq, $mpos]);
                        }
                    } else {
                        $sites{$ref}{$site}{$side}{'perf'}{'unshift'} += $count;
                        $num_aligned_perfect += $count;
                        $potential_mm{$ref}{$site}{$side}{16} += $count if $leng == 16;
                        push @{$read_seqs{$ref}{$site}{$side}}, ([$count, $leng, $seq, "X"]) if $assem;
                    }
                } else {
                    if ($tail =~ m/T$/ and !$mm){ #if the sequence is otherwise perfect, but apparently missing the terminal A, will shift the site and add it to the count
                        $site++;
                        $sites{$ref}{$site}{$side}{'perf'}{'shift'} += $count;
                        $num_aligned_shifted += $count;
                        push @{$read_seqs{$ref}{$site}{$side}}, ([$count, $leng, $seq, "S"]) if $assem;
                    } elsif ($mm){ #if the mismatch is in the last 2 bases, but either the T or the A is there at the end...
                        #do I really want to save these? Will I really have any use for these reads?
                        my $end1 = $leng - 1 . ":A";
                        my $end2 = $leng - 2 . ":T";
                        if (($tail =~ m/T./ and $mm =~ m/^$end1/) or ($tail =~ m/.A/ and $mm =~ m/^$end2/)){
                            $sites{$ref}{$site}{$side}{'imperf'}{'unshift'} += $count;
                            $num_aligned_to_TA += $count;
                            $num_aligned_mismatch += $count;
                            if ($assem){
                                my $mpos = (split(":", $mm))[0];
                                push @{$read_seqs{$ref}{$site}{$side}}, ([$count, $leng, $seq, $mpos]);
                            }
                        }
                    }
                }
            } else {
                my $side = "right";
                my $site = $pos + 1;
                my $head = substr($seq, 0, 2);
                my $seq_rc;
                if ($assem){
                    $seq_rc = reverse($seq);
                    $seq_rc =~ tr/ACTG/TGAC/;
                }
                if ($head eq "TA"){
                    $num_aligned_to_TA += $count;
                    if ($mm){
                        $sites{$ref}{$site}{$side}{'imperf'}{'unshift'} += $count;
                        $num_aligned_mismatch += $count;
                        my $mpos;
                        if ($leng == 17){
                            $mpos = (split(":", $mm))[0];
                            $potential_mm{$ref}{$site}{$side}{17} += $count if $mpos == 0;
                        }
                        if ($assem){
                            $mpos = (split(":", $mm))[0] if !$mpos;
                            push @{$read_seqs{$ref}{$site}{$side}}, ([$count, $leng, $seq_rc, $mpos]);
                        }
                    } else {
                        $sites{$ref}{$site}{$side}{'perf'}{'unshift'} += $count;
                        $num_aligned_perfect += $count;
                        $potential_mm{$ref}{$site}{$side}{16} += $count if $leng == 16;
                        push @{$read_seqs{$ref}{$site}{$side}}, ([$count, $leng, $seq_rc, "X"]) if $assem;
                    }
                } else {
                    if ($head =~ m/^A/ and !$mm){ #if the sequence is otherwise perfect, but apparently missing the beginning T, will shift the site and add it to the count
                        $site--;
                        $sites{$ref}{$site}{$side}{'perf'}{'shift'} += $count;
                        $num_aligned_shifted += $count;
                        push @{$read_seqs{$ref}{$site}{$side}}, ([$count, $leng, $seq_rc, "S"]) if $assem;
                    } elsif ($mm){ #if the mismatch is in the last 2 bases, but either the T or the A is there at the end...
                        #do I really want to save these? Will I really have any use for these reads?
                        my $end1 = $leng - 1 . ":T";
                        my $end2 = $leng - 2 . ":A";
                        if (($head =~ m/T./ and $mm =~ m/^$end2/) or ($head =~ m/.A/ and $mm =~ m/^$end1/)){
                            $sites{$ref}{$site}{$side}{'imperf'}{'unshift'} += $count;
                            $num_aligned_to_TA += $count;
                            $num_aligned_mismatch += $count;
                            if ($assem){
                                my $mpos = (split(":", $mm))[0];
                                push @{$read_seqs{$ref}{$site}{$side}}, ([$count, $leng, $seq, $mpos]);
                            }
                        }
                    }
                }
            }
        }
        close ($in);
        print "Total reads aligned: $num_aligned\n";
        print "Total reads aligned to TA sites: $num_aligned_to_TA\n";
        print "Total reads perfectly aligned to TA sites: $num_aligned_perfect\n";
        print "Total reads aligned to TA sites with mismatch: $num_aligned_mismatch\n";
        print "Total reads shifted to TA sites: $num_aligned_shifted\n\n";
        
        $align_stats{$bar} = ([$num_aligned, $num_aligned_to_TA, $num_aligned_perfect, $num_aligned_mismatch, $num_aligned_shifted]);
        
        #process and correct transposon insertion site counts
        foreach my $ref (sort keys %sites){
            $refs{$ref}++;
            
            my $outfile = "$id\_$bar.$ref.wiggle";
            $outfile = "allta_split_$id\_$bar.$ref.counts.txt" if $outform eq "essen";
            $outfile = "INSEQ_experiment.scarf_$id\_$bar.bowtiemapped_processed.txt_$ref" if $outform eq "inseq";
            open (my $wig_out, "> $outfile") or die "ERROR: Can't open output file '$outfile': $!\n";
            
            #for status update
            my $num_sites = keys %{$sites{$ref}};
            my $site_count = 0;
            my $site_gte_count = 0;
            print STDERR "Processing site #$site_count (of $num_sites)";
            
            my @site_type_counts;
            my @corrected_site_counts; #for keeping track of how many sites end up with modified output counts
            #set each base value to 0;
            for my $i (1 .. 8){
                $site_type_counts[$i] = 0;
                $corrected_site_counts[$i] = 0;
            }
            my $beg_mm_correct_count = 0; #for keeping track of how many sites are corrected for a beginning mismatch
            
            foreach my $site (sort {$a <=> $b} keys %{$sites{$ref}}){
                $site_count++;
                print STDERR "\rProcessing site #$site_count (of $num_sites)" if $site_count % 100 == 0;
                
                my ($l_perf_unshift, $l_perf_shift, $l_imperf_unshift) = (0) x 3;
                my ($r_perf_unshift, $r_perf_shift, $r_imperf_unshift) = (0) x 3;
                $l_perf_unshift = $sites{$ref}{$site}{'left'}{'perf'}{'unshift'} if $sites{$ref}{$site}{'left'}{'perf'}{'unshift'};
                $l_perf_shift = $sites{$ref}{$site}{'left'}{'perf'}{'shift'} if $sites{$ref}{$site}{'left'}{'perf'}{'shift'};
                $l_imperf_unshift = $sites{$ref}{$site}{'left'}{'imperf'}{'unshift'} if $sites{$ref}{$site}{'left'}{'imperf'}{'unshift'};
                $r_perf_unshift = $sites{$ref}{$site}{'right'}{'perf'}{'unshift'} if $sites{$ref}{$site}{'right'}{'perf'}{'unshift'};
                $r_perf_shift = $sites{$ref}{$site}{'right'}{'perf'}{'shift'} if $sites{$ref}{$site}{'right'}{'perf'}{'shift'};
                $r_imperf_unshift = $sites{$ref}{$site}{'right'}{'imperf'}{'unshift'} if $sites{$ref}{$site}{'right'}{'imperf'}{'unshift'};
                my $l_imp_tot = $l_perf_shift + $l_imperf_unshift;
                my $r_imp_tot = $r_perf_shift + $r_imperf_unshift;
                my $tot = $l_perf_unshift + $l_imp_tot + $r_perf_unshift + $r_imp_tot;
                
                next if $tot < $minrd; #skip any sites with less than 3 (or user input) reads
                $site_gte_count++;
                
                #correct counts at cases of read beginning mismatches
                my %r_beg_mm_sites;
                foreach my $side (qw(left right)){
                    my ($perf_unshift, $imperf_unshift) = ($l_perf_unshift, $l_imperf_unshift);
                    ($perf_unshift, $imperf_unshift) = ($r_perf_unshift, $r_imperf_unshift) if $side eq "right";
                    if ($potential_mm{$ref}{$site}{$side}{16} and $potential_mm{$ref}{$site}{$side}{17}){
                        my $sixteen_count = $potential_mm{$ref}{$site}{$side}{16};
                        my $seventeen_count = $potential_mm{$ref}{$site}{$side}{17};
                        if ($seventeen_count / $imperf_unshift > 0.8 and $seventeen_count > 2 and $perf_unshift - $sixteen_count < 3){ #if the majority (more than 80%) of the total mismatch reads have the mismatch at the first base and there are 3 or more of these reads and all or almost all the perfect reads were these 16 bp reads
                            $perf_unshift -= $sixteen_count; #subtract the number of perfect 16bp reads from the total perfect read count
                            $imperf_unshift += $sixteen_count; #and add that total to the imperfect read count
                            $beg_mm_correct_count ++;
                            $r_beg_mm_sites{$site}++;
                        } #...otherwise, unable to know how many of those 16 bp reads belong to perfect vs imperfect, so we'll just count them all as perfect since we have a "significant" number of perfect 17bp reads
                    }
                }
                
                #assign site to site type and correct counts if possible
                my $assigned;
                my ($l_out, $r_out) = (0,0);
                if ($l_perf_unshift > 0 and $r_perf_unshift > 0){ #AT LEAST ONE PERFECT ON BOTH SIDES
                    if ($l_imp_tot < $l_perf_unshift and $r_imp_tot < $r_perf_unshift){ #MAJORITY PERFECT ON BOTH SIDES
                        if ($l_imp_tot < 1 and $r_imp_tot < 1){ #NO IMPERFECT ON EITHER SIDE
                            $site_type_counts[1]++;
                            $assigned = 1;
                            ($l_out, $r_out) = ($l_perf_unshift, $r_perf_unshift);
                            if (abs($l_perf_unshift - $r_perf_unshift) >= $maxdiff){ 
                                my $l_seq;
                                foreach (@{$read_seqs{$ref}{$site}{"left"}}){
                                    my ($count, $leng, $seq, $mm) = @{$_};
                                    next unless $mm eq "X";
                                    $l_seq = $seq unless ($l_seq and $leng == 16); #try to get 17bp sequence, but will take 16bp if no 17bp available
                                }
                                my $r_seq;
                                foreach (@{$read_seqs{$ref}{$site}{"right"}}){
                                    my ($count, $leng, $seq, $mm) = @{$_};
                                    next unless $mm eq "X";
                                    $r_seq = $seq unless ($r_seq and $leng == 16); #try to get 17bp sequence, but will take 16bp if no 17bp available
                                }
                                if ($l_seq and $r_seq){
                                    my @order = ([$l_perf_unshift, $l_seq, "L"], [$r_perf_unshift, $r_seq, "R"]);
                                    @order = sort{$b->[0] <=> $a->[0]} @order;
                                    my ($top, $bot) = flank_assess($order[0][0], $order[0][1], $order[1][0], $order[1][1], $bar, $ref, 1, 0);
                                    $corrected_site_counts[$assigned]++ if ($top != $order[0][0] or $bot != $order[1][0]); 
                                    ($l_out, $r_out) = ($top, $bot);
                                    ($l_out, $r_out) = ($bot, $top) if $order[0][2] eq "R";
                                }
                            }
                        } else { #MINORITY IMPERFECT ON ONE OR BOTH SIDES
                            $site_type_counts[2]++;
                            $assigned = 2;
                            ($l_out, $r_out) = ($l_perf_unshift, $r_perf_unshift);
                            if (abs($l_perf_unshift - $r_perf_unshift) >= $maxdiff){
                                my $l_seq;
                                foreach (@{$read_seqs{$ref}{$site}{"left"}}){
                                    my ($count, $leng, $seq, $mm) = @{$_};
                                    next unless $mm eq "X";
                                    $l_seq = $seq unless ($l_seq and $leng == 16); #try to get 17bp sequence, but will take 16bp if no 17bp available
                                }
                                my $r_seq;
                                foreach (@{$read_seqs{$ref}{$site}{"right"}}){
                                    my ($count, $leng, $seq, $mm) = @{$_};
                                    next unless $mm eq "X";
                                    $r_seq = $seq unless ($r_seq and $leng == 16); #try to get 17bp sequence, but will take 16bp if no 17bp available
                                }
                                if ($l_seq and $r_seq){
                                    my @order = ([$l_perf_unshift, $l_seq, "L"], [$r_perf_unshift, $r_seq, "R"]);
                                    @order = sort{$b->[0] <=> $a->[0]} @order;
                                    my ($top, $bot) = flank_assess($order[0][0], $order[0][1], $order[1][0], $order[1][1], $bar, $ref, 2, 0);
                                    ($l_out, $r_out) = ($top, $bot);
                                    ($l_out, $r_out) = ($bot, $top) if $order[0][2] eq "R";
                                }
                            }
                            unless ($opt_h){
                                unless ($l_perf_shift > 1 and $r_perf_shift > 1){
                                    if ($l_perf_shift > 1){
                                        ($r_out, $l_out) = shift_assess($r_out, $l_out, $l_perf_shift);
                                    } elsif ($r_perf_shift > 1){
                                        ($l_out, $r_out) = shift_assess($l_out, $r_out, $r_perf_shift);
                                    }
                                }
                            }
                            $corrected_site_counts[$assigned]++ if ($l_out != $l_perf_unshift or $r_out != $r_perf_unshift);
                        }
                    } else { #MAJORITY IMPERFECT ON ONE OR BOTH SIDES
                        $site_type_counts[3]++;
                        $assigned = 3;
                        if ($assem){
                            my @l_reads = @{$read_seqs{$ref}{$site}{"left"}};
                            my @r_reads = @{$read_seqs{$ref}{$site}{"right"}};
                            my (@test, @query);
                            if ($l_imp_tot < $l_perf_unshift){
                                @test = @l_reads;
                                @query = @r_reads;
                            } elsif ($r_imp_tot < $r_perf_unshift){
                                @test = @r_reads;
                                @query = @l_reads;
                            } else { #if perfect counts on BOTH sides are less than imperfect counts on BOTH sides
                                #don't do anything, just output the perfects on both sides.
                            }
                            my %opps;
                            my $top_perf;
                            foreach my $slice (@test){
                                my ($count, $leng, $seq, $mpos) = @{$slice};
                                if ($mpos eq "X"){
                                    if ($top_perf){
                                        $top_perf = $seq if $leng == 17;
                                    } else {
                                        $top_perf = $seq;
                                    }
                                    #my $shortseq = substr($seq, 0, $leng - 2);
                                    #my $shortseq_rc = reverse($shortseq);
                                    #$shortseq_rc =~ tr/ACTG/TGAC/;
                                    #my @rmatches = $assem_string =~ /$shortseq(TA[ACGT]{15})/g;
                                    #my @fmatches = $assem_string =~ /([ACGT]{15}TA)$shortseq_rc/g;
                                    #foreach (@fmatches){
                                    #    $opps{$_}++;
                                    #    $opps{substr($_, 1)}++;
                                    #}
                                    #foreach (@rmatches){
                                    #    my $rc = reverse($_);
                                    #    $rc =~ tr/ACTG/TGAC/;
                                    #    $opps{$rc}++;
                                    #    $opps{substr($rc, 1)}++;
                                    #}
                                    
                                    my @fmatches;
                                    @fmatches = @{$assem_flanks{$seq}} if $assem_flanks{$seq};
                                    foreach (@fmatches){
                                        $opps{$_}++;
                                        $opps{substr($_, 1)}++;
                                    }
                                }
                            }
                            my ($num_x, $num_mm) = (0) x 2;
                            my ($top_x, $top_mm);
                            foreach my $slice (@query){
                                my ($count, $leng, $seq, $mpos) = @{$slice};
                                next if $mpos eq "S";
                                if ($opps{$seq}){
                                    if ($mpos eq "X"){
                                        if ($r_beg_mm_sites{$site} and $leng == 16){ #if this is a read-beginning-mismatch site, count perfect 16bp reads as imperfect
                                            $num_mm += $count;
                                            $top_mm = $seq if !$top_mm;
                                        } else {
                                            $num_x += $count;
                                            if ($top_x){
                                                $top_x = $seq if $leng == 17
                                            } else {
                                                $top_x = $seq;
                                            }
                                        }
                                    } else {
                                        $num_mm += $count;
                                        if ($top_mm){
                                            $top_mm = $seq if $leng == 17
                                        } else {
                                            $top_mm = $seq;
                                        }
                                    }
                                }
                            }
                            ($l_out, $r_out) = ($l_perf_unshift, $r_perf_unshift);
                            if ($num_x == 0 and $num_mm > 0){ #output perfect on one side and imperfect on the other if there are no perfect reads opposite the majority perfect side
                                if ($l_imp_tot < $l_perf_unshift){
                                    ($l_out, $r_out) = ($l_perf_unshift, $num_mm);
                                    if (abs($l_out - $r_out) >= $maxdiff){ ##assess flanks if read difference is too great
                                        my @order = ([$l_out, $top_perf, "L"], [$r_out, $top_mm, "R"]);
                                        @order = sort{$b->[0] <=> $a->[0]} @order;
                                        my ($top, $bot) = flank_assess($order[0][0], $order[0][1], $order[1][0], $order[1][1], $bar, $ref, 3, 0);
                                        ($l_out, $r_out) = ($top, $bot);
                                        ($l_out, $r_out) = ($bot, $top) if $order[0][2] eq "R";
                                    }
                                    if ($l_perf_shift > 1 and !$opt_h){
                                        ($r_out, $l_out) = shift_assess($num_mm, $l_perf_unshift, $l_perf_shift);
                                    }
                                } elsif ($r_imp_tot < $r_perf_unshift){
                                    ($l_out, $r_out) = ($num_mm, $r_perf_unshift);
                                    if (abs($l_out - $r_out) >= $maxdiff){ ##assess flanks if read difference is too great
                                        my @order = ([$l_out, $top_mm, "L"], [$r_out, $top_perf, "R"]);
                                        @order = sort{$b->[0] <=> $a->[0]} @order;
                                        my ($top, $bot) = flank_assess($order[0][0], $order[0][1], $order[1][0], $order[1][1], $bar, $ref, 3, 0);
                                        ($l_out, $r_out) = ($top, $bot);
                                        ($l_out, $r_out) = ($bot, $top) if $order[0][2] eq "R";
                                    }
                                    if ($r_perf_shift > 1 and !$opt_h){
                                        ($l_out, $r_out) = shift_assess($num_mm, $r_perf_unshift, $r_perf_shift);
                                    }
                                }
                            } elsif ($num_x > 0 and $num_mm > 0) { #output minority perfect x 2 if there are perfect and imperfect opposite the majority perfect sequence
                                if ($l_imp_tot < $l_perf_unshift){
                                    ($l_out, $r_out) = ($r_perf_unshift, $r_perf_unshift);
                                    if ($r_perf_shift > 1 and !$opt_h){
                                        my ($side1, $side2) = shift_assess($l_perf_shift, ($r_perf_unshift + $num_mm), $r_perf_shift); #check to see if shifted reads brings the total read count on one side closer to the perfects on the other
                                        unless ($side2 == $r_perf_unshift + $num_mm){ #add shifted
                                            ($l_out, $r_out) = (($r_perf_unshift + $r_perf_shift), ($r_perf_unshift + $r_perf_shift));
                                        }
                                    }
                                } elsif ($r_imp_tot < $r_perf_unshift){
                                    ($l_out, $r_out) = ($l_perf_unshift, $l_perf_unshift);
                                    if ($l_perf_shift > 1 and !$opt_h){
                                        my ($side1, $side2) = shift_assess($r_perf_shift, ($l_perf_unshift + $num_mm), $l_perf_shift); #check to see if shifted reads brings the total read count on one side closer to the perfects on the other
                                        unless ($side2 == $l_perf_unshift + $num_mm){ #add shifted
                                            ($l_out, $r_out) = (($l_perf_unshift + $l_perf_shift), ($l_perf_unshift + $l_perf_shift));
                                        }
                                    }
                                }
                            } elsif ($num_x == 0 and $num_mm == 0){ #if no aligned reads match the flanking site(s), search for reads that do (likely indels)
                                if (keys %opps == 0){ #if the majority perfect sequence was not found in the assembly sequence (? misassembly)
                                    #output perfect on both sides
                                    ($l_out, $r_out) = ($l_perf_unshift, $r_perf_unshift);
                                    unless ($l_perf_shift > 2 and $r_perf_shift > 2){
                                        if ($l_perf_shift > 1 and !$opt_h){
                                            ($r_out, $l_out) = shift_assess($r_perf_unshift, $l_perf_unshift, $l_perf_shift);
                                        }
                                        if ($r_perf_shift > 1 and !$opt_h){
                                            ($l_out, $r_out) = shift_assess($l_perf_unshift, $r_perf_unshift, $r_perf_shift);
                                        }
                                    }
                                } else {
                                    my (@bp16, @bp17);
                                    foreach my $flank (keys %opps){
                                        my $flankleng = length($flank);
                                        my $val = 0;
                                        $val = $sorted{$bar}{$flank}{'count'} if ($sorted{$bar}{$flank} and $sorted{$bar}{$flank}{'count'});
                                        if ($val > 1){ #allow a 1 read slop
                                            push @bp16, ([$val, $flank]) if $flankleng == 16;
                                            push @bp17, ([$val, $flank]) if $flankleng == 17;
                                        }
                                    }
                                    if (scalar @bp16 <= 1 and scalar @bp17 <= 1){
                                        $num_mm += $bp16[0][0] if $bp16[0];
                                        $num_mm += $bp17[0][0] if $bp17[0];
                                        my $top_flank = $bp16[0][1] if $bp16[0];
                                        $top_flank = $bp17[0][1] if $bp17[0];
                                        if ($l_imp_tot < $l_perf_unshift){
                                            ($l_out, $r_out) = ($l_perf_unshift, $num_mm);
                                            if (abs($l_out - $r_out) >= $maxdiff and $num_mm > 0){ ##assess flanks if read difference is too great
                                                my @order = ([$l_out, $top_perf, "L"], [$r_out, $top_flank, "R"]);
                                                @order = sort{$b->[0] <=> $a->[0]} @order;
                                                my ($top, $bot) = flank_assess($order[0][0], $order[0][1], $order[1][0], $order[1][1], $bar, $ref, 3, 0);
                                                ($l_out, $r_out) = ($top, $bot);
                                                ($l_out, $r_out) = ($bot, $top) if $order[0][2] eq "R";
                                            } 
                                        } elsif ($r_imp_tot < $r_perf_unshift) {
                                            ($l_out, $r_out) = ($num_mm, $r_perf_unshift);
                                            if (abs($l_out - $r_out) >= $maxdiff and $num_mm > 0){ ##assess flanks if read difference is too great
                                                my @order = ([$l_out, $top_flank, "L"], [$r_out, $top_perf, "R"]);
                                                @order = sort{$b->[0] <=> $a->[0]} @order;
                                                my ($top, $bot) = flank_assess($order[0][0], $order[0][1], $order[1][0], $order[1][1], $bar, $ref, 3, 0);
                                                ($l_out, $r_out) = ($top, $bot);
                                                ($l_out, $r_out) = ($bot, $top) if $order[0][2] eq "R";
                                            }
                                        }
                                        unless ($l_perf_shift > 2 and $r_perf_shift > 2){
                                            if ($l_perf_shift > 1 and !$opt_h){
                                                ($r_out, $l_out) = shift_assess($r_out, $l_out, $l_perf_shift);
                                            }
                                            if ($r_perf_shift > 1 and !$opt_h){
                                                ($l_out, $r_out) = shift_assess($l_out, $r_out, $r_perf_shift);
                                            }
                                        }
                                    } else {
                                        #if there are multiple flank sites and reads at the multiple flank sites, won't be clear where the reads belong and we should ignore the site
                                        ($l_out, $r_out) = (0,0);
                                    }
                                }
                            } else {
                                #if there are only perfect reads opposite the majority perfect side, default to outputting only perfect read counts from either side
                                if (abs($l_perf_unshift - $r_perf_unshift) >= $maxdiff){ 
                                    my $l_seq;
                                    foreach (@{$read_seqs{$ref}{$site}{"left"}}){
                                        my ($count, $leng, $seq, $mm) = @{$_};
                                        next unless $mm eq "X";
                                        $l_seq = $seq unless ($l_seq and $leng == 16); #try to get 17bp sequence, but will take 16bp if no 17bp available
                                    }
                                    my $r_seq;
                                    foreach (@{$read_seqs{$ref}{$site}{"right"}}){
                                        my ($count, $leng, $seq, $mm) = @{$_};
                                        next unless $mm eq "X";
                                        $r_seq = $seq unless ($r_seq and $leng == 16); #try to get 17bp sequence, but will take 16bp if no 17bp available
                                    }
                                    if ($l_seq and $r_seq){
                                        my @order = ([$l_perf_unshift, $l_seq, "L"], [$r_perf_unshift, $r_seq, "R"]);
                                        @order = sort{$b->[0] <=> $a->[0]} @order;
                                        my ($top, $bot) = flank_assess($order[0][0], $order[0][1], $order[1][0], $order[1][1], $bar, $ref, 3, 0);
                                        ($l_out, $r_out) = ($top, $bot);
                                        ($l_out, $r_out) = ($bot, $top) if $order[0][2] eq "R";
                                    }
                                }
                                unless ($l_perf_shift > 1 and $r_perf_shift > 1){
                                    if ($l_perf_shift > 1 and !$opt_h){
                                        ($r_out, $l_out) = shift_assess($r_out, $l_out, $l_perf_shift);
                                    } elsif ($r_perf_shift > 1 and !$opt_h){
                                        ($l_out, $r_out) = shift_assess($l_out, $r_out, $r_perf_shift);
                                    }
                                }
                            }
                        } else { #if no assembly is available, will default to outputting double the minority perfect side, which seems safest
                            ($l_out, $r_out) = ($l_perf_unshift, $r_perf_unshift); #default to perfects on both sides if perfects are minority on both sides
                            unless ($l_perf_shift > 2 and $r_perf_shift > 2){ ##### is this OK? should I be checking both sides when both sides have a lot of shifted?
                                if ($l_perf_shift > 1 and !$opt_h){
                                    ($r_out, $l_out) = shift_assess($r_perf_unshift, $l_perf_unshift, $l_perf_shift);
                                } elsif ($r_perf_shift > 1 and !$opt_h){
                                    ($l_out, $r_out) = shift_assess($l_perf_unshift, $r_perf_unshift, $r_perf_shift);
                                }
                            }
                            if ($l_imp_tot < $l_perf_unshift){
                                ($l_out, $r_out) = ($r_perf_unshift, $r_perf_unshift);
                                if ($r_perf_shift > 1 and !$opt_h){
                                    my ($side1, $side2) = shift_assess($l_perf_shift, ($r_perf_unshift + $r_imperf_unshift), $r_perf_shift); #check to see if shifted reads brings the total read count on one side closer to the perfects on the other
                                    unless ($side2 == $r_perf_unshift + $r_imperf_unshift){ #add shifted
                                        ($l_out, $r_out) = (($r_perf_unshift + $r_perf_shift), ($r_perf_unshift + $r_perf_shift));
                                    }
                                }
                            } elsif ($r_imp_tot < $r_perf_unshift){
                                ($l_out, $r_out) = ($l_perf_unshift, $l_perf_unshift);
                                if ($l_perf_shift > 1 and !$opt_h){
                                    my ($side1, $side2) = shift_assess($r_perf_shift, ($l_perf_unshift + $l_imperf_unshift), $l_perf_shift); #check to see if shifted reads brings the total read count on one side closer to the perfects on the other
                                    unless ($side2 == $l_perf_unshift + $l_imperf_unshift){ #add shifted
                                        ($l_out, $r_out) = (($l_perf_unshift + $l_perf_shift), ($l_perf_unshift + $l_perf_shift));
                                    }
                                }
                            }
                        }
                        $corrected_site_counts[$assigned]++ if ($l_out != $l_perf_unshift or $r_out != $r_perf_unshift);
                    }
                } else { #ONE OR BOTH SIDES WITHOUT ANY PERFECT
                    if (($l_perf_unshift > 0 and $r_imp_tot > 0) or ($r_perf_unshift > 0 and $l_imp_tot > 0)){ #PERFECT ON ONE SIDE, IMPERFECT ON THE OTHER
                        $site_type_counts[4]++;
                        $assigned = 4;
                        if ($assem){
                            my @perf_reads = @{$read_seqs{$ref}{$site}{"left"}};
                            @perf_reads = @{$read_seqs{$ref}{$site}{"right"}} if $r_perf_unshift > 1;
                            my %opps;
                            my %s_opps;
                            my @perf_seqs;
                            my $top_perf;
                            foreach my $slice (@perf_reads){
                                my ($count, $leng, $seq, $mpos) = @{$slice};
                                if ($mpos eq "X"){
                                    if ($top_perf){
                                        $top_perf = $seq if $leng == 17;
                                    } else {
                                        $top_perf = $seq;
                                    }
                                    my $shortseq = substr($seq, 0, $leng - 2);
                                    my $shortseq_rc = reverse($shortseq);
                                    $shortseq_rc =~ tr/ACTG/TGAC/;
                                    push @perf_seqs, ([$shortseq, $shortseq_rc, $leng]); #### need to figure out shifted. Do I need them at all?
                                    #my @rmatches = $assem_string =~ /$shortseq(TA[ACGT]{15})/g;
                                    $shortseq .= "T";
                                    my @rshiftmatches = $assem_string =~ /$shortseq(A[ACGT]{16})/g; #look for shifted sequences
                                    #my @fmatches = $assem_string =~ /([ACGT]{15}TA)$shortseq_rc/g;
                                    my $astring = "A";
                                    $astring .= $shortseq_rc;
                                    my @fshiftmatches = $assem_string =~ /([ACGT]{16}T)$astring/g; #look for shifted sequences
                                    #foreach (@fmatches){
                                    #    $opps{$_}{$leng}++;
                                    #    $opps{substr($_, 1)}{$leng}++;
                                    #}
                                    #foreach (@rmatches){
                                    #    my $rc = reverse($_);
                                    #    $rc =~ tr/ACTG/TGAC/;
                                    #    $opps{$rc}{$leng}++;
                                    #    $opps{substr($rc, 1)}{$leng}++;
                                    #}
                                    foreach (@fshiftmatches){
                                        $s_opps{$_}{$leng}++;
                                        $s_opps{substr($_, 1)}{$leng}++;
                                    }
                                    foreach (@rshiftmatches){
                                        my $rc = reverse($_);
                                        $rc =~ tr/ACTG/TGAC/;
                                        $s_opps{$rc}{$leng}++;
                                        $s_opps{substr($rc, 1)}{$leng}++;
                                    }
                                    
                                    my @fmatches;
                                    @fmatches = @{$assem_flanks{$seq}} if $assem_flanks{$seq};
                                    foreach (@fmatches){
                                        $opps{$_}{$leng}++;
                                        $opps{substr($_, 1)}{$leng}++;
                                    }
                                    
                                }
                            }                            
                            my $fs = scalar keys %opps;
                            if ($fs == 0){ #no flanking sequence found, either because the sequence had degenerate bases (Ns) or was near the end of a contig
                                #default to outputting double the perfect side (which seems safest)
                                ($l_out, $r_out) = ($l_perf_unshift, $l_perf_unshift);
                                ($l_out, $r_out) = ($r_perf_unshift, $r_perf_unshift) if ($r_perf_unshift > $l_perf_unshift);
                            } elsif ($fs > 0 and $fs <= 2){ #2 is the minimum (16bp flanking sequence and 17bp flanking sequence)
                                my $fcount = 0;
                                my $fs_count = 0;
                                my $top_flank;
                                foreach my $fseq (keys %opps){
                                    #print STDERR "fseq: $fseq\n";
                                    $fcount += $sorted{$bar}{$fseq}{'count'} if ($sorted{$bar}{$fseq} and $sorted{$bar}{$fseq}{'count'});
                                    if ($top_flank){
                                        $top_flank = $fseq if length($fseq) == 17;
                                    } else {
                                        $top_flank = $fseq;
                                    }
                                }
                                foreach my $sseq (keys %s_opps){
                                    $fs_count += $sorted{$bar}{$sseq}{'count'} if ($sorted{$bar}{$sseq} and $sorted{$bar}{$sseq}{'count'});
                                }
                                if ($l_perf_unshift > 1){
                                    ($l_out, $r_out) = ($l_perf_unshift, $fcount);
                                    if (abs($l_out - $r_out) >= $maxdiff){
                                        my @order = ([$l_perf_unshift, $top_perf, "L"], [$fcount, $top_flank, "R"]);
                                        @order = sort{$b->[0] <=> $a->[0]} @order;
                                        my ($top, $bot) = flank_assess($order[0][0], $order[0][1], $order[1][0], $order[1][1], $bar, $ref, 4, 0); 
                                        ($l_out, $r_out) = ($top, $bot);
                                        ($l_out, $r_out) = ($bot, $top) if $order[0][2] eq "R";
                                        $r_out += $fs_count; #add back the shifted read counts
                                        if ($l_perf_shift > 1 and !$opt_h){
                                            ($r_out, $l_out) = shift_assess($r_out, $l_perf_unshift, $l_perf_shift);
                                        }
                                    }
                                } elsif ($r_perf_unshift > 1){
                                    ($l_out, $r_out) = ($fcount, $r_perf_unshift);
                                    if (abs($l_out - $r_out) >= $maxdiff){
                                        my @order = ([$fcount, $top_flank, "L"], [$r_perf_unshift, $top_perf, "R"]);
                                        @order = sort{$b->[0] <=> $a->[0]} @order;
                                        my ($top, $bot) = flank_assess($order[0][0], $order[0][1], $order[1][0], $order[1][1], $bar, $ref, 4, 0); 
                                        ($l_out, $r_out) = ($top, $bot);
                                        ($l_out, $r_out) = ($bot, $top) if $order[0][2] eq "R";
                                        $l_out += $fs_count; #add back the shifted read counts
                                        if ($r_perf_shift > 1 and !$opt_h){
                                            ($l_out, $r_out) = shift_assess($l_out, $r_perf_unshift, $r_perf_shift);
                                        }
                                    }   
                                }
                            } elsif ($fs > 2) { #If the sequence on the majority perfect side shows up more than once in the assembly contigs, we should look to see if one of those times it is opposite a "perfect" sequence matching the reference.
                                #get the "perfect" sequence of the opposite side from the reference
                                my ($o_perf, $o_perf_shift);
                                foreach my $slice (@perf_seqs){
                                    last if $o_perf;
                                    my ($shortseq, $shortseq_rc, $leng) = @{$slice};
                                    (my $rmatch) = $fullseq{$ref} =~ /$shortseq(TA[ACGT]{15})/;
                                    if ($rmatch){
                                        $o_perf = reverse($rmatch);
                                        $o_perf =~ tr/ACTG/TGAC/;
                                        $shortseq .= "T";
                                        (my $rshiftmatch) = $fullseq{$ref} =~ /$shortseq(A[ACGT]{16})/; #look for shifted sequences
                                        $o_perf_shift = reverse($rshiftmatch);
                                        $o_perf_shift =~ tr/ACTG/TGAC/;
                                    } else{
                                        ($o_perf) = $fullseq{$ref} =~ /([ACGT]{15}TA)$shortseq_rc/;
                                        my $astring = "A";
                                        $astring .= $shortseq_rc;
                                        ($o_perf_shift) = $fullseq{$ref} =~ /([ACGT]{16}T)$astring/; #look for shifted sequences
                                    }
                                }
                                ($l_out, $r_out) = (0) x 2;
                                if ($opps{$o_perf} or $opps{substr($o_perf, 1)}){ #check whether any of the sequences opposite the perfect sequence in the assembly are also "perfect" with regards to the reference
                                    #if so, count the number of actual reads aligning to the perfect sequence
                                    my @query = @{$read_seqs{$ref}{$site}{"right"}};
                                    @query = @{$read_seqs{$ref}{$site}{"left"}} if $r_perf_unshift > $l_perf_unshift;
                                    my $perfcount = 0;
                                    foreach my $slice (@query){
                                        my ($count, $leng, $seq, $mpos) = @{$slice};
                                        next unless $mpos eq "X";
                                        $perfcount += $count;
                                    }
                                    if ($perfcount > 0){ #if reads aligning to the perfect site were found, output double that number
                                        ($l_out, $r_out) = ($perfcount, $perfcount);
                                    }
                                    #if no reads aligning to the perfect site were found, then output 0 reads for the site
                                } else {
                                    #if none of the sequences opposite the perfect sequence were "perfect" with regards to the reference, it will be too difficult to determine which of the "imperfect" sequences is most likely to correspond to the site in the reference
                                    #of the three possible choices (1: output no reads for the site or 2: output the perfect reads on one side and imperfect reads on the other side, or 3: output twice the perfect reads), I think perhaps not outputting any reads is the safest since I can't tell which reads belong where
                                }
                            }
                            
                            ## some code here to deal with shifted reads (%s_opps), or just ignore?
                            
                        } else { #if no assembly is available, default to outputting double the perfect side (which seems safest)
                            ($l_out, $r_out) = ($l_perf_unshift, $l_perf_unshift);
                            ($l_out, $r_out) = ($r_perf_unshift, $r_perf_unshift) if ($r_perf_unshift > $l_perf_unshift);
                        }
                        $corrected_site_counts[$assigned]++ if ($l_out != $l_perf_unshift or $r_out != $r_perf_unshift);  
                    } elsif ($l_imp_tot > 0 and $r_imp_tot > 0){ #IMPERFECT ON BOTH SIDES
                        $site_type_counts[5]++;
                        $assigned = 5;
                        ($l_out, $r_out) = (0,0); #if assembly is not available, don't output these sites. Can't know if they're real and if both sides are at the same site.
                        if ($assem){
                            # take major sequences of both sides, see if they are opposite each other in the assembly, output both if yes.
                            my @left = @{$read_seqs{$ref}{$site}{"left"}};
                            my @right = @{$read_seqs{$ref}{$site}{"right"}};
                            my %l_test_seqs;
                            my %l_test_seq_counts;
                            my %l_max;
                            my %r_test_seqs;
                            my %r_test_seq_counts;
                            my %r_max;
                            foreach my $lshift (@left){
                                my ($lcount, $lleng, $lseq, $lmpos) = @{$lshift};
                                $lseq .= "A" if $lmpos eq "S";
                                my $l_test_seq = $lseq;
                                if ($lmpos eq "S"){
                                    $l_test_seq = substr($lseq, 1); #if lleng is 16
                                    $l_test_seq = substr($lseq, 2) if $lleng == 17; 
                                } else {
                                    $l_test_seq = substr($lseq, 1) if $lleng == 17;
                                }
                                $l_test_seq_counts{$l_test_seq} += $lcount;
                                if ($l_max{$l_test_seq}){
                                    my $val = $l_max{$l_test_seq};
                                    if (length($lseq) > $val){
                                        $l_test_seqs{$l_test_seq} = $lseq;
                                        $l_max{$l_test_seq} = length($lseq) 
                                    }
                                } else {
                                    $l_test_seqs{$l_test_seq} = $lseq;
                                    $l_max{$l_test_seq} = length($lseq);
                                }
                            }
                            foreach my $rshift (@right){
                                my ($rcount, $rleng, $rseq, $rmpos) = @{$rshift};
                                $rseq .= "A" if $rmpos eq "S";
                                my $r_test_seq = $rseq;
                                if ($rmpos eq "S"){
                                    $r_test_seq = substr($rseq, 1); #if rleng is 16
                                    $r_test_seq = substr($rseq, 2) if $rleng == 17; 
                                } else {
                                    $r_test_seq = substr($rseq, 1) if $rleng == 17;
                                }
                                $r_test_seq_counts{$r_test_seq} += $rcount;
                                if ($r_max{$r_test_seq}){
                                    my $val = $r_max{$r_test_seq};
                                    if (length($rseq) > $val){
                                        $r_test_seqs{$r_test_seq} = $rseq;
                                        $r_max{$r_test_seq} = length($rseq) 
                                    }
                                } else {
                                    $r_test_seqs{$r_test_seq} = $rseq;
                                    $r_max{$r_test_seq} = length($rseq);
                                }
                            }
                            my @hits;
                            foreach my $l_test_seq (keys %l_test_seqs){
                                my $l_long_seq = $l_test_seqs{$l_test_seq};
                                foreach my $r_test_seq (keys %r_test_seqs){
                                    my $r_long_seq = $r_test_seqs{$r_test_seq};
                                    my $rc = reverse($r_long_seq);
                                    $rc = substr($rc, 2);
                                    $rc =~ tr/ACTG/TGAC/;
                                    my $full = $l_test_seq . $rc;
                                    my $full_rc = reverse($full);
                                    $full_rc =~ tr/ACTG/TGAC/;
                                    #look for the full sequence in the assembly;
                                    if ($assem_string =~ m/$full|$full_rc/){ # should I count the number of matches?
                                        my $lcount = $l_test_seq_counts{$l_test_seq};
                                        my $rcount = $r_test_seq_counts{$r_test_seq};
                                        my $ratio = ($lcount / $rcount);
                                        $ratio = (1/$ratio) if $ratio < 1;
                                        next if $ratio >= 20; #not sure if the ratio approach is best, but I can't see a better approach to weed out errors and maintain consistency between pools than this right now.
                                        push @hits, ([$lcount, $rcount, $lcount + $rcount]);
                                    }
                                }
                            }
                            if (@hits){
                                @hits = sort{$b->[2] <=> $a->[2]}@hits;
                                ($l_out, $r_out) = ($hits[0][0], $hits[0][1]);
                            }
                        }
                        $corrected_site_counts[$assigned]++ if ($l_out != 0 or $r_out != 0);
                    } else {
                        if ($l_perf_unshift > 0 or $r_perf_unshift > 0){ #ONE-SIDED PERFECT (COMPLETE)
                            $site_type_counts[6]++;
                            $assigned = 6;
                            ($l_out, $r_out) = ($l_perf_unshift, $r_perf_unshift);
                            if (abs($l_perf_unshift - $r_perf_unshift) >= $maxdiff/2){ #set a difference cutoff of half of the maximum difference
                                my $l_seq;
                                foreach (@{$read_seqs{$ref}{$site}{"left"}}){
                                    my ($count, $leng, $seq, $mm) = @{$_};
                                    next unless $mm eq "X";
                                    $l_seq = $seq unless ($l_seq and $leng == 16); #try to get 17bp sequence, but will take 16bp if no 17bp available
                                }
                                my $r_seq;
                                foreach (@{$read_seqs{$ref}{$site}{"right"}}){
                                    my ($count, $leng, $seq, $mm) = @{$_};
                                    next unless $mm eq "X";
                                    $r_seq = $seq unless ($r_seq and $leng == 16); #try to get 17bp sequence, but will take 16bp if no 17bp available
                                }
                                if ($l_out == 0){
                                    $l_seq = substr($fullseq{$ref}, ($site - 16), 17);
                                }
                                if ($r_out == 0){
                                    $r_seq = substr($fullseq{$ref}, ($site - 1), 17);
                                    $r_seq = reverse($r_seq);
                                    $r_seq =~ tr/ACTG/TGAC/;
                                }
                                if ($l_seq and $r_seq){
                                    my @order = ([$l_perf_unshift, $l_seq, "L"], [$r_perf_unshift, $r_seq, "R"]);
                                    @order = sort{$b->[0] <=> $a->[0]} @order;
                                    my ($top, $bot) = flank_assess($order[0][0], $order[0][1], $order[1][0], $order[1][1], $bar, $ref, 6, 0);
                                    ($l_out, $r_out) = ($top, $bot);
                                    ($l_out, $r_out) = ($bot, $top) if $order[0][2] eq "R";
                                }
                            }
                            $corrected_site_counts[$assigned]++ if ($l_out != $l_perf_unshift or $r_out != $r_perf_unshift);
                        } else { #ONE-SIDED IMPERFECT (COMPLETE)
                            $site_type_counts[7]++;
                            $assigned = 7;
                            ($l_out, $r_out) = (0, 0);
                            if (abs($l_imperf_unshift - $r_imperf_unshift) >= $maxdiff/2){ #set a difference cutoff of half of maximum difference
                                my %l_seqs;
                                foreach (@{$read_seqs{$ref}{$site}{"left"}}){
                                    my ($count, $leng, $seq, $mm) = @{$_};
                                    next unless $mm =~ m/\d+/;
                                    $seq = substr($seq, 1) if $leng == 17;
                                    $l_seqs{$seq} += $count;
                                }
                                my $l_seq;
                                if (%l_seqs){
                                    my @array;
                                    foreach my $seq (keys %l_seqs){
                                        push @array, ([$seq, $l_seqs{$seq}]);
                                    }
                                    @array = sort{$b->[1] <=> $a->[1]}@array;
                                    $l_out = $array[0][1];
                                    foreach (@{$read_seqs{$ref}{$site}{"left"}}){
                                    my ($count, $leng, $seq, $mm) = @{$_};
                                        next unless $mm =~ m/\d+/;
                                        if ($seq =~ m/$array[0][0]$/){
                                            $l_seq = $seq unless ($l_seq and $leng == 16); #try to get 17bp sequence, but will take 16bp if no 17bp available
                                        }
                                    }
                                }
                                
                                my %r_seqs;
                                foreach (@{$read_seqs{$ref}{$site}{"right"}}){
                                    my ($count, $leng, $seq, $mm) = @{$_};
                                    next unless $mm =~ m/\d+/;
                                    $seq = substr($seq, 1) if $leng == 17;
                                    $r_seqs{$seq} += $count;
                                }
                                my $r_seq;
                                if (%r_seqs){
                                    my @array;
                                    foreach my $seq (keys %r_seqs){
                                        push @array, ([$seq, $r_seqs{$seq}]);
                                    }
                                    @array = sort{$b->[1] <=> $a->[1]}@array;
                                    $r_out = $array[0][1];
                                    foreach (@{$read_seqs{$ref}{$site}{"right"}}){
                                    my ($count, $leng, $seq, $mm) = @{$_};
                                        next unless $mm =~ m/\d+/;
                                        if ($seq =~ m/$array[0][0]$/){
                                            $r_seq = $seq unless ($r_seq and $leng == 16); #try to get 17bp sequence, but will take 16bp if no 17bp available
                                        }
                                    }
                                }
                                
                                if ($l_out == 0){
                                    $l_seq = substr($fullseq{$ref}, ($site - 16), 17);
                                }
                                if ($r_out == 0){
                                    $r_seq = substr($fullseq{$ref}, ($site - 1), 17);
                                    $r_seq = reverse($r_seq);
                                    $r_seq =~ tr/ACTG/TGAC/;
                                }
                                
                                if ($l_seq and $r_seq){
                                    my @order = ([$l_out, $l_seq, "L"], [$r_out, $r_seq, "R"]);
                                    @order = sort{$b->[0] <=> $a->[0]} @order;
                                    my ($top, $bot) = flank_assess($order[0][0], $order[0][1], $order[1][0], $order[1][1], $bar, $ref, 7, 0);
                                    ($l_out, $r_out) = ($top, $bot);
                                    ($l_out, $r_out) = ($bot, $top) if $order[0][2] eq "R";
                                    unless ($l_out > 0 and $r_out > 0){
                                        ($l_out, $r_out) = (0,0);
                                    }
                                }
                            }
                            
                            unless ($l_perf_shift > 1 and $r_perf_shift > 1){
                                if ($l_perf_shift > 1 and $l_out > 1 and !$opt_h){
                                    ($r_out, $l_out) = shift_assess($r_out, $l_out, $l_perf_shift);
                                } elsif ($r_perf_shift > 1 and $r_out > 1 and !$opt_h){
                                    ($l_out, $r_out) = shift_assess($l_out, $r_out, $r_perf_shift);
                                }
                            }
                            $corrected_site_counts[$assigned]++ if ($l_out != 0 or $r_out != 0);
                        }
                    }
                }
                unless ($assigned){
                    $site_type_counts[8]++;
                    $assigned = 8;
                }
                
                #output read counts in format requested
                my $tot_out = $l_out + $r_out;
                my $outline = "-$site\t$l_out\n$site\t$r_out\n";
                $outline = "$l_out\t-$site\n$r_out\t$site\n" if $outform eq "essen";
                $outline = ">$ref\t$site\t$l_out\t$r_out\t$tot_out\n" if $outform eq "inseq";
                print $wig_out "$outline";
            }
            close ($wig_out);
            print STDERR "\rProcessing site #$site_count (of $num_sites) ... Done!\n";
            print "$id, $bar, $ref\n";
            print "# sites total: $site_count, # sites with at least $minrd reads: $site_gte_count\n";
            print "# of sites of each type (# of sites modified/filtered):\n";
            print "\t#1 both perfect, no imperfect (mismatch or shifted): $site_type_counts[1] ($corrected_site_counts[1])\n";
            print "\t#2 both perfect, minor imperfect on one or both sides: $site_type_counts[2] ($corrected_site_counts[2])\n";
            print "\t#3 both perfect, major imperfect on one or both sides: $site_type_counts[3] ($corrected_site_counts[3])\n";
            print "\t#4 perfect on one side, imperfect on the other: $site_type_counts[4] ($corrected_site_counts[4])\n";
            print "\t#5 imperfect only on both sides: $site_type_counts[5] ($corrected_site_counts[5])\n";
            print "\t#6 perfect on one side, 0 reads on other: $site_type_counts[6] ($corrected_site_counts[6])\n";
            print "\t#7 imperfect on one side, 0 reads on other: $site_type_counts[7] ($corrected_site_counts[7])\n";
            print "\t#8 unassigned: $site_type_counts[8] ($corrected_site_counts[8])\n\n";
            $site_type_counts[0] = $site_gte_count;
            $all_st{$ref}{$bar} = [@site_type_counts];
            $all_corr{$ref}{$bar} = [@corrected_site_counts];
        }
    }
}

#output summary
print "ID\tbarcode\treads_aligned\treads_align_TA\tperfect\tmismatch\tshifted\n";
foreach my $slice (@bc_list){
    my ($id, $bar) = @{$slice};
    if ($align_stats{$bar}){
        my @array = @{$align_stats{$bar}};
        print "$id\t$bar\t", join("\t", @array), "\n";
    }
}
print "\nreference\tID\tbarcode\ttotal_sites\ttype1\t1_corrected\ttype2\t2_corrected\ttype3\t3_corrected\ttype4\t4_corrected\ttype5\t5_corrected\ttype6\t6_corrected\ttype7\t7_corrected\ttype8\t8_corrected\n";
foreach my $ref (sort keys %refs){
    foreach my $slice (@bc_list){
        my ($id, $bar) = @{$slice};
        print "$ref\t$id\t$bar";
        my @st = @{$all_st{$ref}{$bar}};
        my @corr = @{$all_corr{$ref}{$bar}};
        for my $i (0 .. 8){
            print "\t$st[$i]";
            next if $i == 0;
            print "\t$corr[$i]";
        }
        print "\n";
    }
}


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