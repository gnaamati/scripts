#!/bin/env perl
use 5.12.0;
use warnings;
use Data::Dumper;

## Easy manipulation of sets of integers (arbitrary intervals)
use Set::IntRange;

## Command line...
use Getopt::Long;

## We have the 'query/target' concept here
my $side = 'query';

## Get commad line options
GetOptions( 
    ## Either query or target
    'side=s'            => \$side,
)
    or die "failure to communicate\n";

die "side must be either query or target\n"
    unless
        $side eq 'query'  ||
        $side eq 'target';

warn "using $side\n";

die "pass two range-files\n"
    unless @ARGV == 2;

my $rf1 = $ARGV[0];
my $rf2 = $ARGV[1];


## This 'first pass' is needed for efficieny below. Basically the 'set
## resize' operations are relatively costly, so by finding the maximum
## set size required per ID up front, we save a lot of work later...
warn "Parsing pass 1/2\n";

my %length;

for my $range ($rf1, $rf2){
    warn "$range\n";
    
    open( my $RANGE, '<', $range )
        or die "Failed to open '$range': $!\n";
    
    while(<$RANGE>){
        chomp;
        
        my ($id1, $st1, $en1, 
            $id2, $st2, $en2) = split "\t";
        
        my ($id, $st, $en) = ($id1, $st1, $en1);
        
        if ($side eq 'target'){
            if (defined $id2){
                ($id, $st, $en) = ($id2, $st2, $en2);
            }
        }
        
        $length{$id} = $en
            if $en > ($length{$id} || 0)
    }
    
    warn scalar keys %length, "\n";
}



warn "Parsing pass 2/2\n";

my %range_sets;
my %range_space;

for my $range ($rf1, $rf2){
    warn "$range\n";
    
    open( my $RANGE, '<', $range )
        or die "Failed to open '$range': $!\n";
    
    my $ranges = 0;
    my $range_space = 0;
    
    while(<$RANGE>){
        chomp;
        
        my ($id1, $st1, $en1, 
            $id2, $st2, $en2) = split "\t";
        
        my ($id, $st, $en) = ($id1, $st1, $en1);
        if ($st !~ /\d+/){
            next;
        }
        
        if ($side eq 'target'){
            if (defined $id2){
                ($id, $st, $en) = ($id2, $st2, $en2);
            }
            if ($st !~ /\d+/){
                next;
            }
        }
        
        #if ($id ne 'IWGSC_CSS_7DS_scaff_3923931'){
        #    next;
        #}

        $st++;
        #$en++;
        
        $ranges++;
        $range_space += $en-$st+1;

    
        $range_sets{$range}{$id} = Set::IntRange->new(1, $length{$id})
            unless exists $range_sets{$range}{$id};
        if ($st>$en){
            $st = $en-1;
            #$st--;
        }
        $range_sets{$range}{$id}->Interval_Fill($st, $en);
    }
    
    warn "loaded $ranges ranges on ",
    scalar keys %{$range_sets{$range}}, " seq-regions\n";
    
    my $range_space_no_overlap;
    $range_space_no_overlap += $_->Norm
        for values %{$range_sets{$range}};
    
    warn "range space is $range_space total ".
        "($range_space_no_overlap no-overlap) bp\n";
    
    ## Hide this value away...
    $range_space{$range} = $range_space_no_overlap;
}



## What is the range-range coverage?

my $coverage_sr = 0;
my $coverage_bp = 0;

for my $id (keys %{$range_sets{$rf1}}){
    next unless exists $range_sets{$rf2}{$id};
    
    $coverage_sr++;
    
    my $intersection = Set::IntRange->new(1, $length{$id});
    
    $intersection->Intersection($range_sets{$rf1}{$id},
                                $range_sets{$rf2}{$id});
    
    $coverage_bp += $intersection->Norm;
}

warn "coverage of xxxxx ranges on $coverage_sr seq-regions\n";
warn "coverage (non-overlapping overlap) is ", $coverage_bp, " bp\n";

warn $coverage_bp / $range_space{$rf1} * 100,
    " of non-overlapping range space\n";
warn $coverage_bp / $range_space{$rf2} * 100,
    " of non-overlapping range space\n";

print
    join("\t",
         (scalar keys %{$range_sets{$rf1}}),
         $range_space{$rf1},
         
         (scalar keys %{$range_sets{$rf2}}),
         $range_space{$rf2},

         $coverage_sr,
         $coverage_bp,
         
         $coverage_bp / $range_space{$rf1} * 100,
         $coverage_bp / $range_space{$rf2} * 100,
         
         $rf1,
         $rf2,
    ), "\n";
