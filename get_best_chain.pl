#========================================= 
#Filters the hsps by score and then selects the highest scoring HSPs 
# for the filtered results (as long as the match is possible
#========================================= 

use strict;
use warnings;
use feature qw/say/;
use Data::Dumper;
use lib '/nfs/production/panda/ensemblgenomes/development/gnaamati/lib';
use FileReader qw(slurp read_file file2hash file2hash_tab line2hash);
use Set::IntRange;
use constant {
    PASS => 1,
    FAIL => 0,
};

{
    my ($file) = @ARGV;
    if (@ARGV < 1){
        usage();
    }
    my @lines = slurp($file);
    my $header = shift(@lines);
    say $header;

    my $count = 0;
    my $no_match_count = 0;

    ##Get the hsps
    my @hsps = get_hsps($header, @lines);
    warn 'hsps parsed - sorting';

    my $chain_hash;
    my $count;
    ##Create all chains with scores
    for my $h (@hsps){
        $count++;
        
        my ($name1,$name2)   = ($h->{'name1'},$h->{'name2'});
        my ($score, $strand) = ($h->{'#score'}, $h->{'strand2'});
        
        ##Update the chain and score according to names and strand
        push (@{$chain_hash->{$name2}{"$name1:$strand"}{'hsps'}}, $h);
        $chain_hash->{$name2}{"$name1:$strand"}{'score'} += $score;
       
        if ($count % 500000 == 0){
            warn $count;
        }
    }

    ##Get best chains
    for my $key (keys %$chain_hash){

        my @chains;
        my $seq_hash = $chain_hash->{$key};
        for my $k (keys %$seq_hash){
            push (@chains, $seq_hash->{$k});
        }

        ##Sort each chain according to the score and choose the best chain
        my @sorted_chains = sort {$b->{score} <=> $a->{score}} @chains;
        my $best_chain = $sorted_chains[0];
        
        ##print out the best chain
        for my $hsp (@{ $best_chain->{'hsps'} }){
            say $hsp->{to_string};
        }
    }

}
        
#======================================== 
sub get_match {
#======================================== 
    my ($hash,$name,$size, $set) = @_;
    my $existing_set = $hash->{$name};

    my $match;
    #If it doesn't exist, add new source set to hash
    if (!$existing_set){
        $hash->{$name} = $set; 
        $match = PASS;
    }
    else{
        #Set exists, look for intersection, if found FAIL)
        my $set1 = Set::IntRange->new(0, $size);
        $set1->Intersection($set,$existing_set);
        
        #No intersection - PASS
        if ($set1->is_empty){
            $match = PASS;
            $set1->Union($set, $existing_set);
            $hash->{$name} = $set1;
        }
        else{
            #Intersection found - FAIL
            $match = FAIL;
        }
    }
    return ($hash,$match);
}

#======================================== 
sub get_hsps {
#======================================== 
    my ($header, @lines) = @_;
    my @hsps;
    my $c = 0;
    for my $line (@lines){
        my $hsp = line2hash($header, $line);
        $hsp->{to_string} = $line;
        push (@hsps, $hsp);
        if ($c % 500000 == 0){
            warn $c;
        }
        $c++;
    }
    return @hsps;
}

sub usage {
    say "Usage perl get_best_chain.pl [chained_file]";
    exit 0;
}

=comment
        #Create new target set for range
        my $target_set  = Set::IntRange->new(0,$h->{size2});
        $target_set->Interval_Fill($h->{zstart2}, $h->{end2});
        
        #Create new source set for range
        my $new_source_set = Set::IntRange->new(0, $h->{'size1'});
        $new_source_set->Interval_Fill($h->{'zstart1'}, $h->{'end1'});
        
        #Get existing source set
        my $existing_source_set = $source_hash->{$name};
        my $existing_source_set = $target_hash->{$name};

        #If it doesn't exist, add new source set to hash (PASS)
        if (!$existing_source_set){
            $source_hash->{$name} = $new_source_set; 
            $source_hash->{$name} = $new_source_set; 
            $match = PASS;
        }
        else{
            #Set exists, look for intersection, if found FAIL)
            my $set1 = Set::IntRange->new(0, $h->{'size1'});
            $set1->Intersection($new_source_set, $existing_source_set);
            
            #No intersection - PASS
            if ($set1->is_empty){
                $match = PASS;
                $set1->Union($new_source_set, $existing_source_set);
                $source_hash->{$name} = $set1;
            }
            else{
                #Intersection found - FAIL
                $match = FAIL;
                $no_match_count++;
                say $name;
                say $existing_source_set->to_Enum();
                die('');
            }
        }
        
        #Print if match is OK
        if ($match){
        $target_set->Interval_Fill($h->{zstart2}, $h->{end2});
            say "$h->{'name2'}\t$h->{zstart2}\t$h->{end2}";
        }
        if ($count % 1000 == 0){
            warn $count;
        }
    }

    say "no matches = $no_match_count\n";

}

=cut

 
