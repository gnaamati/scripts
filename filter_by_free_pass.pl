#Filters the best HSPs according to score
#HSPs completely overlapping a gene will get a bonus score and will be ranked higher
use strict;
use warnings;
use feature qw/say/;
use Data::Dumper;
use lib '/nfs/production/panda/ensemblgenomes/development/gnaamati/lib';
use FileReader qw(slurp read_file file2hash file2hash_tab line2hash get_slice_adaptor);
use Set::IntRange;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception;

use constant {
    PASS  => 1,
    ON    => 1,
    FAIL  => 0,
    RANGE => 100_000,
};

my $hsp_hash;

##using operator overload see
## http://search.cpan.org/~stbey/Set-IntRange-5.2/IntRange.pm
{
    my ($gene_file,$lastz_file) = @ARGV;
    if (@ARGV != 2){
        usage();
    }
    
    ##Create a gene hash (where the genes are mapped)
    my $gene_hash = create_gene_hash($gene_file);

    ##Get hsps (high scoring alignment pairs between CSS and TGAC) 
    ## giving a free pass to ones that completely overlap genes
    my @hsps = get_hsps($lastz_file, $gene_hash);

    my $tgac_hash;
    my $tgac_match;
    my $css_hash;
    my $css_match;
    my $count;

    ##Get slice adaptor to convert name to ID
    my $registry = '/homes/gnaamati/registries/prod1.reg';
    my $slice_adaptor = get_slice_adaptor($registry,'triticum_aestivum');

    ##Get all the HSP matches
    for my $h (sort {$b->{'#score'} <=> $a->{'#score'}} @hsps){
        ($tgac_match, $css_match) = (FAIL,FAIL);

        ##Get vals from hsp
        my ($name1,$name2) = ($h->{'name1'},$h->{'name2'});
        my ($size1,$size2) = ($h->{'size1'},$h->{'size2'});
        
        ##Print and continue if free pass is on
        if ($h->{'free_pass'}){
            print_hsp($h, $slice_adaptor);
            next;
        }

        ##Create new TGAC set for range
        my $tgac_set = Set::IntRange->new(0,$size1);
        $tgac_set->Interval_Fill($h->{'zstart1'}, $h->{'end1'});

        ##Create new css set for range
        my $css_set  = Set::IntRange->new(0,$size2);
        $css_set->Interval_Fill($h->{'zstart2'}, $h->{'end2'});

        ##Get the tgac_match
        ($tgac_hash, $tgac_match) = get_match($tgac_hash,$name1,$size1,$tgac_set);

        ##Get the css_match
        ($css_hash, $css_match) = get_match($css_hash,$name2,$size2,$css_set);

        ##Print the results if there is a match (with seq ids)
        if ($css_match and $tgac_match){
            print_hsp($h, $slice_adaptor);
        }
        
        ##increment count and print if needed
        $count++;
        if ($count % 500000 == 0){
            warn $count;
        }
        
    }
}

#======================================== 
sub print_hsp{
#======================================== 
    my ($hsp, $slice_adaptor) = @_;
    my ($name1,$name2) = ($hsp->{'name1'},$hsp->{'name2'});
    my $strand = $hsp->{strand};
    
    my ($start1, $start2) = ($hsp->{zstart1}, $hsp->{zstart2});
    my ($end1, $end2)     = ($hsp->{end1}, $hsp->{end2});


    ##Fix the strand, and update start bacause start is origin-zero (see LASTZ doc online)
    my $strand = $hsp->{'strand2'};
    if ($strand eq '+'){
        $strand = 1;
    }
    $start1++;
    $start2++;

    ##Print the results if there is a match (with seq ids)
    my ($asm_seq, $cmp_seq) = get_seq_ids($slice_adaptor, $name1,$name2);

    ##Make sure this isn't a duplicate
    my $hsp_string = "$asm_seq-$cmp_seq-$start1-$end1-$start2-$end2-$strand";
    if ($hsp_hash->{$hsp_string}){
        return;
    }
    $hsp_hash->{$hsp_string}++;

    ##SQL for insert
    my $sql = qq{
    insert into assembly (asm_seq_region_id,cmp_seq_region_id,asm_start,asm_end,cmp_start,cmp_end,ori)
    values ($asm_seq, $cmp_seq, $start1, $end1, $start2, $end2,$strand);
    };

    say $sql;
}

#======================================== 
sub get_seq_ids {
#======================================== 
    my ($slice_adaptor, $name1, $name2) = @_;
    my $slice1   = $slice_adaptor->fetch_by_region('toplevel',$name1);
    my $asm_seq  = $slice_adaptor->get_seq_region_id($slice1);
    
    my $slice2   = $slice_adaptor->fetch_by_region('toplevel',$name2);
    my $cmp_seq  = $slice_adaptor->get_seq_region_id($slice2);

    return ($asm_seq, $cmp_seq);
}


#======================================== 
sub get_match {
#======================================== 
    my ($hash,$name,$size, $set) = @_;
    my $existing_set = $hash->{$name};

    my $match;
    #If it doesn't exist, add new TGAC set to hash
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
sub create_gene_hash {
#========================================  
    my ($gene_file) = @_;   
    
    my @lines = slurp($gene_file);
    my $gene_hash = {};
    my ($match, $count) = 0;
    for my $line (@lines){
        my ($name, $start, $end) = split(/\t/, $line);
        
        ##Create gene set
        my $gene_set  = Set::IntRange->new(0,RANGE);
        $gene_set->Interval_Fill($start,$end);

        push (@{$gene_hash->{$name}}, $gene_set);
    }

    return $gene_hash;
}

#======================================== 
sub get_hsps {
#======================================== 
    my ($file,$gene_hash) = @_;
    my @lines = slurp($file);
    my $header = shift(@lines);
    my @hsps;
    my $c = 0;
    my $gene_count = 0;
    for my $line (@lines){
        my $hsp = line2hash($header, $line);
        $hsp->{to_string} = $line;

        ##Special case for negative strand: 
        ##Start = size - end; End = size - start 
        if ($hsp->{'strand2'} eq '-'){
            my $start = $hsp->{'size2'} - $hsp->{'end2'};
            my $end   = $hsp->{'size2'} - $hsp->{'zstart2'};
            $hsp->{'zstart2'} = $start;
            $hsp->{'end2'}  = $end;
        }

        #Create new css set for range
        my $name   = $hsp->{'name2'};
        my $css_set  = Set::IntRange->new(0,RANGE);
        $css_set->Interval_Fill($hsp->{'zstart2'}, $hsp->{'end2'});
        

        ##Get the gene sets for the given name (if exists)
        my @gene_sets = ();
        my $set_count = 0;
        if ($gene_hash->{$name}){
            @gene_sets = @{$gene_hash->{$name}};
        }

        ##Try to find a match for the gene
        for my $gene_set (@gene_sets){
            ##Ignore empty (deleted before)
            if (!$gene_set){
                next;
            }
            ##If we find a match, give hsp free pass and remove gene set
            if ($css_set >= $gene_set){
                $hsp->{'free_pass'} = ON;
                $gene_count++;
                
                ##Remove gene set and update gene hash
                delete($gene_sets[$set_count]);
                $gene_hash->{$name} = \@gene_sets;
                last;
            }
            $set_count++;
        }

        push (@hsps, $hsp);
        if ($c % 500000 == 0){
            warn $c;
        }
        $c++;
    }

    return @hsps;

}

sub usage {
    say "Usage perl get_gene_match [GENE_FILE] [LASTZ_FILE]";
    exit 0;
}


