#!/usr/bin/env perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Bio::SeqIO;
use feature qw /say/;

my @dbids = (1..100_000);
#my @dbids = (1..5);

my $registry = 'Bio::EnsEMBL::Registry';
my $reg_conf = '/homes/gnaamati/registries/prod3_gn.reg';
$registry->load_all($reg_conf);

##fetch a gene by its stable identifier
my $gene_adaptor = $registry->get_adaptor("triticum_aestivum", "core", "gene");

#my $gene = $gene_adaptor->fetch_by_stable_id('Traes_4DS_4BD85B5C7');
#my $genes = $gene_adaptor->fetch_all_by_biotype('protein_coding');
my $genes  = $gene_adaptor->fetch_all_by_dbID_list(\@dbids);
#my $count = scalar (@$genes);
#say $count;
my $count;
for my $gene (@$genes){
    if ($gene->stable_id !~ /Traes/){
        next;
    }
    $count++;
    my $transcripts = $gene->get_all_Transcripts;
    for my $t (@$transcripts){
        print ">",$gene->stable_id,"\n";
        
        my $translation = $t->translation;
        print $translation->seq, "\n";
    }
}
say $count;

#my $tgene = $gene->transform('');



