##Parses a fasta file for various purposes
use strict;
use warnings;
use feature qw/say/;
use Data::Dumper;
use lib '/nfs/production/panda/ensemblgenomes/development/gnaamati/lib';
use FileReader qw(slurp read_file file2hash file2hash_tab);
use constant {
    TRUE  => 1,
    FALSE => 2,
};

my $fasta_dir = '/nfs/nobackup/ensemblgenomes/grabmuel/ensgen-sequences/ftp';

{
    my ($fasta_file) = @ARGV;
    if (@ARGV != 1){
        usage();
    }
    my $full_file = "$fasta_dir/$fasta_file";
    
    my $print_status = FALSE;
    open IN, "<", $full_file or die "can't open $full_file\n";
    while (my $line = <IN>){
        chomp($line);
        ##Check which component we are using and print if correct component
        #if ($line =~ />TGAC/){
        #if ($line =~ /IWGSC_CSS/){
        if ($line =~ /^>TGACv1_scaffold/){
            #$print_status = TRUE;
            if ($line =~ /^>TGACv1_scaffold_\d+_1/){
                $line =~ s/ spam:spam spam:TGACv1//g;
                $print_status = TRUE; 
            }
            else{
                $print_status = FALSE;
            }
        }
        if ($print_status == TRUE){
            say $line;
        }
    }
}

sub usage {
    say "Usage perl create_css_fasta [contig_file]";
    exit 0;
}
 

