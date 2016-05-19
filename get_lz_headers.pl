##Get the css headers for the overlap test
use 5.12.0;
use autodie qw(open close);
use lib '/nfs/production/panda/ensemblgenomes/development/gnaamati/lib';
use FileReader qw(slurp read_file file2hash file2hash_tab);

{
    my ($a) = @ARGV;
    #my @lines = slurp($a);
    open IN, $a;


    while (my $line=<IN>) {
        chomp ($line);
        if ($line =~ /name/){
            next;
        }
        my @cols = split(/\t/, $line);
        say "$cols[6]\t$cols[9]\t$cols[10]";
    }
}

sub usage {
    say "Usage perl get_css_headers.pl [a] [b]";
    exit 0;
}
 
