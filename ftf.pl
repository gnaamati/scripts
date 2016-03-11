use strict;
use warnings;
use feature qw/say/;

my ($path) = @ARGV;
$path =~ s/::/\//g;
$path .= '.pm';

require $path;
say $INC{$path};
