
# Example 2-2 Mutate DNA
# using a random number generator to randomly select bases to mutate
use strict;
use warnings;
# Declare the variables
# The DNA is chosen to make it easy to see mutations:
open OUTPUT,">output.txt";
my $DNA = 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC';
# $i is a common name for a counter variable, short for "integer"
my $i;
my $mutant;
# Seed the random number generator.
# time|$$ combines the current time with the current process id
srand(time|$$);
$mutant = mutate($DNA);
print "\nMutateDNA\n\n";
print "\nHereis the original DNA:\n\n";
print "$DNA\n";
print "\nHereis the mutant DNA:\n\n";
print "$mutant\n";
# Let's put it in a loop and watch that bad boy accumulate mutations:
print "\nHereare 10 more successive mutations:\n\n";
for ($i=0 ; $i < 10 ; ++$i) {
$mutant = mutate($DNA);
print "$mutant\n";
}
sub mutate {
my($dna) = @_;
my(@nucleotides) = ('A', 'C', 'G', 'T');
# Pick a random position in the DNA
my($position) = randomposition($dna);
#print $position."\n";
# Pick a random nucleotide
my($newbase) = randomnucleotide(@nucleotides);
# Insert the random nucleotide into the random position in the DNA
# The substrarguments mean the following:
# In the string $dnaat position $position change 1 character to
# the string in $newbase
substr($dna,$position,1,$newbase);
return $dna;
}
sub randomelement{
my(@array) = @_;
# Here the code is succinctly represented rather than
# “return $array[intrand scalar @array];”
return $array[rand @array];
}
# randomnucleotide
#
# A subroutine to select at random one of the four nucleotides
#
# WARNING: make sure you call srandto seed the
# random number generator before you call this function.
sub randomnucleotide{
my(@nucleotides) = ('A', 'C', 'G', 'T');
# scalar returns the size of an array.
# The elements of the array are numbered 0 to size-1
return randomelement(@nucleotides);
}
sub randomposition{
my($string) = @_;
# Notice the "nested" arguments:
#
# $string is the argument to length
# length($string) is the argument to rand
# rand(length($string))) is the argument to int
# int(rand(length($string))) is the argument to return
#
# rand returns a decimal number between 0 and its argument.
# intreturns the integer portion of a decimal number.
#
# The whole expression returns a random number
# between 0 and length-1,
# which is how the positions in a string are numbered in Perl.
#
#return intrand length $string;
return int(rand(length($string)));
}
print OUTPUT "";
close OUTPUT;
exit;