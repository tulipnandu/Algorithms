README

K-median clustering
Paul Hale and Tulip Nandu

To execute the file you can first compile the code with javac KClustering4.java and then
to execute in mode one:
java KClustering4 [input file] [k]
to execute in mode two:
java KClustering4 [input file] [k] [seed file]

Program uses the same input method as the first assignment as well as the hamming distance method.  
Once the program has the input file, it first tries to find the seed file and determine
the initial seeds.  If it cannot find it, the program generates k-random seeds.
The program then clusters the sequences based on hamming distance.  At the end of each
round of clustering it finds the sequence that is most representative of the cluster by 
finding the one with the least hamming distance to all other sequences.

The program again uses a modified bubble sort to arrange the sequences in order of cluster
ID.  The result is outputted to a file.

The folder also includes the java class file, our test sequence(mutatedsequence.txt), the program
used to generate that sequence written in python, a sample seed file to work with our test file(4 clusters).

The main clustering method was written togeather, the output and sort written by Paul, and the test sequence code 
written by Tulip.