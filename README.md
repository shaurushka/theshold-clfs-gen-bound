# theshold-clfs-gen-bound
Here is presented the C++ algorithm for calculating the tight generalization bounds of threshold classifiers.

It calculates one of the generalization bounds CCV, EOF, and Qeps for one of the methods of learning ERM or MAXD. 
The type of the calculated bound and method of learning are specified in "main.cpp" in variables functional, mu_type and epsilon for Qeps. 

Threshold classifiers are specified in txt file. 
First line contains parameters of the set, separated by space: 
<full sample size> <train sample size> <noise> <output_name>
where noise is the error on neutral sample where all classifiers return the same response.
Then line by line are specified the theshold classifier errors. Errors on one line are separated by space. Each line corresponds to the specific distribution of class samples.

To compile project run makefile. It generates executable gen_bounds_calc. To run it write down in console 
./gen_bounds_calc <path_to_file> <input_filename>
where <input_filename> is text file with the threshold classifiers errors, <path_to_file> is the path to input filename. Output filename is saved there.
  

