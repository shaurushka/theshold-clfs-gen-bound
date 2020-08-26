# theshold-clfs-gen-bound
Here is presented the C++ algorithm for calculating the tight generalization bounds of threshold classifiers.

It calculates one of the generalization bounds CCV, EOF, and Qeps for one of the methods of learning ERM or MAXD. 
The type of the calculated bound and method of learning are specified in "main.cpp" in variables functional, mu_type and epsilon for Qeps. 

**Input**

Threshold classifiers are specified in input filename (would be better if it was of txt extension). 
First line contains parameters of the set, separated by space: 

   *\<full sample size\> \<train sample size\> \<noise\> \<output_name\>,*

where *noise* is the error on neutral sample where all classifiers return the same response, *output_name* is the name of output filename (would be better if it was of txt extension). 

Then line by line are specified the theshold classifier errors. Errors on one line are separated by space. Each line corresponds to the specific distribution of class samples.

**Compilation**
1. To compile project run makefile: write down in console
   
   *make*
2. It generates executable gen_bounds_calc. To run it write down in console 
   
   *./gen_bounds_calc \<path_to_file\> \<input_filename\>,* 

    where <input_filename> is text file with the threshold classifiers errors, and <path_to_file> is the path to input filename. Output filename is saved there.

**Output**

General bounds of threshold classifiers for each distribution are written line by line.

**Example**

Example of input and output filenames are input.txt and output.txt
  

