#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include "base_functions.h"
#include "gen_bounds_calculation.h"

using std::string;
using std::vector;
using std::istringstream;
using std::copy;


void GetGenBoundOfChainsList(string path, string input_name, 
                             const MU_TYPE& mu_type,
                             const Functional& functional, 
                             long double eps = 0.05) {
  string input_filename(path);
  input_filename.append(input_name);

  std::ifstream input_file;
  input_file.open(input_filename);
  
  string line;
  getline(input_file, line);
  istringstream linestream(line);
  
  int L, l, m;
  string output_name;
  linestream >> L >> l >> m >> output_name;
  
  string output_filename(path);
  output_filename.append(output_name);
  std::ofstream output_file;
  output_file.open(output_filename);
  output_file.close();
  
  int max_block_size = 10;
  vector<long double> gen_bounds;
  gen_bounds.reserve(max_block_size);

  int chains_count = 0;
  do {
    getline (input_file, line);
    if (!line.empty()) {
      vector<int> chain_errors(ReadChainErrors(line));
      long double gen_bound = GetChainGenBound(L, l, m, chain_errors, mu_type, functional, eps);
      gen_bounds.push_back(gen_bound);
      ++chains_count;
    }
    if ((chains_count == max_block_size) or (input_file.eof())) {
      output_file.open(output_filename, std::ios::app);
      copy(gen_bounds.begin(), gen_bounds.end(), std::ostream_iterator<long double>(output_file, "\n"));
      output_file.close();
      gen_bounds.clear();
      gen_bounds.reserve(max_block_size);
      chains_count = 0;
    }
    std::cout << gen_bounds.size() << "\n";
  } 
  while (!input_file.eof());
  input_file.close();
}


int main(int nargs, char**args) {
  srand(145);
  Functional functional = QEPS;
  MU_TYPE mu_type = ERM;
  long double eps = 0.05;
  GetGenBoundOfChainsList(args[1], args[2], mu_type, functional, eps);
  return 0;
}
