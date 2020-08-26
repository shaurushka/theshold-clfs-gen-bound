#include <vector>
#include <iterator>
#include <cmath> // fmin, fabs
#include <string>
#include <sstream> 
#include <queue>

using std::vector;
using std::queue;
using std::min;
using std::max;
using std::string;

enum MU_TYPE {ERM, MAXD};

enum Functional {QEPS, CCV, EOFF};

long double GetChainGenBound(int L, int l, int m, 
                             const vector<int>& chain_errors,
                             const MU_TYPE& mu_type,
                             const Functional& functional,
                             long double eps = 0.05,
                             vector<long double>* contributions = NULL);

vector<int> ReadChainErrors(string line);

