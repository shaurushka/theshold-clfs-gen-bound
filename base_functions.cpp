#include "base_functions.h"

long double CalcCombination(int of_count, int from_count) {
  if ((from_count < of_count)||(of_count < 0) || (from_count < 0)) {
    return 0;
  }
  
  /* use symmetry of CalcCombination to optimize calculation */
  if (from_count - of_count < of_count) {
    of_count = from_count - of_count;
  }
  
  long double current_combination = 1;
  for (size_t of = 1; of <= of_count; ++of) {
    current_combination = current_combination * (from_count - of + 1) / of;
  }
  return current_combination;
}

long double H(int L, int l, int m, int s) {
  if ((s < 0) || (L < l) || (l < 0) || (L < 0)) {
    return 0;
  }
  long double sum = 0;
  
  for (int i = fmax(0, m - L + l); i <= fmin(fmin(l, m), s); ++i) {
    sum += CalcCombination(i, m) * CalcCombination(l - i, L - m);
  }
  return sum / CalcCombination(l, L);
}

long double TildeH(int L, int l, int D, int m, int z_0,  int s) {
  if ((s < 0)||(L < l)||(l < 0) || (L < 0)) {
    return 0;
  }
  return CalcCombination(l - z_0, L - D) * H(L - D, l - z_0, m, s) / CalcCombination(l, L);
}

long double HH(int L, int l, int D, int m, int z_0,  int s) {
  if ((s < 0)||(L < l)||(l < 0) || (L < 0)) {
    return 0;
  }
  return CalcCombination(l - z_0, L - D) * H(L - D, l - z_0, m, s);
}

int getIntOfValue(long double value) {
  int ceil_of_value = ceil(value);
  if (ceil_of_value - value < 1e-9) {
    return ceil_of_value;
  }
  else {
    return floor(value);
  }
}
