#include "gen_bounds_calculation.h"
#include "base_functions.h"

typedef vector<vector<int> > ChainMatrixType;

typedef vector<vector<long double> > LargeValuesMatrixType;

enum CHAIN_DIRECTION {LEFT = -1, RIGHT = 1};

enum EDGE_DIR {UP, DOWN};

class AlgContribution {
  long double contribution_qeps;
  long double contribution_ccv;
  long double contribution_eof;

public:
  AlgContribution() : contribution_qeps(0), contribution_ccv(0), contribution_eof(0) {}

  AlgContribution(long double contribution_qeps, 
                  long double contribution_ccv, 
                  long double contribution_eof)
    : contribution_qeps(contribution_qeps), 
      contribution_ccv(contribution_ccv),
      contribution_eof(contribution_eof) {}

  void Add(long double qeps, long double ccv, long double eof) {
    contribution_qeps += qeps;
    contribution_ccv += ccv;
    contribution_eof += eof;
  }

  long double GetContributionQeps() const {
    return contribution_qeps;
  }

  long double GetContributionCCV() const {
    return contribution_ccv;
  }

  long double GetContributionEOF() const {
    return contribution_eof;
  }
};

class GenBound {
  long double value;
  vector<long double> contributions;

public:
  GenBound() {}
    
  GenBound(size_t chain_size) : value(0) {
    contributions.assign(chain_size, 0);
  }

  void AddContribution(long double contribution, int alg_index) {
    contributions[alg_index] += contribution; 
    value += contribution;
  }

  long double GetValue() const {
    return value;
  }

  vector<long double> GetContributions() const {
    return contributions;
  }
};

class GenBoundsPack {
  GenBound qeps;
  GenBound ccv;
  GenBound eof;

public:
  GenBoundsPack() {}

  GenBoundsPack(size_t chain_size) 
    : qeps(GenBound(chain_size)), ccv(GenBound(chain_size)), eof(GenBound(chain_size)) {}

  void AddContribution(const AlgContribution& alg, double denominator, int alg_index) {
    qeps.AddContribution(alg.GetContributionQeps() / denominator, alg_index); 
    ccv.AddContribution(alg.GetContributionCCV() / denominator, alg_index); 
    eof.AddContribution(alg.GetContributionEOF() / denominator, alg_index); 
  }

  GenBound GetQeps() const {
    return qeps;
  }

  GenBound GetCCV() const {
    return ccv;
  }

  GenBound GetEOF() const {
    return eof;
  }
};

class ChainWithCorrectedGaps {
public:
  ChainWithCorrectedGaps() {}
  
  ChainWithCorrectedGaps(const ChainMatrixType& chain_with_gaps) {
    errors_vectors_with_gaps_ = chain_with_gaps;
    vector<int> chain_errors_count = GetChainErrors(chain_with_gaps);
    
    for (int index = 0; index + 1 < chain_errors_count.size(); ++index) {
      chain_with_no_gaps_errors_count_.push_back(chain_errors_count[index]);
      
      if (HammingDistance(chain_with_gaps[index + 1], chain_with_gaps[index]) > 1) {
        int left_classifier_correct_samples_count = 0;
        int left_classifier_incorrect_samples_count = 0;
        
        FindCorrectAndIncorrectSamplesCount(chain_with_gaps[index],
                                            chain_with_gaps[index + 1],
                                            &left_classifier_correct_samples_count,
                                            &left_classifier_incorrect_samples_count);
        
        if (left_classifier_incorrect_samples_count > 0) {
          --left_classifier_incorrect_samples_count;
        }
        else {
          --left_classifier_correct_samples_count;
        }
        
        int added_chain_begin_index = chain_with_no_gaps_errors_count_.size();
        
        AddMonotoneChain(left_classifier_correct_samples_count, INCREASE);
        AddMonotoneChain(left_classifier_incorrect_samples_count, DECREASE);
        
        int added_chain_end_index = chain_with_no_gaps_errors_count_.size();
        begin_end_of_added_chains_.push_back(ChainDescription(added_chain_begin_index,
                                                              added_chain_end_index));
      }
    }
    chain_with_no_gaps_errors_count_.push_back(chain_errors_count.back());
    GetIsInInitialChainIndicator();
    MapNewIndicesToOld();
  }
  
  vector<int> GetChainWithNoGapsErrorsCount() const {
    return chain_with_no_gaps_errors_count_;
  }
  
  int GetOldIndex(int index_in_chain_with_gaps) const {
    return new_index_to_old_map_[index_in_chain_with_gaps];
  }
  
  vector<bool> IsInInitialChainIndicator() const {
    return is_in_initial_chain_;
  }
  
private:
  enum MonotoneChainType {DECREASE, INCREASE};
  
  struct ChainDescription {
    int begin_index;
    int end_index;
    
    ChainDescription() : begin_index(0), end_index(0) {}

    ChainDescription(int begin_index_, int end_index_)
      : begin_index(begin_index_), end_index(end_index_) {}

  };

  ChainMatrixType errors_vectors_with_gaps_;
  vector<int> chain_with_no_gaps_errors_count_;
  vector<ChainDescription> begin_end_of_added_chains_;
  vector<ChainDescription> begin_end_of_initial_chains_;
  vector<bool> is_in_initial_chain_;
  vector<int> new_index_to_old_map_;

  int ErrorsCount(const ChainMatrixType& chain, int d) {
    int errors = 0;
    for (int index = 0; index < chain[d].size(); ++index) {
      errors += chain[d][index];
    }
    return errors;
  }

  vector<int> GetChainErrors(const ChainMatrixType& chain) {
    vector<int> chain_errors(chain.size());
    for (int index = 0; index < chain.size(); ++index) {
      chain_errors[index] = ErrorsCount(chain, index);
    }
    return chain_errors;
  }

  int HammingDistance(const vector<int>& first_vector,
                    const vector<int>& second_vector) {
    if (first_vector.size() != second_vector.size()) {
      throw std::bad_exception();
    }
    int distance = 0;
    for (int index = 0; index < first_vector.size(); ++index) {
      distance += std::abs(first_vector[index] - second_vector[index]);
    }
    return distance;
  }

  void AddMonotoneChain(size_t chain_size, const MonotoneChainType& type) {
    int errors_count = chain_with_no_gaps_errors_count_.back();
    for (int index = 0; index < chain_size; ++index) {
      if (type == INCREASE) {
        ++errors_count;
      }
      else {
        --errors_count;
      }
      chain_with_no_gaps_errors_count_.push_back(errors_count);
    }
  }
  
  void FindCorrectAndIncorrectSamplesCount(const vector<int>& left_classifier,
                                           const vector<int>& right_classifier,
                                           int* left_classifier_correct_samples_count,
                                           int* left_classifier_incorrect_samples_count) const {
    for (size_t sample_index = 0; sample_index < left_classifier.size(); ++sample_index) {
      if (left_classifier[sample_index] < right_classifier[sample_index]) {
        ++(*left_classifier_correct_samples_count);
      }
      else if (left_classifier[sample_index] > right_classifier[sample_index]) {
        ++(*left_classifier_incorrect_samples_count);
      }
    }
  }
  
  void AddInitialChainEdgeIndicators(int edge_begin_index, int edge_end_index) {
    for (int initial_chain_index = edge_begin_index;
         initial_chain_index < edge_end_index;
         ++initial_chain_index) {
      is_in_initial_chain_[initial_chain_index] = true;
    }
  }
  
  void GetIsInInitialChainIndicator() {
    begin_end_of_initial_chains_.clear();
    int initial_chain_begin = 0;
    int initial_chain_end;
    is_in_initial_chain_.resize(chain_with_no_gaps_errors_count_.size());
    
    
    for (int index = 0; index < begin_end_of_added_chains_.size(); ++index) {
      initial_chain_end = begin_end_of_added_chains_[index].begin_index;
      AddInitialChainEdgeIndicators(initial_chain_begin, initial_chain_end);
      begin_end_of_initial_chains_.push_back(ChainDescription(initial_chain_begin,
                                                              initial_chain_end));
      
      initial_chain_begin = begin_end_of_added_chains_[index].end_index;
    }
    initial_chain_end = chain_with_no_gaps_errors_count_.size();
    AddInitialChainEdgeIndicators(initial_chain_begin, initial_chain_end);
    begin_end_of_initial_chains_.push_back(ChainDescription(initial_chain_begin,
                                                            initial_chain_end));
  }
  
  void MapNewIndicesToOld() {
    size_t chain_size = is_in_initial_chain_.size();
    new_index_to_old_map_.assign(chain_size, -1);
    int index_in_chain_with_gaps = 0;
    for (int index = 0; index < chain_size; ++index) {
      if (is_in_initial_chain_[index]) {
        new_index_to_old_map_[index] = index_in_chain_with_gaps;
        ++index_in_chain_with_gaps;
      }
    }
  }
};

long double GetPathsCountFromLastDomain(int delta_from_middle,
                                        int train_error,
                                        const LargeValuesMatrixType& last_domain) {
  if (train_error < last_domain.size()) {
    int last_domain_size = last_domain[train_error].size();
    int last_domain_middle = (last_domain_size - 1) / 2;
    int last_domain_point = last_domain_middle + delta_from_middle;
    if (last_domain_point < last_domain_size && last_domain_point >= 0) {
      return last_domain[train_error][last_domain_point];
    }
  }
  return 0;
}

int GetMinDelta(int max_delta, bool is_in_initial_chain, const MU_TYPE& mu_type) {
  if ((!is_in_initial_chain) || (mu_type == MAXD)){
    return -max_delta;
  }
  return 0;
}

bool PessimististicConditionIsHeld(int first_alg_err, int current_alg_err, const CHAIN_DIRECTION& chain_dir) {
  return ((chain_dir == RIGHT && current_alg_err < first_alg_err) ||
          (chain_dir == LEFT && current_alg_err <= first_alg_err));
}

bool DeltaIsNotAcceptable(int delta, bool is_in_initial_chain, 
                          int first_alg_err, int current_alg_err, 
                          int l, int L, 
                          const CHAIN_DIRECTION& chain_dir, 
                          const MU_TYPE& mu_type) {
  if (mu_type == ERM) {
    return (delta == 0 &&
            is_in_initial_chain &&
            !PessimististicConditionIsHeld(first_alg_err, current_alg_err, chain_dir));
  }
  // (mu_type == MAXD) 
  else {
    long double lower_bound = (l * 1.0 / L) * (current_alg_err - first_alg_err);
    return (is_in_initial_chain && 
            ((delta < lower_bound) || 
             ((chain_dir == RIGHT) && (fabs(delta - lower_bound) < 1e-10))));
  }
}

LargeValuesMatrixType GetPathsIntoMiddleDomainCount(int first_alg_err,
                                                    int current_alg_err,
                                                    const LargeValuesMatrixType& last_domain,
                                                    bool is_in_initial_chain,
                                                    int l, 
                                                    int L, 
                                                    const EDGE_DIR& edge_dir,
                                                    const CHAIN_DIRECTION& chain_dir, 
                                                    const MU_TYPE& mu_type) {
  int max_train_err = last_domain.size() - 1;
  LargeValuesMatrixType current_domain(max_train_err + 1);
  
  for (int train_err = 0; train_err <= max_train_err; ++train_err) {
    current_domain[train_err].assign(2 * max_train_err + 1, 0);
    
    int min_delta = GetMinDelta(max_train_err, is_in_initial_chain, mu_type);
    for (int delta = min_delta; delta <= max_train_err; ++delta) {
      if (DeltaIsNotAcceptable(delta, is_in_initial_chain, first_alg_err, current_alg_err, 
                               l, L, chain_dir, mu_type)) {
        current_domain[train_err][max_train_err + delta] = 0;
      }
      else {
        if (edge_dir == UP) {
          current_domain[train_err][max_train_err + delta] +=
            GetPathsCountFromLastDomain(delta - 1, train_err, last_domain) +
            GetPathsCountFromLastDomain(delta, train_err, last_domain);
        }
        else {
          current_domain[train_err][max_train_err + delta] +=
            GetPathsCountFromLastDomain(delta, train_err, last_domain) +
            GetPathsCountFromLastDomain(delta + 1, train_err - 1, last_domain);
        }
      }
    }
  }
  return current_domain;
}

size_t GetSubChainEdgesCount(const CHAIN_DIRECTION& dir, int first_alg_index, int full_size) {
  if (dir == LEFT) {
    return first_alg_index;
  }
  else {
    return full_size - 1 - first_alg_index;
  }
}

LargeValuesMatrixType GetPathsIntoFinalDomainCount(const ChainWithCorrectedGaps& chain_with_corrected_gaps,
                                                   int first_alg_index,
                                                   int l, 
                                                   int L,
                                                   const CHAIN_DIRECTION& step_in_dir, 
                                                   const MU_TYPE& mu_type) {
  vector<int> chain_errors = chain_with_corrected_gaps.GetChainWithNoGapsErrorsCount();
  /* subchain_size - это количество ребер в цепи! */
  size_t subchain_size = GetSubChainEdgesCount(step_in_dir, first_alg_index, chain_errors.size());
  
  LargeValuesMatrixType last_domain(subchain_size + 1, vector<long double>(2 * subchain_size + 1));
  last_domain[0][subchain_size] = 1;
  
  for (int shift = 1; shift <= subchain_size; ++shift) {
    int current_alg_index = first_alg_index + step_in_dir * shift;
    bool is_in_initial_chain = chain_with_corrected_gaps.IsInInitialChainIndicator()[current_alg_index];
    
    EDGE_DIR edge_dir = DOWN;
    if (chain_errors[current_alg_index] > chain_errors[current_alg_index - step_in_dir]) {
      edge_dir = UP;
    }
    last_domain = GetPathsIntoMiddleDomainCount(chain_errors[first_alg_index],
                                                chain_errors[current_alg_index],
                                                last_domain,
                                                is_in_initial_chain,
                                                l, 
                                                L, 
                                                edge_dir,
                                                step_in_dir, 
                                                mu_type);
  }
  return last_domain;
}

LargeValuesMatrixType GetSubChainSplitsCount(const ChainWithCorrectedGaps& chain_with_corrected_gaps,
                                             int first_alg_index,
                                             int l, 
                                             int L, 
                                             const CHAIN_DIRECTION& chain_dir, 
                                             const MU_TYPE& mu_type) {
  LargeValuesMatrixType pathsToLastDomain(GetPathsIntoFinalDomainCount(chain_with_corrected_gaps,
                                                                       first_alg_index,
                                                                       l, 
                                                                       L, 
                                                                       chain_dir, 
                                                                       mu_type));
  int subchain_size = pathsToLastDomain.size() - 1;
  LargeValuesMatrixType splits_count(subchain_size + 1,
                                     vector<long double>(subchain_size + 1, 0));
  for (int train_edges = 0; train_edges <= subchain_size; ++train_edges) {
    for (int train_err = 0; train_err <= train_edges; ++train_err) {
      splits_count[train_edges][train_err] = 
        pathsToLastDomain[train_err][subchain_size + train_edges - 2 * train_err];
    }
  }
  return splits_count;
}

LargeValuesMatrixType GetNeutralSplitsCount(int train_edges_max, int alg_chain_err, int alg_neutral_err,
                                            int L, int l, int chain_edges_count, long double eps) {
  vector<vector<long double > > splits_count(train_edges_max + 1);
  for (int train_edges = 0; train_edges <= train_edges_max; ++train_edges) {
    int train_err_max = fmin(train_edges, alg_chain_err);
    splits_count[train_edges].resize(train_err_max + 1);
    for (int train_err = 0; train_err <= train_err_max; ++train_err) {
      int s = getIntOfValue(l * (alg_neutral_err + alg_chain_err - eps * (L - l)) / L) - train_err;
      if (s >= 0) {
        splits_count[train_edges][train_err] = HH(L, l, chain_edges_count, alg_neutral_err, train_edges, s);
      }
    }
  }
  return splits_count;
}

LargeValuesMatrixType GetNeutralValidError(int train_edges_max, int alg_chain_err, int m,
                                           int L, int l, int chain_edges_count) {
  LargeValuesMatrixType neutral_valid_err(train_edges_max + 1, vector<long double>(train_edges_max + 1));

  for (int train_edges = 0; train_edges <= train_edges_max; ++train_edges)
    for (int train_err = 0; train_err <= fmin(alg_chain_err, train_edges); ++train_err)
      for (int s = 0; s <= fmin(l - train_edges, m); ++s) 
        neutral_valid_err[train_edges][train_err] += CalcCombination(s, m) * 
                                                     CalcCombination(l - train_edges - s, L - chain_edges_count - m) * 
                                                     (alg_chain_err + m - s - train_err) * 1.0 / (L - l);
  return neutral_valid_err;
}

LargeValuesMatrixType GetDiscrepancy(int train_edges_max, 
                                     int alg_chain_err, int m,
                                     int L, int l, int chain_edges_count) {
  LargeValuesMatrixType discrepancy(train_edges_max + 1, vector<long double>(train_edges_max + 1));

  for (int train_edges = 0; train_edges <= train_edges_max; ++train_edges)
    for (int train_err = 0; train_err <= fmin(alg_chain_err, train_edges); ++train_err)
      for (int s = 0; s <= fmin(l - train_edges, m); ++s) 
        discrepancy[train_edges][train_err] += 
          CalcCombination(s, m) * CalcCombination(l - train_edges - s, L - chain_edges_count - m) * 
          fmax(alg_chain_err + m - (s + train_err) * 1.0 * L / l, 0) * 1.0 / (L - l);
  return discrepancy;
}

LargeValuesMatrixType GetFunctionalValues(int train_edges_max, int alg_chain_err, int alg_neutral_err,
                                          int L, int l, int chain_edges_count, long double eps, 
                                          const Functional& functional) {
  LargeValuesMatrixType functional_values;
  switch (functional) {
    case QEPS: 
      functional_values = GetNeutralSplitsCount(train_edges_max, alg_chain_err, alg_neutral_err,
                                                L, l, chain_edges_count, eps);
      break;
    case CCV: 
      functional_values = GetNeutralValidError(train_edges_max, alg_chain_err, alg_neutral_err,
                                               L, l, chain_edges_count);
      break;
    case EOFF: 
      functional_values = GetDiscrepancy(train_edges_max, alg_chain_err, alg_neutral_err,
                                         L, l, chain_edges_count);
      break;
  }
  return functional_values;
}

int LeftSampleErrorCount(const ChainMatrixType& chain,
                         int d,
                         const vector<int>& graph_edges) {
  int err_count = 0;
  for (int index = 0; index < d; ++index) {
    err_count += chain[d][graph_edges[index]];
  }
  return err_count;
}

void GetAlgContribution(int alg_index,
                        const ChainWithCorrectedGaps& corrected_chain,
                        const vector<int> &corrected_chain_errors,
                        const ChainMatrixType &corrected_chain_matrix,
                        vector<int> &graph_edges,
                        int chain_edges_count,
                        int L,
                        int l,
                        int m,
                        long double eps, 
                        const MU_TYPE& mu_type, 
                        AlgContribution* contribution) {
  LargeValuesMatrixType left_splits = GetSubChainSplitsCount(corrected_chain, alg_index, l, L, LEFT, mu_type);
  LargeValuesMatrixType right_splits = GetSubChainSplitsCount(corrected_chain, alg_index, l, L, RIGHT, mu_type);
  int alg_chain_err = corrected_chain_errors[alg_index] - m;
  int train_edges_max = fmin(l, chain_edges_count);
  LargeValuesMatrixType neutral_splits(GetFunctionalValues(train_edges_max, alg_chain_err, m,
                                                           L, l, chain_edges_count, eps, QEPS));
  LargeValuesMatrixType valid_err(GetFunctionalValues(train_edges_max, alg_chain_err, m,
                                                      L, l, chain_edges_count, eps, CCV));
  LargeValuesMatrixType discrepancy(GetFunctionalValues(train_edges_max, alg_chain_err, m,
                                                        L, l, chain_edges_count, eps, EOFF));
  int left_chain_err = LeftSampleErrorCount(corrected_chain_matrix, alg_index, graph_edges);
  
  for (int left_train_edges = 0; left_train_edges <= fmin(l, alg_index); ++left_train_edges) {
    for (int right_train_edges = 0; 
         right_train_edges <= fmin(l, chain_edges_count - alg_index); 
         ++right_train_edges) {
      for (int left_train_err = 0; 
           left_train_err <= fmin(left_train_edges, left_chain_err); 
           ++left_train_err) {
        for (int right_train_err = 0; right_train_err <= fmin(right_train_edges, alg_chain_err - left_chain_err);
             ++right_train_err) {
          int train_edges = left_train_edges + right_train_edges;
          int train_err = left_train_err + right_train_err;
          if ((train_edges <= train_edges_max) and (train_err <= fmin(train_edges, alg_chain_err))) {
            long double chain_splits = left_splits[left_train_edges][left_train_err] *
                                       right_splits[right_train_edges][right_train_err];
            contribution->Add(chain_splits * neutral_splits[train_edges][train_err], 
                              chain_splits * valid_err[train_edges][train_err], 
                              chain_splits * discrepancy[train_edges][train_err]); 
          }
        }
      }
    }
  }
}

ChainMatrixType GenerateChainMatrixFromErrorVector(int L,
                                                   const vector<int>& errors) {
  size_t chain_size = errors.size();
  ChainMatrixType chain(chain_size);
  
  chain[0].resize(L);
  queue<int> zero_indices;
  queue<int> one_indices;
  for (int index = 0; index < errors[0]; ++index) {
    one_indices.push(index);
    chain[0][index] = 1;
  }
  for (int index = errors[0]; index < L; ++index) {
    zero_indices.push(index);
    chain[0][index] = 0;
  }
  int prev_algor_error = errors[0];
  for (int index = 1; index < chain_size; ++index) {
    chain[index] = chain[index - 1];
    if (errors[index] < prev_algor_error) {
      int diff = prev_algor_error - errors[index];
      while (diff > 0) {
        int one_index = one_indices.front();
        chain[index][one_index] = 0;
        one_indices.pop();
        zero_indices.push(one_index);
        --diff;
      }
    }
    else {
      int diff = errors[index] - prev_algor_error;
      while (diff > 0) {
        int zero_index = zero_indices.front();
        chain[index][zero_index] = 1;
        zero_indices.pop();
        one_indices.push(zero_index);
        --diff;
      }
    }
    prev_algor_error = errors[index];
  }
  return chain;
}

vector<int> GetGraphEdges(const ChainMatrixType& chain) {
  size_t chain_size = chain.size();
  if (chain_size < 1) {
    return vector<int>(0);
  }
  
  vector<int> edges(chain_size - 1);
  for (int index = 1; index < chain_size; ++index) {
    int element_index = 0;
    while (chain[index][element_index] == chain[index - 1][element_index]) {
      ++element_index;
    }
    edges[index - 1] = element_index;
  }
  return edges;
}

void GetChainGenBounds(int L, int l, int m, long double eps,
                       const vector<int>& chain_errors,
                       GenBoundsPack* gen_bounds, 
                       const MU_TYPE& mu_type) {
  ChainMatrixType chain_with_gaps(GenerateChainMatrixFromErrorVector(L, chain_errors));
  if (chain_with_gaps.size() <= 1) {
    throw std::bad_exception();
  }
  ChainWithCorrectedGaps chain_with_corrected_gaps(chain_with_gaps);
  vector<int> corrected_chain_errors(chain_with_corrected_gaps.GetChainWithNoGapsErrorsCount());
  ChainMatrixType corrected_chain_matrix(GenerateChainMatrixFromErrorVector(
                                              L, 
                                              corrected_chain_errors));
  vector<int> graph_edges = GetGraphEdges(corrected_chain_matrix);
  
  long double denominator = CalcCombination(l, L);
  (*gen_bounds) = GenBoundsPack(chain_with_gaps.size());
  
  size_t chain_size = corrected_chain_errors.size();
  for (int alg_index = 0; alg_index < chain_size; ++alg_index) {
    if (chain_with_corrected_gaps.IsInInitialChainIndicator()[alg_index]) {
      AlgContribution contribution;
      GetAlgContribution(alg_index,
                         chain_with_corrected_gaps,
                         corrected_chain_errors,
                         corrected_chain_matrix,
                         graph_edges,
                         chain_size - 1,
                         L, l, m, eps, 
                         mu_type,
                         &contribution);
      gen_bounds->AddContribution(contribution, 
                                  denominator, 
                                  chain_with_corrected_gaps.GetOldIndex(alg_index));
    }
  }
}

long double GetChainGenBound(int L, int l, int m, 
                             const vector<int>& chain_errors,
                             const MU_TYPE& mu_type,
                             const Functional& functional,
                             long double eps,
                             vector<long double>* contributions) {
  GenBoundsPack gen_bounds;
  GetChainGenBounds(L, l, m, eps, chain_errors, &gen_bounds, mu_type); 
  GenBound gen_bound;
  switch (functional) {
    case QEPS: 
      gen_bound = gen_bounds.GetQeps();
      break;
    case CCV: 
      gen_bound = gen_bounds.GetCCV();
      break;
    case EOFF: 
      gen_bound = gen_bounds.GetEOF();
      break;
  }
  if (contributions) {
    (*contributions) = gen_bound.GetContributions();
  }
  return gen_bound.GetValue();
}

vector<int> ReadChainErrors(string line) {
  vector<int> chains_errors;
  
  std::istringstream linestream(line);
  copy(std::istream_iterator<int>(linestream), 
       std::istream_iterator<int>(), 
       std::back_inserter(chains_errors));
  return chains_errors;
}

