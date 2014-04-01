#ifndef _COMMONFUNCTION_H_
#define _COMMONFUNCTION_H_

#include <algorithm>
#include <set>
#include "gsl/gsl_cdf.h"
#include "Utils.h"

//////////////////////////////////////////////////
// Statistics functions
//////////////////////////////////////////////////
template<class T>
class OrderFunction {
public:
OrderFunction(T& t): v(t) {};
  bool operator() (int i, int j)  {
    return v[i] < v[j];
  };
  const T& v;
};

/**
 * @param freq: 0.3, 0.2, 0.1, 0.4
 * will return 0-based rank, does not handle tie
 * @param ord:  2, 1, 0, 3
 *
 * algorithm:
 * make a pair like this: (0.3, 0), (0.2, 1), (0.1, 2), (0.4, 3)
 * sort by first element:
 * (0.1, 2), (0.2, 1), (0.3, 0), (0.4, 3)
 * extract second element:
 * 2, 1, 0, 3
 */
inline void order(std::vector<double>& freq, std::vector<int>* ord) {
  ord->resize(freq.size());
  for (size_t i = 0; i < freq.size(); ++i)
    (*ord)[i] = i;

  OrderFunction< std::vector<double> > func(freq);
  std::sort(ord->begin(), ord->end(), func);
};

inline void order(std::vector<int>& freq, std::vector<int>* ord) {
  ord->resize(freq.size());
  for (size_t i = 0; i < freq.size(); ++i)
    (*ord)[i] = i;

  OrderFunction< std::vector<int> > func(freq);
  std::sort(ord->begin(), ord->end(), func);
};

/**
 * Calculate rank, using average for ties.
 */
inline void calculateRank(std::vector<double>& freq, std::vector<int>* res) {
  std::vector<int> ord;
  order(freq, &ord);
  order(ord, res);
  
  // calculate the mean of rank if there are ties
};

inline double pnorm(double x) {
  return gsl_cdf_ugaussian_P(x);
};
inline double qnorm(double x) {
  return gsl_cdf_ugaussian_Pinv(x);
};
inline double qnorm(double x, double sigma) {
  return gsl_cdf_gaussian_Pinv(x, sigma);
};

inline double calculateMean(const std::vector<double>& v ){
  double s = 0.0;
  if (v.empty()) return 0.0;

  for (unsigned int i = 0; i < v.size(); ++i) {
    s += v[i];
  }
  return s / v.size();
}

// this calculate standard deviation by dividing N
// sd ^ 2 = \frac{1}{N}  (v_i - mean of v_i) ^ 2
inline double calculateSD(const std::vector<double>& v ){
  double s = 0.0;
  if (v.empty()) return 0.0;

  double mean = calculateMean(v);
  for (unsigned int i = 0; i < v.size(); ++i) {
    const double t = v[i] - mean;
    s += t * t;
  }
  return sqrt(s / v.size());
}

// this calculate standard deviation by dividing N
// sd ^ 2 = \frac{1}{N-1}  (v_i - mean of v_i) ^ 2
inline double calculateSampleSD(const std::vector<double>& v ){
  double s = 0.0;
  if (v.size() <= 1) return 0.0;

  double mean = calculateMean(v);
  for (unsigned int i = 0; i < v.size(); ++i) {
    const double t = v[i] - mean;
    s += t * t;
  }
  return sqrt(s / ( v.size() - 1));
}

inline void inverseNormalizeLikeMerlin(std::vector<double>* y){
  if (!y || !y->size()) return;
  const int n = y->size();
  std::vector<int> yRank;
  calculateRank(*y, &yRank);

  // change order to 1-based rank
  for (int i = 0; i < n; ++i ) {
    (*y)[i] = qnorm( ( 0.5 + yRank[i]) / n );
    // fprintf(stderr, "%zu - %g - %d \n", i,  (*y)[i], ord[i]);
  }
  /* double m = calculateMean(*y); */
  /* double sd = calculateSD(*y); */
  // fprintf(stderr, "mean = %g, sd = %g\n", m, sd);
}

// inverse normal a vector AND center it.
inline void inverseNormalizeLikeR(std::vector<double>* y){
  if (!y || !y->size()) return;
  const int n = y->size();
  std::vector<int> ord;
  order(*y, &ord);

  for (int i = 0; i < n; i++)
    (*y)[i] = ord[i];
  order(*y, &ord);

  double a;
  if ( n <= 10) {
    a = 0.375;
  } else {
    a = 0.5;
  }
  for (int i = 0; i < n ; i++)
    (*y)[i] = qnorm( ( 1.0 + ord[i] - a) / ( n + (1 - a) - a));
}

inline void standardize(std::vector<double>* y) {
  // center
  size_t n = y->size();
  double m = calculateMean(*y);
  double sd = calculateSD(*y);
  for (size_t i = 0; i < n ; i++) {
    (*y)[i] -= m;
    (*y)[i] /= sd;
  }
  // logger->info("Done: centering to 0.0 and scaling to 1.0 finished.");
};

//////////////////////////////////////////////////
// STL based functions
//////////////////////////////////////////////////

/**
 * Convert a @param string separated by @param sep to set (stored in @param s)
 */
inline void makeSet(const std::string& str, char sep, std::set<std::string>* s) {
  s->clear();
  if (str.empty())
    return;

  std::vector<std::string> fd;
  stringTokenize(str, ",", &fd);
  for (size_t i = 0; i < fd.size(); i++)
    s->insert(fd[i]);
}

inline void makeSet(const std::vector<std::string>& in, std::set<std::string>* s) {
  s->clear();
  if (in.empty())
    return;

  for (size_t i = 0; i < in.size(); i++)
    s->insert(in[i]);
}

/**
 * make a vector to map.
 *   when there is no duplciation: key is vector[i], value is i
 *   when there is duplication: key is vector[i], value is the index of the first appearance of vector[i]
 */
inline void makeMap(const std::vector<std::string>& in, std::map<std::string, int>* s) {
  s->clear();
  if (in.empty())
    return;

  for (size_t i = 0; i < in.size(); i++) {
    if (s->find(in[i]) != s->end()) continue;
    (*s)[in[i]] = i;
  }
}

/**
 * Test whether x contain unique elements
 */
inline bool isUnique(const std::vector<std::string>& x) {
  std::set<std::string> s;
  for (size_t i = 0; i < x.size(); i++) {
    s.insert(x[i]);
    if (s.size() != i + 1) {
      return false;
    }
  }
  return true;
}

#endif /* _COMMONFUNCTION_H_ */
