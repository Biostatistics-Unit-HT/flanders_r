#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List compute_cs99_cpp(IntegerVector p, IntegerVector ii, NumericVector xx, int nc, double threshold) {
  List out(nc);
  for (int j = 0; j < nc; ++j) {
    int start = p[j];
    int endp = p[j+1];
    if (start >= endp) { out[j] = IntegerVector(); continue; }
    int len = endp - start;
    std::vector<std::pair<double,int>> vals;
    vals.reserve(len);
    for (int k = 0; k < len; ++k) {
      double v = std::exp(xx[start + k]);
      int row = ii[start + k] + 1; // convert to 1-based for R
      vals.emplace_back(v, row);
    }
    std::sort(vals.begin(), vals.end(), [](const std::pair<double,int>& a, const std::pair<double,int>& b){
      return a.first > b.first;
    });
    double cs = 0.0;
    std::vector<int> keep;
    keep.reserve(vals.size());
    for (size_t t = 0; t < vals.size(); ++t) {
      cs += vals[t].first;
      keep.push_back(vals[t].second);
      if (cs >= threshold) break;
    }
    out[j] = IntegerVector(keep.begin(), keep.end());
  }
  return out;
}
