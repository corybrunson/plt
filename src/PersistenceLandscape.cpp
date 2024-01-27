
#include "PersistenceLandscape.h"

// Section: List operations

PersistenceLandscape PLsum(List pl_list) {
  
  PersistenceLandscape sum_out = as<PersistenceLandscape>(pl_list[0]);
  
  for (int i = 1; i < pl_list.size(); i++) {
    sum_out = sum_out.add(as<PersistenceLandscape>(pl_list[i]));
  }
  
  return sum_out;
}

List PLdiff(List pl_list) {
  
  List diff_out;
  
  for (int i = 1; i < pl_list.size(); i++) {
    PersistenceLandscape
    diff_i = as<PersistenceLandscape>(pl_list[i]);
    diff_i = diff_i.add(as<PersistenceLandscape>(
      pl_list[i - 1]).scale(-1));
    diff_out.push_back(diff_i);
  }
  
  return diff_out;
}

PersistenceLandscape PLmean(List pl_list) {
  
  PersistenceLandscape mean_out = PLsum(pl_list);
  
  return mean_out.scale(1.0 / pl_list.size());
}

NumericMatrix PLdist(List pl_list, unsigned p) {
  
  // empty matrix
  unsigned n = pl_list.size();
  NumericMatrix dist_out(n, n);
  
  // not assuming symmetric distance calculation
  for (int i = 0; i != n; i++) {
    PersistenceLandscape
    pl_i = as<PersistenceLandscape>(pl_list[i]);
    for (int j = 0; j != n; j++) {
      if (j == i) {
        // same landscape
        dist_out(i, j) = 0.;
      } else {
        // different landscape
        PersistenceLandscape
        pl_j = as<PersistenceLandscape>(pl_list[j]);
        dist_out(i, j) = pl_i.distance(pl_j, p);
      }
    }
  }
  
  return dist_out;
}

double PLvar(List pl_list, unsigned p) {
  
  // average landscape
  PersistenceLandscape avg = PLmean(pl_list);
  
  // sum-squared distance
  double ssd = 0;
  
  for (size_t i = 0; i != pl_list.size(); ++i) {
    
    PersistenceLandscape
    pl_i = as<PersistenceLandscape>(pl_list[i]);
    
    double d = avg.distance(pl_i, p);
    
    ssd += d * d;
  }
  
  // sample standard deviation
  double var_out = ssd / pl_list.size();
  // double var_out = ssd / (pl_list.size() - 1.0);
  return var_out;
}

double PLsd(List pl_list, unsigned p) {
  
  double sd_out = PLvar(pl_list, p);
  
  sd_out = sqrt(sd_out);
  return sd_out;
}

// Section: Module

RCPP_EXPOSED_CLASS(PersistenceLandscape)

RCPP_MODULE(persistence_landscape_module) {
  using namespace Rcpp;
  
  class_<PersistenceLandscape>("PersistenceLandscape")
    
    // .constructor<NumericMatrix>(""
    // "Creates an exact PL from a PD as a 2-column numeric matrix")
    // .constructor<NumericMatrix, double, double, double>(""
    // "Creates a discrete PL from a PD as a 2-column numeric matrix")
    .constructor<NumericMatrix, bool, double, double, double>(""
    "Creates a PL from a PD as a 2-column numeric matrix")
    .constructor<NumericVector, NumericMatrix>(""
    "Creates a discrete PL from a time vector and level matrix")
    
    .method("isExact",
    &PersistenceLandscape::isExact,
    "Queries whether the underlying PL representation is exact")
    .method("xMin",
    &PersistenceLandscape::xMin,
    "Returns the infimum (left endpoint) of the PL support")
    .method("xMax",
    &PersistenceLandscape::xMax,
    "Returns the supremum (right endpoint) of the PL support")
    .method("xBy",
    &PersistenceLandscape::xBy,
    "Returns the resolution of the discrete representation")
    
    .method("getInternal",
    &PersistenceLandscape::getInternal,
    "Returns the internal tensor representation of the PL")
    .method("toDiscrete",
    &PersistenceLandscape::toDiscrete,
    "Casts an exact PL to a discrete one")
    .method("discretize",
    &PersistenceLandscape::discretize,
    "Casts an exact PL to a discrete one")
    
    .method("scale",
    &PersistenceLandscape::scale,
    "Multiplies this PL by a scalar")
    .method("add",
    &PersistenceLandscape::add,
    "Adds this PL to another")
    .method("abs",
    &PersistenceLandscape::abs,
    "Takes the absolute value of this PL")
    .method("inner",
    &PersistenceLandscape::inner,
    "Takes the inner product of this PL with another")
    
    .method("minimum",
    &PersistenceLandscape::minimum,
    "Finds the minimum value of one level of this PL")
    .method("maximum",
    &PersistenceLandscape::maximum,
    "Finds the maximum value of one level of this PL")
    
    .method("moment",
    &PersistenceLandscape::moment,
    "Computes the n^th moment of one level of this PL")
    .method("integrate",
    &PersistenceLandscape::integrate,
    "Computes the integral of this PL")
    .method("distance",
    &PersistenceLandscape::distance,
    "Takes the p-distance between this PL and another")
    .method("norm",
    &PersistenceLandscape::norm,
    "Computes the p-norm of this PL")
    .method("indicator",
    &PersistenceLandscape::indicator,
    "Multiplies this PL by a level-indexed set of indicator functions")
    .method("indicator_form",
    &PersistenceLandscape::indicator_form,
    "Computes the integral of the productof this PL with an indicator")
    ;
  
  Rcpp::function("PLsum", &PLsum);
  Rcpp::function("PLdiff", &PLdiff);
  Rcpp::function("PLmean", &PLmean);
  Rcpp::function("PLdist", &PLdist);
  Rcpp::function("PLvar", &PLvar);
  Rcpp::function("PLsd", &PLsd);
}
