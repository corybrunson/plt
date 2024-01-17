
#include <Rcpp.h>
using namespace std;
using namespace Rcpp;

// Section: Helpers

// TODO: Make this a user-settable option through {Rcpp}.
double epsi = 0.000005;
inline bool almostEqual(double a, double b) {
  if (fabs(a - b) < epsi)
    return true;
  return false;
}

// compare birth-death pairs
bool compareFeatures(
    std::pair<double, double> x,
    std::pair<double, double> y
) {
  if (x.first < y.first) {
    return true;
  }
  // x.first >= y.first
  if (x.first > y.first) {
    return false;
  }
  // x.first == y.first
  if (x.second > y.second) {
    return true;
  } else {
    return false;
  }
}

// named functions for infix operators
inline double addOp(double x, double y) { return x + y; }
inline double subOp(double x, double y) { return x - y; }

// get linear function value at a point
double interpolatePoint(
    std::pair<double, double> land1,
    std::pair<double, double> land2,
    double x
) {
  // we assume that land1.first <= x <= land2.first and that land1 and land2
  // are points between which we will put the line segment
  double a = (land2.second - land1.second) / (land2.first - land1.first);
  double b = land1.second - a * land1.first;
  
  return (a * x + b);
}

// construct a vector of pairs from a NumericMatrix
std::vector<std::pair<double, double>> persistencePairs(NumericMatrix x) {
  
  // restructure the NumericMatrix to a vector of pairs, imputing bounds
  std::vector<std::pair<double, double>> pd;
  for (int i = 0; i < x.nrow(); i++)
    pd.push_back(std::make_pair(x(i,0), x(i,1)));
  
  // ensure that each birth (death) is in the first (second) column
  // FIXME: Harmonize this step with extended persistence data? -JCB
  for (size_t i = 0; i != pd.size(); ++i) {
    if (pd[i].second < pd[i].first) {
      double sec = pd[i].second;
      pd[i].second = pd[i].first;
      pd[i].first = sec;
    }
  }
  
  return pd;
}

// retrieve births and deaths from characteristic points (in landmark space)
double birth(std::pair<double, double> a) { return a.first - a.second; }
double death(std::pair<double, double> a) { return a.first + a.second; }

int TDIndex2(int X, int Y, int x, int y) {
  return x + X * y;
}
std::vector<NumericVector> getInternalExact(
    std::vector<std::vector<std::pair<double, double>>> input) {
  
  std::vector<NumericVector> out_d;
  
  for (int j = 0; j < input.size(); j++) {
    Dimension d(input[j].size(), 2);
    NumericVector out(input[j].size() * 2);
    
    for (int i = 0; i < input[j].size(); i++) {
      out[TDIndex2(input[j].size(), 2, i, 0)] =
        input[j][i].first;
      out[TDIndex2(input[j].size(), 2, i, 1)] =
        input[j][i].second;
      
      if (input[j][i].first == INT_MAX)
        out[TDIndex2(input[j].size(), 2, i, 0)] = R_PosInf;
      
      if (input[j][i].first == INT_MIN)
        out[TDIndex2(input[j].size(), 2, i, 0)] = R_NegInf;
    }
    
    out.attr("dim") = d;
    out_d.push_back(out);
  }
  
  return out_d;
}

int TDIndex(int X, int Y, int Z, int x, int y, int z) {
  return x + X * (y + Y * z);
}
NumericVector getInternalDiscrete(
    std::vector<std::vector<std::pair<double, double>>> input) {
  
  Dimension d(input.size(), input[0].size(), 2);
  NumericVector out(d);
  
  for (int j = 0; j < input.size(); j++) {
    for (int i = 0; i < input[0].size(); i++) {
      out[TDIndex(input.size(), input[0].size(), 2, j, i, 0)] =
        input[j][i].first;
      out[TDIndex(input.size(), input[0].size(), 2, j, i, 1)] =
        input[j][i].second;
      
      // // REVIEW: Why is this needed for discrete landscapes? -JCB
      // if (input[j][i].first == INT_MAX)
      //   out[TDIndex(input.size(), input[0].size(), 2, j, i, 0)] = R_PosInf;
      // if (input[j][i].first == INT_MIN)
      //   out[TDIndex(input.size(), input[0].size(), 2, j, i, 0)] = R_NegInf;
    }
  }
  
  return out;
}

// Section: Class

//' @name PersistenceLandscape
//' @aliases Rcpp_PersistenceLandscape-class
//' @title Exported C++ Class 'PersistenceLandscape'
//' @description Export, and create and manipulate objects of, the
//'   PersistenceLandscape' C++ class.
//' @details The C++ class 'PersistenceLandscape' is exposed as the S4 class
//'   'Rcpp_PersistenceLandscape' via the `RCPP_MODULE()` macro provided by
//'   **[Rcpp][Rcpp::Rcpp-package]**. See
//'   <https://github.com/r-pkg-examples/rcpp-modules-student> for an
//'   introduction. New objects should be created from persistence data
//'   (diagrams) using [landscape()].
//' @field new Constructor.
//' @field exact Representation of the underlying PL.
//' @field min_x Infimum (left endpoint) of the support of a discrete PL.
//' @field max_x Supremum (right endpoint) of the support of a discrete PL.
//' @field dx Resolution of the representation of a discrete PL.
//' @example inst/examples/ex-PersistenceLandscape.R
//' @export
class PersistenceLandscape {
  
private:
  
  // fields used to control other behaviors
  // initialize inline if uniform across constructors
  // https://stackoverflow.com/a/11594963/4556798
  bool exact;
  double min_x;
  double max_x;
  double dx;
  
public:
  
  // Subsection: Constructors
  
  // REVIEW: I don't understand this. -JCB
  PersistenceLandscape(){}
  
  // // construct an exact persistence landscape
  // PersistenceLandscape(const NumericMatrix &x);
  
  // // construct a discrete persistence landscape
  // PersistenceLandscape(
  //   const NumericMatrix &x,
  //   // imposed range of landscape abscissa
  //   double min_x = 0, double max_x = 1, double dx = 0.001
  // );
  
  // construct a persistence landscape, either exact or discrete
  PersistenceLandscape(
    const NumericMatrix &x,
    // type of landscape encoding
    bool exact = true,
    // imposed range of landscape abscissa (discrete only)
    double min_x = 0, double max_x = 1, double dx = 0.001
  );
  
  // // construct a discrete landscape from an exact one
  // PersistenceLandscape(
  //   const PersistenceLandscape &pl,
  //   double min_x = 0, double max_x = 1, double dx = 0.001
  // );
  
  // Subsection: Internals
  
  // landscape data
  std::vector<std::vector<std::pair<double, double>>> land;
  // level counter
  size_t size() const { return this->land.size(); }
  
  double landscapeValue(
      unsigned level,
      double x
  ) const;
  
  // getters
  bool isExact() const { return exact; }
  double xMin() const { return min_x; }
  double xMax() const { return max_x; }
  double xBy() const { return dx; }
  
  SEXP getInternal() {
    if (exact)
      return wrap(getInternalExact(land));
    else
      return wrap(getInternalDiscrete(land));
  }
  
  NumericVector discretize() {
    PersistenceLandscape disc;
    if (exact)
      disc = discretizeExactLandscape(*this, min_x, max_x, dx);
    else
      disc = *this;
    return wrap(getInternalDiscrete(disc.land));
  }
  
  // Subsection: Friendzone
  
  friend PersistenceLandscape discretizeExactLandscape(
    const PersistenceLandscape &pl,
    double min_x, double max_x, double dx
  );
  
  friend PersistenceLandscape discretizeLandscape(
      const PersistenceLandscape &pl,
      double min_x, double max_x, double dx
  );
  
  // TODO: Find more uses for this or replace with addition only. -JCB
  // This is a general algorithm to perform linear operations on persistence
  // landscapes. It perform it by doing operations on landscape points.
  friend PersistenceLandscape operateExactLandscapes(
      const PersistenceLandscape &pl1,
      const PersistenceLandscape &pl2,
      double (*oper)(double, double)
  );
  friend PersistenceLandscape addExactLandscapes(
      const PersistenceLandscape &pl1,
      const PersistenceLandscape &pl2
  ) {
    return operateExactLandscapes(pl1, pl2, addOp);
  }
  friend PersistenceLandscape addDiscreteLandscapes(
      const PersistenceLandscape &pl1,
      const PersistenceLandscape &pl2
  );
  
  // friend PersistenceLandscape operator*(
  //     const PersistenceLandscape &first,
  //     double con) {
  //   return first.scale(con);
  // }
  // friend PersistenceLandscape operator*(
  //     double con,
  //     const PersistenceLandscape &first) {
  //   return first.scale(con);
  // }
  // friend PersistenceLandscape operator+(
  //     const PersistenceLandscape &pl1,
  //     const PersistenceLandscape &pl2) {
  //   return addExactLandscapes(pl1, pl2);
  // }
  
  friend double innerProductExactLandscapes(
      const PersistenceLandscape &pl1,
      const PersistenceLandscape &pl2
  );
  friend double innerProductDiscreteLandscapes(
      const PersistenceLandscape &pl1,
      const PersistenceLandscape &pl2
  );
  
  // Subsection: Functions
  
  double minimum(unsigned level) const;
  
  double maximum(unsigned level) const;
  
  PersistenceLandscape scale(double x) const;
  
  PersistenceLandscape add(const PersistenceLandscape &other);
  
  PersistenceLandscape abs();
  
  double inner(const PersistenceLandscape &other);
  
  double computeMoment(
      unsigned p,
      double center,
      unsigned level
  ) const;
  
};

// Section: Constructors

// TODO: Write separate exact and discrete constructors. This may require
// writing submodules called by both gnostic and agnostic constructors.

// PersistenceLandscape::PersistenceLandscape(
//   const NumericMatrix &x
// ) : exact(true), min_x(R_NegInf), max_x(R_PosInf), dx(0) {
//   // move conditional part of `PersistenceLandscape()` here
// }

// PersistenceLandscape::PersistenceLandscape(
//   const NumericMatrix &x,
//   double min_x, double max_x, double dx
// ) : exact(false), min_x(min_x), max_x(max_x), dx(dx) {
//   // move conditional part of `PersistenceLandscape()` here
// }

PersistenceLandscape::PersistenceLandscape(
  const NumericMatrix &x,
  bool exact,
  double min_x, double max_x, double dx
) : exact(exact), min_x(min_x), max_x(max_x), dx(dx) {
  
  std::vector<std::pair<double, double>> pd = persistencePairs(x);
  
  if (exact) {
    
    std::vector<std::pair<double, double>> pds;
    pds.insert(pds.begin(), pd.begin(), pd.end());
    std::sort(pds.begin(), pds.end(), compareFeatures);
    
    // TODO: Choose a shorter name than `characteristicPoints()`.
    std::vector<std::pair<double, double>> characteristicPoints(pds.size());
    
    for (size_t i = 0; i != pds.size(); ++i) {
      characteristicPoints[i] =
        std::make_pair((pds[i].first + pds[i].second) / 2.0,
                       (pds[i].second - pds[i].first) / 2.0);
    }
    
    while (! characteristicPoints.empty()) {
      std::vector<std::pair<double, double>> lambda_n;
      lambda_n.push_back(std::make_pair(INT_MIN, 0.));
      lambda_n.push_back(std::make_pair(birth(characteristicPoints[0]), 0.));
      lambda_n.push_back(characteristicPoints[0]);
      
      int i = 1;
      std::vector<std::pair<double, double>> newCharacteristicPoints;
      while (i < characteristicPoints.size()) {
        int j = 1;
        if ((birth(characteristicPoints[i]) >=
            birth(lambda_n[lambda_n.size() - 1])) &&
            (death(characteristicPoints[i]) >
            death(lambda_n[lambda_n.size() - 1]))) {
          if (birth(characteristicPoints[i]) <
            death(lambda_n[lambda_n.size() - 1])) {
            std::pair<double, double> point =
              std::make_pair((birth(characteristicPoints[i]) +
              death(lambda_n[lambda_n.size() - 1])) /
                2.,
                (death(lambda_n[lambda_n.size() - 1]) -
                  birth(characteristicPoints[i])) /
                    2.);
            lambda_n.push_back(point);
            
            while ((i + j < characteristicPoints.size()) &&
                   (almostEqual(birth(point),
                                birth(characteristicPoints[i + j]))) &&
                                  (death(point) <=
                                  death(characteristicPoints[i + j]))) {
              newCharacteristicPoints.push_back(characteristicPoints[i + j]);
              
              ++j;
            }
            
            newCharacteristicPoints.push_back(point);
            
            while ((i + j < characteristicPoints.size()) &&
                   (birth(point) <= birth(characteristicPoints[i + j])) &&
                   (death(point) >= death(characteristicPoints[i + j]))) {
              newCharacteristicPoints.push_back(characteristicPoints[i + j]);
              ++j;
            }
            
          } else {
            lambda_n.push_back(
              std::make_pair(death(lambda_n[lambda_n.size() - 1]), 0.));
            lambda_n.push_back(
              std::make_pair(birth(characteristicPoints[i]), 0.));
          }
          lambda_n.push_back(characteristicPoints[i]);
        } else {
          newCharacteristicPoints.push_back(characteristicPoints[i]);
        }
        i = i + j;
      }
      lambda_n.push_back(
        std::make_pair(death(lambda_n[lambda_n.size() - 1]), 0.));
      lambda_n.push_back(std::make_pair(INT_MAX, 0.));
      
      // CHANGE
      characteristicPoints = newCharacteristicPoints;
      // characteristicPoints.swap(newCharacteristicPoints);
      
      lambda_n.erase(std::unique(lambda_n.begin(), lambda_n.end()),
                     lambda_n.end());
      
      this->land.push_back(lambda_n);
    }
    
  } else {
    
    size_t numberOfBins = 2 * ((max_x - min_x) / dx) + 1;
    
    // The first element of a pair `std::pair< double, std::vector<double> >`
    // is an x-value. The second element is a vector of values of landscapes.
    std::vector<std::pair<double, std::vector<double>>>
      criticalValuesOnPointsOfGrid(numberOfBins);
    
    // Filling up the bins:
    
    // Now, the idea is to iterate on `this->land[lambda-1]` and use only points
    // over there. The problem is at the very beginning, when there is nothing
    // in `this->land`. That is why over here, we make a fake `this->land[0]`.
    // It will be later deleted before moving on.
    std::vector<std::pair<double, double>> aa;
    double x = min_x;
    for (size_t i = 0; i != numberOfBins; ++i) {
      std::vector<double> v;
      std::pair<double, std::vector<double>> p = std::make_pair(x, v);
      criticalValuesOnPointsOfGrid[i] = p;
      aa.push_back(std::make_pair(x, 0.));
      x += 0.5 * dx;
    }
    
    // For every persistent interval, sample on the grid.
    for (size_t intervalNo = 0; intervalNo != pd.size(); ++intervalNo) {
      size_t beginn = 0;
      
      while (beginn < criticalValuesOnPointsOfGrid.size()) {
        if (fabs(criticalValuesOnPointsOfGrid[beginn].first >
                   pd[intervalNo].first) &&
                   fabs(criticalValuesOnPointsOfGrid[beginn].first <
                     pd[intervalNo].second)) {
          criticalValuesOnPointsOfGrid[beginn].second.push_back(
              std::min(fabs(criticalValuesOnPointsOfGrid[beginn].first -
                pd[intervalNo].first),
                fabs(criticalValuesOnPointsOfGrid[beginn].first -
                  pd[intervalNo].second)));
        } else
          criticalValuesOnPointsOfGrid[beginn].second.push_back(0.0);
        
        ++beginn;
      }
    }
    
    // Now, the basic structure is created. We need to translate it to a
    // persistence landscape data structure. To do so, first we need to sort all
    // the vectors in `criticalValuesOnPointsOfGrid[i].second`.
    size_t maxNonzeroLambda = 0;
    for (size_t i = 0; i != criticalValuesOnPointsOfGrid.size(); ++i) {
      std::sort(criticalValuesOnPointsOfGrid[i].second.begin(),
                criticalValuesOnPointsOfGrid[i].second.end(),
                std::greater<double>());
      if (criticalValuesOnPointsOfGrid[i].second.size() > maxNonzeroLambda) {
        maxNonzeroLambda = criticalValuesOnPointsOfGrid[i].second.size();
      }
    }
    
    // Initialize to zero
    this->land.resize(maxNonzeroLambda, aa);
    
    // Add values
    for (unsigned int i = 0; i < criticalValuesOnPointsOfGrid.size(); i++) {
      for (size_t lambda = 0;
           lambda < criticalValuesOnPointsOfGrid[i].second.size(); ++lambda) {
        this->land[lambda][i] =
          std::make_pair(criticalValuesOnPointsOfGrid[i].first,
                         criticalValuesOnPointsOfGrid[i].second[lambda]);
      }
    }
    
  }

}

// Section: Converters

// this is O(log(n)) algorithm, where n is number of points in `this->land`.
double PersistenceLandscape::landscapeValue(
    unsigned level,
    double x
) const {
  
  // undefined levels are uniformly zero functions
  if (level > this->land.size())
    return 0;
  
  // the points in `this->land[level]` are ordered by x coordinate, so we can
  // find the point by using bisection (excluding infinite 'endpoints'):
  unsigned coordBegin = 1;
  unsigned coordEnd = this->land[level].size() - 2;
  
  // if `x` is outside the support of this landscape level, then the value is 0
  if (x <= this->land[level][coordBegin].first)
    return 0;
  if (x >= this->land[level][coordEnd].first)
    return 0;
  
  while (coordBegin + 1 != coordEnd) {
    unsigned coordNew = (unsigned)floor((coordEnd + coordBegin) / 2.0);
    
    if (this->land[level][coordNew].first <= x) {
      coordBegin = coordNew;
      if (this->land[level][coordNew].first == x)
        return this->land[level][coordNew].second;
    } else {
      coordEnd = coordNew;
    }
  }
  
  return interpolatePoint(this->land[level][coordBegin],
                          this->land[level][coordEnd], x);
}

PersistenceLandscape discretizeExactLandscape(
    const PersistenceLandscape &pl,
    double min_x, double max_x, double dx
) {
  
  PersistenceLandscape result;
  result.exact = false;
  result.min_x = min_x;
  result.max_x = max_x;
  result.dx = dx;
  
  std::vector<std::vector<std::pair<double, double>>> land_disc;
  std::pair<double, double> currentPoint;
  
  if (pl.land[0][1].first < min_x ||
      pl.land[0][pl.land[0].size() - 2].first > max_x)
    warning("This landscape extends beyond [ `min_x`, `max_x` ].");
  
  for (unsigned int i = 0; i < pl.land.size(); i++) {
    
    auto level = pl.land[i];
    std::vector<std::pair<double, double>> level_disc;
    
    // Start at the level value at `min_x`.
    double x_buff = min_x;
    double y_buff = pl.landscapeValue(i, min_x);
    currentPoint = std::make_pair(x_buff, y_buff);
    level_disc.push_back(currentPoint);
    
    // Iterate over finite critical points...
    for (int j = 1; j < level.size() - 1; j++) {
      
      // Skip to the critical point rightward of `currentPoint` for which the
      // next critical point is just leftward of `currentPoint + dx`.
      if (level[j].first <= x_buff || level[j + 1].first < x_buff + dx)
        continue;
      
      // If the next critical point is at least `dx` rightward, then increment
      // linearly to it; else, compute and increment to the next level value.
      if (level[j].first < x_buff + dx) {
        
        x_buff += dx;
        y_buff = pl.landscapeValue(i, x_buff + dx);
        currentPoint = std::make_pair(x_buff, y_buff);
        level_disc.push_back(currentPoint);
        
      } else {
        std::pair<double, double> nextPoint = level[j];
        
        // If change in x, increment linearly from current to just before next.
        if (nextPoint.first != currentPoint.first) {
          double delta_x = nextPoint.first - currentPoint.first;
          double delta_y = nextPoint.second - currentPoint.second;
          double slope = delta_y / delta_x;
          
          int n_incr = std::floor((std::min(nextPoint.first, max_x) - x_buff) /
                                  dx);
          for (int k = 0; k < n_incr; k++) {
            x_buff += dx;
            y_buff += dx * slope;
            level_disc.push_back(std::make_pair(x_buff, y_buff));
          }
        }
        
        currentPoint = std::make_pair(x_buff, y_buff);
      }
    }
    
    // If `max_x` has not been reached, then increment along zero y values.
    if (x_buff + dx < max_x + epsi) {
      int n_incr = std::floor((max_x + epsi - x_buff) / dx);
      y_buff = 0;
      for (int k = 0; k < n_incr; k++) {
        x_buff += dx;
        level_disc.push_back(std::make_pair(x_buff, y_buff));
      }
    }
    
    land_disc.push_back(level_disc);
  }
  
  result.land = land_disc;
  return result;
}

PersistenceLandscape discretizeLandscape(
    const PersistenceLandscape &pl,
    double min_x, double max_x, double dx
) {
  
  PersistenceLandscape out;
  
  if (pl.exact) {
    out = discretizeExactLandscape(pl, min_x, max_x, dx);
  } else {
    warning("Can not yet re-discretize a discrete PL.");
    out = pl;
  }
  
  return out;
}

// Section: Operations

PersistenceLandscape PersistenceLandscape::scale(
    double x
) const {
  
  std::vector<std::vector<std::pair<double, double>>> result(this->land.size());
  
  for (size_t lev = 0; lev != this->land.size(); ++lev) {
    
    std::vector<std::pair<double, double>> lambda_lev(this->land[lev].size());
    
    for (size_t i = 0; i != this->land[lev].size(); ++i) {
      lambda_lev[i] = std::make_pair(this->land[lev][i].first,
                                     x * this->land[lev][i].second);
    }
    result[lev] = lambda_lev;
  }
  
  // uses empty constructor `PersistenceLandscape(){}`
  PersistenceLandscape product;
  product.exact = this->exact;
  product.min_x = this->min_x;
  product.max_x = this->max_x;
  product.dx = this->dx;
  // CHANGE
  // product.land = result;
  product.land.swap(result);
  
  return product;
}

PersistenceLandscape operateExactLandscapes(
    const PersistenceLandscape &pl1,
    const PersistenceLandscape &pl2,
    double (*oper)(double, double)
) {
  if (! pl1.exact || ! pl2.exact)
    stop("`operateExactLandscapes()` requires both landscapes to be exact.");
  
  PersistenceLandscape result;
  result.exact = true;
  result.min_x = min(pl1.min_x, pl2.min_x);
  result.max_x = max(pl1.max_x, pl2.max_x);
  std::vector<std::vector<std::pair<double, double>>> land(
      std::max(pl1.land.size(), pl2.land.size()));
  result.land = land;
  
  for (size_t i = 0; i != std::min(pl1.land.size(), pl2.land.size()); ++i) {
    
    std::vector<std::pair<double, double>> lambda_n;
    
    int p = 0;
    int q = 0;
    while ((p + 1 < pl1.land[i].size()) && (q + 1 < pl2.land[i].size())) {
      
      if (pl1.land[i][p].first < pl2.land[i][q].first) {
        lambda_n.push_back(std::make_pair(
            pl1.land[i][p].first,
            oper(pl1.land[i][p].second,
                 interpolatePoint(pl2.land[i][q - 1], pl2.land[i][q],
                                  pl1.land[i][p].first))));
        ++p;
        continue;
      }
      if (pl1.land[i][p].first > pl2.land[i][q].first) {
        lambda_n.push_back(std::make_pair(
            pl2.land[i][q].first,
            oper(interpolatePoint(pl1.land[i][p], pl1.land[i][p - 1],
                                  pl2.land[i][q].first),
                                  pl2.land[i][q].second)));
        ++q;
        continue;
      }
      if (pl1.land[i][p].first == pl2.land[i][q].first) {
        lambda_n.push_back(std::make_pair(
            pl2.land[i][q].first,
            oper(pl1.land[i][p].second, pl2.land[i][q].second)));
        ++p;
        ++q;
      }
    }
    
    while ((p + 1 < pl1.land[i].size()) && (q + 1 >= pl2.land[i].size())) {
      lambda_n.push_back(std::make_pair(pl1.land[i][p].first,
                                        oper(pl1.land[i][p].second, 0)));
      ++p;
    }
    while ((p + 1 >= pl1.land[i].size()) && (q + 1 < pl2.land[i].size())) {
      lambda_n.push_back(std::make_pair(pl2.land[i][q].first,
                                        oper(0, pl2.land[i][q].second)));
      ++q;
    }
    
    // REVIEW: Try to prevent operations from infinitizing endpoints. -JCB
    if (pl1.land[i][p].first == pl2.land[i][q].first) {
      if (pl2.land[i][q].first == INT_MAX) {
        lambda_n.push_back(std::make_pair(INT_MAX, 0.));
      } else {
        lambda_n.push_back(std::make_pair(pl1.land[i][p].first,
                                          oper(pl1.land[i][p].second,
                                               pl2.land[i][q].second)));
      }
    } else {
      // REVIEW: Is there a better option in this case? -JCB
      lambda_n.push_back(std::make_pair(INT_MAX, 0.));
    }
    // CHANGE
    // result.land[i] = lambda_n;
    result.land[i].swap(lambda_n);
  }
  
  if (pl1.land.size() > std::min(pl1.land.size(), pl2.land.size())) {
    for (size_t i = std::min(pl1.land.size(), pl2.land.size());
         i != std::max(pl1.land.size(), pl2.land.size()); ++i) {
      std::vector<std::pair<double, double>> lambda_n(pl1.land[i]);
      for (size_t nr = 0; nr != pl1.land[i].size(); ++nr) {
        lambda_n[nr] = std::make_pair(pl1.land[i][nr].first,
                                      oper(pl1.land[i][nr].second, 0));
      }
      // CHANGE
      // result.land[i] = lambda_n;
      result.land[i].swap(lambda_n);
    }
  }
  
  if (pl2.land.size() > std::min(pl1.land.size(), pl2.land.size())) {
    for (size_t i = std::min(pl1.land.size(), pl2.land.size());
         i != std::max(pl1.land.size(), pl2.land.size()); ++i) {
      std::vector<std::pair<double, double>> lambda_n(pl2.land[i]);
      for (size_t nr = 0; nr != pl2.land[i].size(); ++nr) {
        lambda_n[nr] = std::make_pair(pl2.land[i][nr].first,
                                      oper(0, pl2.land[i][nr].second));
      }
      // CHANGE
      // result.land[i] = lambda_n;
      result.land[i].swap(lambda_n);
    }
  }
  
  return result;
}

// REVIEW: Assume only that the `dx` are equal and that they divide the
// difference between the `min_x`. -JCB
PersistenceLandscape addDiscreteLandscapes(
    const PersistenceLandscape &pl1,
    const PersistenceLandscape &pl2
) {
  if (pl1.exact || pl2.exact)
    stop("`addDiscreteLandscapes()` requires two discrete PLs.");
  if (! almostEqual(pl1.dx, pl2.dx))
    stop("`addDiscreteLandscapes()` requires PLs with same `dx`.");
  
  PersistenceLandscape sum;
  sum.exact = false;
  sum.min_x = min(pl1.min_x, pl2.min_x);
  sum.max_x = max(pl1.max_x, pl2.max_x);
  sum.dx = pl1.dx;
  std::vector<std::vector<std::pair<double, double>>> land_sum;
  
  // insert sums of shared levels
  int min_level = std::min(pl1.land.size(), pl2.land.size());
  for (int i = 0; i < min_level; i++) {
    
    std::vector<std::pair<double, double>> level_sum;
    
    double lead_x = pl1.land[i][1].first - pl2.land[i][1].first;
    double lag_x = pl1.land[i][pl1.land[i].size()].first -
      pl2.land[i][pl2.land[i].size()].first;
    int lead_index = std::round(lead_x / pl1.dx);
    int lag_index = std::round(lag_x / pl1.dx);
    
    // insert preceding indices from whichever landscape has them
    if (lead_index < 0) {
      // `pl1` begins first
      for (int j = 0; j < -lead_index; j++)
        level_sum.push_back(pl1.land[i][j]);
    } else if (lead_index > 0) {
      // `pl2` begins first
      for (int j = 0; j < lead_index; j++)
        level_sum.push_back(std::make_pair(
          pl1.land[i][0].first - (lead_index - j) * pl1.dx,
          pl2.land[i][j].second
        ));
    }
    
    // insert sums of shared indices
    int min_index = std::max(0, lead_index);
    int max_index = std::min(pl1.land[i].size() + lead_index,
                             pl2.land[i].size());
    for (int j = min_index; j < max_index; j++)
      level_sum.push_back(std::make_pair(
        pl1.land[i][j - lead_index].first,
        pl1.land[i][j].second + pl2.land[i][j].second
      ));
    
    // insert succeeding indices from whichever landscape has them
    if (lag_index < 0) {
      // `pl1` ends first
      for (int j = max_index; j < pl2.land[i].size(); j++)
        level_sum.push_back(std::make_pair(
          pl1.land[i][pl1.land[i].size()].first + (j - max_index + 1) * pl1.dx,
          pl2.land[i][j].second
        ));
    } else if (lag_index > 0) {
      // `pl2` ends first
      for (int j = max_index - lead_index; j < pl1.land[i].size(); j++)
        level_sum.push_back(pl1.land[i][j]);
    }
    
    land_sum.push_back(level_sum);
  }
  
  // insert succeeding levels from whichever landscape has them
  for (; min_level < pl1.land.size(); min_level++)
    land_sum.push_back(pl1.land[min_level]);
  for (; min_level < pl2.land.size(); min_level++)
    land_sum.push_back(pl2.land[min_level]);
  
  sum.land = land_sum;
  return sum;
}

PersistenceLandscape PersistenceLandscape::add(
    const PersistenceLandscape &other
) {
  
  // uses empty constructor `PersistenceLandscape(){}`
  PersistenceLandscape sum;
  
  // both landscapes are exact
  if (this->exact && other.exact)
    sum = addExactLandscapes(*this, other);
  
  // both landscapes are discrete
  else if (! this->exact && ! other.exact)
    sum = addDiscreteLandscapes(*this, other);
  
  // one landscape is exact, the other discrete
  else if (this->exact) {
    // only this landscape is exact
    PersistenceLandscape conversion1 = discretizeExactLandscape(
      *this,
      other.min_x, other.max_x, other.dx
    );
    sum = addDiscreteLandscapes(conversion1, other);
  } else {
    // only other landscape is exact
    PersistenceLandscape conversion2 = discretizeExactLandscape(
      other,
      this->min_x, this->max_x, this->dx
    );
    sum = addDiscreteLandscapes(*this, conversion2);
  }
  
  return sum;
}

double locateIntermediateRoot(
    std::pair<double, double> p1,
    std::pair<double, double> p2
) {
  
  if (p1.first == p2.first)
    return p1.first;
  
  // throw error if segment does not cross abscissa
  if (p1.second * p2.second > 0) {
    std::ostringstream errMessage;
    errMessage << "Arguments to `locateIntermediateRoot()` were ("
    << p1.first << "," << p1.second << ") and ("
    << p2.first << "," << p2.second
    << "). The segment between those points does not cross the abscissa.";
    std::string errMessageStr = errMessage.str();
    const char *err = errMessageStr.c_str();
    throw(err);
  }
  
  // assume that `p1.first <= x <= p2.first`
  double a = (p2.second - p1.second) / (p2.first - p1.first);
  double b = p1.second - a * p1.first;
  // cerr << "Line crossing points : (" << p1.first << "," << p1.second << ")
  // oraz (" << p2.first << "," << p2.second << ") : \n"; cerr << "a : " << a <<
  // " , b : " << b << " , x : " << x << endl;
  return -b / a;
}

PersistenceLandscape PersistenceLandscape::abs() {
  
  PersistenceLandscape result;
  result.exact = this->exact;
  result.min_x = this->min_x;
  result.max_x = this->max_x;
  result.dx = this->dx;
  
  for (size_t level = 0; level != this->land.size(); ++level) {
    
    std::vector<std::pair<double, double>> lambda_n;
    
    // REVIEW: Try to prevent operations from infinitizing endpoints. -JCB
    // if (this->land[level][0].first == INT_MIN) {
    //   lambda_n.push_back(std::make_pair(INT_MIN, 0.));
    // } else {
    //   lambda_n.push_back(std::make_pair(this->land[level][0].first,
    //                                     fabs(this->land[level][0].second)));
    // }
    lambda_n.push_back(std::make_pair(
        this->land[level][0].first,
        fabs(this->land[level][0].second)
    ));
    
    for (size_t i = 1; i != this->land[level].size(); ++i) {
      // if a line segment between this->land[level][i-1] and
      // this->land[level][i] crosses the x-axis, then we have to add one
      // landscape point to result
      if ((this->land[level][i - 1].second) *
          (this->land[level][i].second) < 0) {
        
        double zero = locateIntermediateRoot(
          this->land[level][i - 1], this->land[level][i]);
        lambda_n.push_back(std::make_pair(zero, 0.));
        
        lambda_n.push_back(std::make_pair(
            this->land[level][i].first,
            fabs(this->land[level][i].second)
        ));
        
      } else {
        
        lambda_n.push_back(std::make_pair(
            this->land[level][i].first,
            fabs(this->land[level][i].second)
        ));
        
      }
    }
    result.land.push_back(lambda_n);
  }
  
  return result;
}

double innerProductExactLandscapes(
    const PersistenceLandscape &pl1,
    const PersistenceLandscape &pl2
) {
  double result = 0;
  
  for (size_t level = 0; level != std::min(pl1.size(), pl2.size()); ++level) {
    if (pl1.land[level].size() * pl2.land[level].size() == 0)
      continue;
    
    // endpoints of the interval on which we will compute the inner product of
    // two locally linear functions:
    double x1 = INT_MIN;
    double x2;
    if (pl1.land[level][1].first < pl2.land[level][1].first) {
      x2 = pl1.land[level][1].first;
    } else {
      x2 = pl2.land[level][1].first;
    }
    
    // iterators for the landscapes pl1 and pl2
    size_t l1It = 0;
    size_t l2It = 0;
    
    while ((l1It < pl1.land[level].size() - 1) &&
           (l2It < pl2.land[level].size() - 1)) {
      // compute the value of a inner product on an interval [x1,x2]
      
      double a, b, c, d;
      
      a = (pl1.land[level][l1It + 1].second - pl1.land[level][l1It].second) /
        (pl1.land[level][l1It + 1].first - pl1.land[level][l1It].first);
      b = pl1.land[level][l1It].second - a * pl1.land[level][l1It].first;
      c = (pl2.land[level][l2It + 1].second - pl2.land[level][l2It].second) /
        (pl2.land[level][l2It + 1].first - pl2.land[level][l2It].first);
      d = pl2.land[level][l2It].second - c * pl2.land[level][l2It].first;
      
      double contributionFromThisPart =
        (a * c * x2 * x2 * x2 / 3 + (a * d + b * c) * x2 * x2 / 2 +
        b * d * x2) -
        (a * c * x1 * x1 * x1 / 3 + (a * d + b * c) * x1 * x1 / 2 +
        b * d * x1);
      
      result += contributionFromThisPart;
      
      // we have two intervals in which functions are constant:
      //[pl1.land[level][l1It].first , pl1.land[level][l1It+1].first]
      // and
      //[pl2.land[level][l2It].first , pl2.land[level][l2It+1].first]
      // We also have an interval [x1,x2]. Since the intervals in the landscapes
      // cover the whole R, then it is clear that x2 is either
      // pl1.land[level][l1It+1].first of pl2.land[level][l2It+1].first or both.
      // Lets test it.
      if (x2 == pl1.land[level][l1It + 1].first) {
        if (x2 == pl2.land[level][l2It + 1].first) {
          // in this case, we increment both:
          ++l2It;
        }
        ++l1It;
      } else {
        // in this case we increment l2It
        ++l2It;
      }
      // Now, we shift x1 and x2:
      x1 = x2;
      if (pl1.land[level][l1It + 1].first < pl2.land[level][l2It + 1].first) {
        x2 = pl1.land[level][l1It + 1].first;
      } else {
        x2 = pl2.land[level][l2It + 1].first;
      }
    }
  }
  
  return result;
}

double innerProductDiscreteLandscapes(
    const PersistenceLandscape &pl1,
    const PersistenceLandscape &pl2
) {
  
  int min_level = std::min(pl1.land.size(), pl2.land.size());
  double integral_buffer = 0;
  
  for (int i = 0; i < min_level; i++) {
    int min_index = std::min(pl1.land[i].size(), pl2.land[i].size());
    for (int j = 0; j < min_index; j++)
      integral_buffer += pl1.land[i][j].second * pl2.land[i][j].second;
  }
  
  return integral_buffer * pl1.dx;
}

double PersistenceLandscape::inner(
    const PersistenceLandscape &other
) {
  
  double scalar;
  
  // both landscapes are exact
  if (this->exact && other.exact)
    scalar = innerProductExactLandscapes(*this, other);
  
  // both landscapes are discrete
  else if (! this->exact && ! other.exact)
    scalar = innerProductDiscreteLandscapes(*this, other);
  
  // one landscape is exact, the other discrete
  else if (this->exact) {
    // only this landscape is exact
    PersistenceLandscape conversion1 = discretizeExactLandscape(
      *this,
      other.min_x, other.max_x, other.dx
    );
    scalar = innerProductDiscreteLandscapes(conversion1, other);
  } else {
    // only other landscape is exact
    PersistenceLandscape conversion2 = discretizeExactLandscape(
      other,
      this->min_x, this->max_x, this->dx
    );
    scalar = innerProductDiscreteLandscapes(*this, conversion2);
  }
  
  return scalar;
}

// This function finds the minimum value at a given level.
double PersistenceLandscape::minimum(
    unsigned level
) const {
  if (level < 0 || level >= this->land.size())
    return NA_REAL;
  if (this->land.size() < level)
    return 0;
  double level_min = INT_MAX;
  for (size_t i = 0; i != this->land[level].size(); ++i) {
    if (this->land[level][i].second < level_min)
      level_min = this->land[level][i].second;
  }
  return level_min;
}

// This function finds the maximum value at a given level.
double PersistenceLandscape::maximum(
    unsigned level
) const {
  if (level < 0 || level >= this->land.size())
    return NA_REAL;
  if (this->land.size() < level)
    return 0;
  double level_max = INT_MIN;
  for (size_t i = 0; i != this->land[level].size(); ++i) {
    if (this->land[level][i].second > level_max)
      level_max = this->land[level][i].second;
  }
  return level_max;
}

// This function computes the n^th moment of a given level.
double PersistenceLandscape::computeMoment(
    unsigned p,
    double center,
    unsigned level
) const {
  if (p < 1)
    throw("Cannot compute p^th moment for p < 1.\n");
  if (level < 0 || level > this->land.size())
    return NA_REAL;
  
  double result = 0;
  
  if (this->land.size() > level) {
    for (size_t i = 2; i != this->land[level].size() - 1; ++i) {
      if (this->land[level][i].first - this->land[level][i - 1].first == 0)
        continue;
      // Between `this->land[level][i]` and `this->land[level][i-1]`, the
      // `lambda_level` is of the form a x + b. First we need to find a and b.
      double a =
        (this->land[level][i].second - this->land[level][i - 1].second) /
          (this->land[level][i].first - this->land[level][i - 1].first);
      double b =
        this->land[level][i - 1].second - a * this->land[level][i - 1].first;
      
      double x1 = this->land[level][i - 1].first;
      double x2 = this->land[level][i].first;
      
      // double first =
      // b*(pow((x2-center),(double)(p+1))/(p+1)-
      // pow((x1-center),(double)(p+1))/(p+1));
      // double second = a/(p+1)*((x2*pow((x2-center),(double)(p+1))) -
      // (x1*pow((x1-center),(double)(p+1))) )
      //              +
      //              a/(p+1)*( pow((x2-center),(double)(p+2))/(p+2) -
      //              pow((x1-center),(double)(p+2))/(p+2) );
      // result += first;
      // result += second;
      
      double first = a / (p + 2) *
        (pow((x2 - center), (double)(p + 2)) -
        pow((x1 - center), (double)(p + 2)));
      double second = center / (p + 1) *
        (pow((x2 - center), (double)(p + 1)) -
        pow((x1 - center), (double)(p + 1)));
      double third = b / (p + 1) *
        (pow((x2 - center), (double)(p + 1)) -
        pow((x1 - center), (double)(p + 1)));
      
      result += first + second + third;
    }
  }
  
  return result;
}

// Section: List operations

PersistenceLandscape PLsum(List pl_list) {
  
  PersistenceLandscape sum = as<PersistenceLandscape>(pl_list[0]);
  
  for (int i = 1; i < pl_list.size(); i++) {
    sum = sum.add(as<PersistenceLandscape>(pl_list[i]));
  }
  
  return sum;
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
  
  ;
  
  Rcpp::function("PLsum", &PLsum);
}
