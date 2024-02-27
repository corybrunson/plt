
friend PersistenceLandscape addDiscreteLandscapes2(
  const PersistenceLandscape &pl1,
  const PersistenceLandscape &pl2
);

// REVIEW: Assume only that the `dx` are equal and that they divide the
// difference between the `min_x`. -JCB
// NOTE: Follows the convention of using the first PL's `dx` in the result.
// TODO: Add paragraph about safety conventions to documentation.
// FIXME: This function introduces `INT_MIN` & `INT_MAX` into ordinates.
// TODO: Probably best to make & use a separate function to delimit summands.
PersistenceLandscape addDiscreteLandscapes2(
    const PersistenceLandscape &pl1,
    const PersistenceLandscape &pl2
) {
  if (pl1.exact || pl2.exact)
    stop("`addDiscreteLandscapes()` requires two discrete PLs.");
  if (! almostEqual(pl1.dx, pl2.dx))
    stop("`addDiscreteLandscapes()` requires PLs with same `dx`.");
  
  PersistenceLandscape pl_sum;
  pl_sum.exact = false;
  pl_sum.min_x = min(pl1.min_x, pl2.min_x);
  pl_sum.max_x = max(pl1.max_x, pl2.max_x);
  pl_sum.dx = pl1.dx;
  std::vector<std::vector<std::pair<double, double>>> land_sum;
  
  // insert sums of shared levels
  int min_level = std::min(pl1.land.size(), pl2.land.size());
  for (int i = 0; i < min_level; i++) {
    
    std::vector<std::pair<double, double>> level_sum;
    
    double lead_x = pl1.land[i][0].first - pl2.land[i][0].first;
    double lag_x = pl1.land[i][pl1.land[i].size() - 1].first -
      pl2.land[i][pl2.land[i].size() - 1].first;
    int lead_index = std::round(lead_x / pl1.dx);
    int lag_index = std::round(lag_x / pl1.dx);
    
    // insert preceding indices from whichever landscape has them
    if (lead_index < 0) {
      // `pl1` begins first
      for (int j = 0; j < -lead_index; j++) {
        level_sum.push_back(pl1.land[i][j]);
      }
    } else if (lead_index > 0) {
      // `pl2` begins first
      for (int j = 0; j < lead_index; j++) {
        level_sum.push_back(std::make_pair(
            pl1.land[i][0].first - (lead_index - j) * pl1.dx,
            // pl2.land[i][j].first,
            pl2.land[i][j].second
        ));
      }
    }
    
    // insert sums of shared indices
    int min_index = std::max(0, lead_index);
    int max_index = std::min(pl1.land[i].size() + lead_index,
                             pl2.land[i].size());
    for (int j = min_index; j < max_index; j++) {
      level_sum.push_back(std::make_pair(
          pl1.land[i][j - lead_index].first,
          // pl2.land[i][j].first,
          pl1.land[i][j - lead_index].second + pl2.land[i][j].second
      ));
    }
    
    // insert succeeding indices from whichever landscape has them
    if (lag_index < 0) {
      // `pl1` ends first
      for (int j = 0; j < pl2.land[i].size() - max_index; j++) {
        level_sum.push_back(std::make_pair(
            pl1.land[i][pl1.land[i].size() - 1].first + (1 + j) * pl1.dx,
            // pl2.land[i][max_index + j].first,
            pl2.land[i][max_index + j].second
        ));
      }
    } else if (lag_index > 0) {
      // `pl2` ends first
      for (int j = max_index - lead_index; j < pl1.land[i].size(); j++) {
        level_sum.push_back(pl1.land[i][j]);
      }
    }
    
    land_sum.push_back(level_sum);
  }
  
  // insert succeeding levels from whichever landscape has them
  for (; min_level < pl1.land.size(); min_level++)
    land_sum.push_back(pl1.land[min_level]);
  for (; min_level < pl2.land.size(); min_level++)
    land_sum.push_back(pl2.land[min_level]);
  
  pl_sum.land = land_sum;
  return pl_sum;
}
