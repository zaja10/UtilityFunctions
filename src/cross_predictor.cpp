#include <Rcpp.h>
using namespace Rcpp;

//' Calculate Gametic Variance for a Cross (Biparental)
//'
//' Computes the expected genetic variance in the F2 population derived from
//' crossing two inbred parents, accounting for linkage disequilibrium.
//'
//' @param p1 NumericVector. Genotype of Parent 1 (0 or 2).
//' @param p2 NumericVector. Genotype of Parent 2 (0 or 2).
//' @param effects NumericVector. Marker effects (a).
//' @param map NumericVector. Genetic map positions in Morgans.
//' @param max_dist double. Maximum distance (in Morgans) to calculate
//covariance. '        Markers further apart are assumed unlinked. Default 3.0.
//'
//' @return double. The gametic variance.
//'
// [[Rcpp::export]]
double calc_cross_variance_cpp(NumericVector p1, NumericVector p2,
                               NumericVector effects, NumericVector map,
                               double max_dist = 3.0) {
  int m = p1.size();
  if (p2.size() != m || effects.size() != m || map.size() != m) {
    stop("Input vectors must have the same length.");
  }

  // 1. Identify segregating loci
  std::vector<int> seg_idx;
  seg_idx.reserve(m); // Optimization: Reserve memory

  for (int i = 0; i < m; i++) {
    if (std::abs(p1[i] - p2[i]) > 0.5) {
      seg_idx.push_back(i);
    }
  }

  int n_seg = seg_idx.size();
  if (n_seg == 0)
    return 0.0;

  double variance = 0.0;

  // 2. Variance Component (Additive)
  for (int idx : seg_idx) {
    double a = effects[idx];
    variance += 0.5 * a * a;
  }

  // 3. Covariance Component (LD)
  // Optimization: Pre-calculate phase and only loop within max_dist

  for (int i = 0; i < n_seg; i++) {
    int idx_i = seg_idx[i];
    double pos_i = map[idx_i];
    double a_i = effects[idx_i];

    // Inner loop: Look ahead
    for (int j = i + 1; j < n_seg; j++) {
      int idx_j = seg_idx[j];
      double d = std::abs(map[idx_j] - pos_i);

      // OPTIMIZATION: Stop if distance exceeds threshold
      // Since map is sorted, all subsequent j will also be too far.
      // (Assumes map is ordered by position, which is standard)
      if (d > max_dist)
        break;

      // Haldane Mapping Function
      // r = 0.5 * (1 - exp(-2d))
      // linkage parameter (1 - 2r) simplifies to exp(-2d)
      double linkage = std::exp(-2.0 * d);

      // Determine Phase
      // Coupling (+): P1 has same allele (0,0 or 2,2) -> Diff < 0.5
      // Repulsion (-): P1 has diff allele (0,2 or 2,0) -> Diff > 0.5
      double phase = (std::abs(p1[idx_i] - p1[idx_j]) < 0.5) ? 1.0 : -1.0;

      double a_j = effects[idx_j];

      // Add 2 * Cov because matrix is symmetric
      variance += (linkage * phase * a_i * a_j);
      // Simplified: 2.0 * 0.5 * ... cancels out
    }
  }

  return variance;
}
