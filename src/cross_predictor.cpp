#include <Rcpp.h>
using namespace Rcpp;

//' Calculate Gametic Variance for a Cross (Biparental)
//' 
//' Computes the expected genetic variance in the F2 population derived from crossing
//' two inbred parents, accounting for linkage disequilibrium.
//' 
//' @param p1 NumericVector. Genotype of Parent 1 (0 or 2).
//' @param p2 NumericVector. Genotype of Parent 2 (0 or 2).
//' @param effects NumericVector. Marker effects (a).
//' @param map NumericVector. Genetic map positions in Morgans.
//' 
//' @return double. The gametic variance.
//' 
//' @details 
//' Variance = Sum(Var_i) + 2 * Sum(Cov_ij)
//' Var_i = a_i^2 / 2  if segregating (P1 != P2)
//' Cov_ij = (1 - 2r_ij) * (a_i * a_j) / 2
//' where r_ij is recombination frequency calculated from Haldane map function:
//' r_ij = 0.5 * (1 - exp(-2 * |d_i - d_j|))
//' 
//' Assumes parents are fully homozygous inbred lines.
//' 
// [[Rcpp::export]]
double calc_cross_variance_cpp(NumericVector p1, NumericVector p2, NumericVector effects, NumericVector map) {
  int m = p1.size();
  if (p2.size() != m || effects.size() != m || map.size() != m) {
    stop("Input vectors must have the same length.");
  }
  
  // Identify segregating loci
  std::vector<int> seg_idx;
  for (int i = 0; i < m; i++) {
    // If parents differ (e.g. 0 vs 2), locus is segregating
    if (std::abs(p1[i] - p2[i]) > 0.5) { 
        seg_idx.push_back(i);
    }
  }
  
  int n_seg = seg_idx.size();
  if (n_seg == 0) return 0.0;
  
  double variance = 0.0;
  
  // Variance Component (Additive Variance)
  // For an F2 population derived from inbred parents (AA x aa), the allele frequency p = 0.5.
  // Additive Variance Var(A) = 2 * p * (1-p) * a^2 = 0.5 * a^2.
  // Here, 'effects' represents the allele substitution effect (a).
  
  for (int idx : seg_idx) {
    double a = effects[idx];
    variance += 0.5 * a * a;
  }
  
  // Covariance Component (Linkage Disequilibrium)
  // The covariance between two loci i and j depends on the linkage phase in the F1 parent.
  // Cov(g_i, g_j) = 2 * Cov(gamete_i, gamete_j)
  // Cov(gamete_i, gamete_j) = +/- 0.25 * (1 - 2r)
  // Thus, Cov(g_i, g_j) = +/- 0.5 * (1 - 2r) * a_i * a_j
  //
  // Phase Determination:
  // - Coupling (+): Parents have same alleles at both loci (e.g., P1 has AABB, P2 has aabb).
  // - Repulsion (-): Parents have mixed alleles (e.g., P1 has AAbb, P2 has aaBB).
  // Derived by checking if P1 has the same genotype code (0 or 2) at both loci.
  
  for (int i = 0; i < n_seg; i++) {
    for (int j = i + 1; j < n_seg; j++) {
      int idx_i = seg_idx[i];
      int idx_j = seg_idx[j];
      
      double d = std::abs(map[idx_i] - map[idx_j]);
      double r = 0.5 * (1.0 - std::exp(-2.0 * d));
      double linkage = (1.0 - 2.0 * r);
      
      // Determine Phase
      // P1[idx_i] is 0 or 2. P1[idx_j] is 0 or 2.
      // If P1 corresponds to 2 at both, coupling.
      // If P1 corresponds to 0 at both, coupling.
      // If P1 is 2 at i and 0 at j, repulsion.
      double phase = (std::abs(p1[idx_i] - p1[idx_j]) < 0.5) ? 1.0 : -1.0;
      
      double a_i = effects[idx_i];
      double a_j = effects[idx_j];
      
      // Add 2 * Cov because matrix is symmetric
      variance += 2.0 * (0.5 * linkage * phase * a_i * a_j);
    }
  }
  
  return variance;
}
