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
  
  // Variance component (diagonal)
  // Var(x) = p(1-p) * (2a)^2? No.
  // Model: Gene count 0,1,2.
  // F2 from AA x aa: 0.25 AA, 0.50 Aa, 0.25 aa.
  // Values: -a, 0, a (if centered) or 0, a, 2a.
  // Var = 0.5 * a^2.
  // Here inputs 'effects' are usually allele substitution effects (alpha).
  // Genetic Variance = 2 * p * (1-p) * alpha^2.
  // For F2, p=0.5. Var = 2 * 0.25 * alpha^2 = 0.5 * alpha^2.
  
  for (int idx : seg_idx) {
    double a = effects[idx];
    variance += 0.5 * a * a;
  }
  
  // Covariance component (off-diagonal)
  // Cov(g_i, g_j) = D_ij * alpha_i * alpha_j * (something)?
  // For F2 from coupling phase (AB/ab): Cov = 0.5 * (1-2r) * a_i * a_j
  // For F2 from repulsion phase (Ab/aB): Cov = -0.5 * (1-2r) * a_i * a_j
  // Need to determine phase.
  // P1: A A  (0 or 2?) Let's say 2 is Alt.
  // If P1[i] == P1[j], they are in coupling (both Ref or both Alt in P1).
  // If P1[i] != P1[j], one is Ref, one Alt -> Repulsion relative to P1?
  // Wait. P1 contributes gamete G1. P2 contributes G2.
  // F1 is G1/G2.
  // Gametic Disequilibrium D_gamete in F1 = 0.25 * (1 - 2r).
  // Genotypic CoVar in F2 = 2 * Cov_gamete = 0.5 * (1 - 2r).
  // Sign depends on linkage phase in F1.
  // If P1 has ++ at i,j and P2 has --, then F1 is ++/-- (Coupling).
  // If P1 has +- and P2 has -+, then F1 is +-/ -+ (Repulsion).
  
  // Logic:
  // Phase = +1 if (P1[i] == P1[j])
  // Phase = -1 if (P1[i] != P1[j])
  // (Assuming P1 and P2 are opposite at these loci, which they must be to segregate)
  
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
