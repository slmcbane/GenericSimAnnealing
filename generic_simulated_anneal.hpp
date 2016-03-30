#ifndef GENERIC_SIM_ANNEALING_HPP
#define GENERIC_SIM_ANNEALING_HPP

#include <cstdlib>

#ifdef __FAST_MATH__
#include <cmath>
#endif /* __FAST_MATH__ */

using std::size_t;

/* This struct holds the results from a run of simulated annealing.
 * Fields are as follows:
 *   + best: The best solution the algorithm found.
 *   + nextbest: The second best solution - I may remove this field if 
 *     it's not useful.
 *   + finalcost: The final value of the cost function (for best solution).
 *   + function_evals: Number of cost function evaluations.
 *   + iterations: Number of iterations performed (where an iteration is 
 *     defined as the period between one temperature change and the next).
*/
template<class solution, typename cost_t>
struct GenSimAnnealResults {
  solution best;
  cost_t finalcost;
  size_t function_evals;
  size_t iterations;
};

/* This struct holds the parameters for a simulated annealing run.
 * Fields are:
 *   + max_temps: The maximum number of temperatures to consider.
 *   + iters_per_temp: Iterations per temperature to compute.
 *   + alpha: The temperature reduction factor ( 0 < a < 1 )
 *   + cost_reduction_tol: A tolerance for breaking out of the annealing run.
 *     If the ratio of cost of the current solution to cost of the initial 
 *     solution is less than cost_reduction_tol, the algorithm will return.
 *   + verbose: true for diagnostic output to stdout, false for silent operation.
 */
struct GenSimAnnealParams {
  size_t max_temps;
  size_t iters_per_temp;
  double alpha;
  double cost_reduction_tol; 
  bool verbose;
};

/* Default values for the simulation parameters for convenience; some or
 * all of these may not be appropriate for a given problem. Use with caution!
 */
const struct GenSimAnnealParams DEFAULT_SIM_ANNEALING_PARAMS = 
  { 100, 500, 0.9, 0.0001, false };

/* The main routine. Templated on a 'solution' class (must be provided by the
 * user) and a typename 'cost_t' which is the variable type returned by a
 * 'solution''s cost() method. This is the value that the algorithm attempts
 * to minimize. The third template parameter is a function which should accept
 * two arguments of type cost_t, with the second greater than the first, and
 * return a probability that the second is selected over the first. The 
 * implementation of this function is left to the user since a suitable 
 * function may vary with the problem, and it did not seem to belong in the
 * 'solution' class.
 *
 * The 'solution' class must expose two methods; one is 'perturb()', which 
 * should return a slightly and randomly modified copy of the solution object.
 * The second is 'cost()' which returns a single value of type cost_t (must be
 * an integral type).
 */
template <class solution, typename cost_t, 
          double (*acceptance_f)(cost_t, cost_t, double)>
GenSimAnnealResults<solution, cost_t> genericSimulatedAnneal
(const solution& s0, const GenSimAnnealParams& params) {

  // Initialization.
  size_t max_temp = params.max_temps;
  size_t iters_per_temp = params.iters_per_temp;
  double alpha = params.alpha;
  double tol   = params.cost_reduction_tol;
  bool verbose = params.verbose;
  size_t fevals = 0;

  // Make a copy of the solution that was passed in:
  solution s = s0;

  // Keep track of current best solution;
  solution best = s0;

  // Main loop.

  // Iterate over temperatures (main iteration).
  for (size_t outer_iter = 0; outer_iter < max_temps; ++outer_iter) {
    // Inner iterations at each temperature;
    for (size_t inner_iter = 0; inner_iter < iters_per_temp; ++inner_iter) {
      // Perturb the solution:
      solution s1 = s.perturb();

      // Whether or not we keep the new solution;
      bool accepted = false;

      // Compute costs of old and new solution.
      cost_t cost_old = s.cost();
      cost_t cost_new = s1.cost();
      fevals += 2;

      // If new cost is less than old, automatically keep this solution.
      // Also update the stored best solution.
      if (cost_new < cost_old) { 
        accepted = true;
        if (best.cost() > cost_new) best = s1;
        ++fevals;
      }
      
      // Otherwise, with some non-zero probability we accept the increase.
      else {
        double acceptance_prob = acceptance_f(cost_old, cost_new, temp);

        // We don't know what happens inside of acceptance_f; must test for
        // NaN (according to IEEE spec comparisons with NaN are always
        // false, but -ffast-math can mess this up and I don't know how the
        // end user will compile this). However, I can check if the macro is
        // defined and only compile the check if fast math is on. gcc only,
        // I don't know if this can be done for other compilers.

#ifdef __FAST_MATH__
        if (std::isnan(acceptance_prob)) acceptance_prob = 0;
#endif /* __FAST_MATH__ */

        double x = (double(rand()))/RAND_MAX;
        if (acceptance_prob > x) accepted = true;
      }

      // If solution was accepted, update.
      if (accepted) {
        s = s1;

        if (verbose) std::cout << "Updating solution at outer iteration "
          << outer_iter << ", inner iteration " << inner_iter << std::endl
          << "  New cost: " << cost_new << "  Old cost: " << cost_old;

        // Met the cost reduction criterion; break.
        if ((double(cost_new)) / cost_init < tol) {

          if (verbose) std::cout << "Met cost reduction criterion at "
           << "outer iteration " << outer_iter << ", inner iteration "
            << inner_iter << ".\n";

          GenSimAnnealingResults res = 
          { s, s.cost(), fevals, outer_iter };
          return res;
        }
      }
    }
    temp = temp*alpha;
  }

  if (verbose) std::cout << "Completed " << max_temps << " iterations without"
    << " meeting convergence criterion;\nreturning best value found so far."
      << std::endl;

  // When returning, check to see if the stored best solution is better than
  // the current one and if it is, return it instead.
  GenSimAnnealingResults res = { s.cost() < best.cost() ? s : best,
    s.cost() < best.cost() ? s.cost() : best.cost(),
    fevals, outer_iter };
  return res;
}


#endif /* GENERIC_SIM_ANNEALING_HPP */
