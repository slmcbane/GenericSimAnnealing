/* This file copyright (c) 2016 Sean McBane, under the terms of the MIT
 License:

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

 */

#ifndef GENERIC_SIM_ANNEALING_HPP
#define GENERIC_SIM_ANNEALING_HPP

#include <cstdlib>
#include <iostream>

#ifdef __FAST_MATH__
#include <cmath>
#endif /* __FAST_MATH__ */

using std::size_t;

/* This struct holds the results from a run of simulated annealing.
 * Fields are as follows:
 *   + best: The best solution the algorithm found.
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
const struct GenSimAnnealParams DEFAULT_SIM_ANNEALING_PARAMS = { 100, 500, 0.9,
    0.0001, false };

/* The main routine. Templated on a 'solution' class (must be provided by the
 * user) and a class 'cost_t' which is the variable type returned by a
 * 'solution''s cost() method. This is the value that the algorithm attempts
 * to minimize.
 *
 * The third template parameter is a class which can be either a function 
 * pointer, a lambda expression or a functor. This object must take two objects
 * of type cost_t (with the second one being greater) and a double value of 
 * temperature and return a value from 0 to 1, which gives the probability that
 * a solution with the second cost is chosen over the first.
 *
 * It is important to note that this is a different approach than many other 
 * implementations, which use a fixed probability function and let the 
 * temperature function be specified on a per problem basis; in this version 
 * the temperature ALWAYS begins at 1.0 and tends to 0.
 *
 * The 'solution' class must expose two methods; one is 'perturb()', which 
 * should slightly and randomly modify the solution in place.
 * The second is 'cost()' which returns a single value of type cost_t (must be
 * an integral type).
 *
 * There are a couple of practical concerns for the user to consider; first,
 * the cost function needs to be relatively cheap to evaluate for large
 * problems as it will be evaluated millions of times. Second, and less
 * obviously, the routine must make a copy of a solution object at every
 * iteration (eithe updating it or reverting the perturbed copy to the 
 * original). This means that if your 'solution' class
 * has any data that gets deep copied during assignment (such as a std::vector) 
 * there will be a large overhead. Use pointers to any data
 * in the object that does not need to change from iteration to iteration
 * to avoid this.
 */
template<class solution, class cost_t, class acceptance_t>
GenSimAnnealResults<solution, cost_t> genericSimulatedAnneal(const solution& s0,
    const GenSimAnnealParams& params, acceptance_t acceptance_f) {

  // Initialization.
  size_t max_temps = params.max_temps;
  size_t iters_per_temp = params.iters_per_temp;
  double alpha = params.alpha;
  double tol = params.cost_reduction_tol;
  bool verbose = params.verbose;
  size_t fevals = 0;
  cost_t cost_init = s0.cost();
  cost_t cost_old = cost_init;

  // Make a copy of the solution that was passed in:
  solution s = s0;
  solution s1 = s;

  // Keep track of current best solution;
  solution best = s0;
  cost_t best_cost = cost_init;

  // Main loop.
  double temp = 1.0;

  // Iterate over temperatures (main iteration).
  for (size_t outer_iter = 0; outer_iter < max_temps; ++outer_iter) {
    // Inner iterations at each temperature;
    for (size_t inner_iter = 0; inner_iter < iters_per_temp; ++inner_iter) {
      // Perturb the solution;
      s1.perturb();

      // Whether or not we keep the new solution;
      bool accepted = false;

      // Compute cost of new solution.
      cost_t cost_new = s1.cost();
      fevals++;

      // If new cost is less than old, automatically keep this solution.
      // Also update the stored best solution.
      if (cost_new < cost_old) {
        accepted = true;
        if (best_cost > cost_new) {
          best = s1;
          best_cost = cost_new;
        }
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

        double x = (double(rand())) / RAND_MAX;
        if (acceptance_prob > x)
          accepted = true;
      }

      // If solution was accepted, update.
      if (accepted) {
        s = s1;
        cost_old = cost_new;

        if (verbose)
          std::cout << "Updating solution at outer iteration " << outer_iter
              << ", inner iteration " << inner_iter << std::endl
              << "  New cost: " << cost_new << std::endl;

        // Met the cost reduction criterion; break.
        if ((double(cost_new)) / cost_init < tol) {
          if (verbose)
            std::cout << "Met cost reduction criterion at "
                << "outer iteration " << outer_iter << ", inner iteration "
                << inner_iter << ".\n";

          GenSimAnnealResults<solution, cost_t> res = { s, s.cost(), fevals,
              outer_iter };
          return res;
        }
      }

      // If it wasn't accepted, we have to copy s back into s1. There's no way
      // around one copy every iteration.
      else {
        s1 = s;
      }
    }
    temp = temp * alpha;
  }

  if (verbose)
    std::cout << "Completed " << max_temps << " iterations without"
        << " meeting convergence criterion;\nreturning best value found so far."
        << std::endl;

  // When returning, check to see if the stored best solution is better than
  // the current one and if it is, return it instead.
  GenSimAnnealResults<solution, cost_t> res = {
      s.cost() < best.cost() ? s : best,
      s.cost() < best.cost() ? s.cost() : best.cost(), fevals, max_temps };
  return res;
}

#endif /* GENERIC_SIM_ANNEALING_HPP */
