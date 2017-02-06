/* 
 This file copyright (c) 2017 Sean McBane, under the terms of the MIT
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

#include <chrono>
#include <iostream>
#include <random>

using std::chrono::system_clock;

#ifdef __FAST_MATH__
#include <cmath>
#endif /* __FAST_MATH__ */

 /* This struct holds the parameters for a simulated annealing run.
  * Fields are:
  *   + max_temps: The maximum number of temperatures to consider.
  *   + iters_per_temp: Iterations per temperature to compute.
  *   + cost_reduction_tol: A tolerance for breaking out of the annealing run.
  *     If the ratio of cost of the current solution to cost of the initial
  *     solution is less than cost_reduction_tol, the algorithm will return.
  *   + verbose: true for diagnostic output to stdout, false for silent operation.
  */
struct GenSimAnnealParams
{
	unsigned max_temps;
	unsigned iters_per_temp;
	double cost_reduction_tol;
	bool verbose;
};

/*
 * The main simulated annealing routine. The function template has four class
 * parameters:
 *
 * - A class 'solution'. An object of this type must implement 2 public methods:
 *   void solution::perturb(), which modifies the solution slightly in place,
 *   and cost_t solution::cost(), which returns the cost associated with this
 *   particular solution object. This is the value the algorithm attempts to
 *   minimize.
 *
 * - A class 'cost_t'. This is the type returned by solution::cost(). It should
 *   be an integral type convertable to double.
 *
 * - A class acceptance_t. An object of type acceptance_t must be a functor with
 *   signature (cost_t, cost_t, double) -> double. This functor should take as
 *   the first argument the old cost, as a second argument the new (higher) cost,
 *   and as the third argument the current temperature, and return a probability
 *   0 <= p <= 1.0 that the cost increase is accepted.
 *
 * - A class schedule_t. An object of this type should be a functor with
 *   signature (unsigned) -> double; it is called to calculate the temperature
 *   at each *outer* iteration of the algorithm.
 *
 * - A class generator_t. This is an abstract random number generator. It 
 *   should implement operator() to generate a random number >= 0, and max() 
 *   to yield the maximum value of a generated number. Ex.: std::mt19937
 *
 * As a performance consideration, note that one copy of a 'solution' object
 * must be made at each inner iteration. If your solution class has data that
 * isn't changed when perturbing the solution, make sure that copy assignment
 * of the class doesn't copy this data, e.g. by using std::shared_ptr to wrap
 * it.
 */
template<class solution, class cost_t, class acceptance_t, class schedule_t,
	class generator_t>
	solution genericSimulatedAnneal(const solution& s0,
									const GenSimAnnealParams& params,
									const acceptance_t& acceptance_f,
									const schedule_t& cooling_schedule,
									generator_t& generator)
{
	// Initialization.
	cost_t cost_init = s0.cost();
	cost_t cost_old = cost_init;
	// Make a copy of the solution that was passed in:
	solution s = s0;
	solution s1 = s;

	// Keep track of current best solution;
	solution best = s0;
	cost_t best_cost = cost_init;

	// Iterate over temperatures (main iteration).
	for (unsigned outer_iter = 0; outer_iter < params.max_temps; ++outer_iter) {
		double temp = cooling_schedule(outer_iter);
		// Inner iterations at each temperature;
		for (unsigned inner_iter = 0; inner_iter < params.iters_per_temp;
			 ++inner_iter) {

			// Perturb the solution;
			s1.perturb();

			// Whether or not we keep the new solution;
			bool accepted = false;

			// Compute cost of new solution.
			cost_t cost_new = s1.cost();

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

				double x = (double(generator())) / generator.max();
				if (acceptance_prob > x)
					accepted = true;
			}

			// If solution was accepted, update.
			if (accepted) {
				s = s1;
				cost_old = cost_new;

				if (params.verbose)
					std::cout << "Updating solution at outer iteration " << outer_iter
					<< ", inner iteration " << inner_iter << std::endl
					<< "  New cost: " << cost_new << std::endl;

				// Met the cost reduction criterion; break.
				if ((double(cost_new)) / cost_init < params.cost_reduction_tol) {
					if (params.verbose)
						std::cout << "Met cost reduction criterion at "
						<< "outer iteration " << outer_iter << ", inner iteration "
						<< inner_iter << ".\n";

					return s;
				}
			}

			// If it wasn't accepted, we have to copy s back into s1. There's no way
			// around one copy every iteration.
			else {
				s1 = s;
			}
		}
	}

	if (params.verbose)
		std::cout << "Completed " << params.max_temps << " iterations without"
		<< " meeting convergence criterion;\nreturning best value found so far."
		<< std::endl;

	// When returning make sure that we didn't lose a better solution somewhere
	// along the line.
	return s.cost() < best.cost() ? s : best;
}

/*
 * For convenience, this implementation allows the user to ignore specifying
 * a random generator and initializes a Mersenne Twister engine with the 
 * current time by default.
 */
template<class solution, class cost_t, class acceptance_t, class schedule_t>
solution genericSimulatedAnneal(const solution& s0,
								const GenSimAnnealParams& params,
								const acceptance_t& acceptance_f,
								const schedule_t& schedule)
{
	std::mt19937_64 generator(system_clock::to_time_t(system_clock::now()));
	return genericSimulatedAnneal<solution, cost_t, acceptance_t, schedule_t,
		std::mt19937_64>(s0, params, acceptance_f, schedule, generator);
}

#endif /* GENERIC_SIM_ANNEALING_HPP */
